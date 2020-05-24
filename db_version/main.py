def main():
    import pandas as pd
    import sys
    import importlib
    import logging
    from handlers import BarcodeHandler, EnzymeHandler, GenomicHandler, DataHelper, MongoHandler

    config_file = importlib.import_module(sys.argv[1][0:-3])

    # set up logger
    logger = logging.getLogger(__name__)
    log_handler = logging.FileHandler(config_file.out_output + '_logfile.txt')
    logging.basicConfig(level=logging.DEBUG, format='%(message)s') # everything above debug
    log_format = logging.Formatter('%(asctime)s.%(msecs)03d - %(levelname)s: %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    log_handler.setFormatter(log_format)
    logger.addHandler(log_handler)

    # load hadlers
    bh = BarcodeHandler()
    eh = EnzymeHandler()
    gh = GenomicHandler()
    dh = DataHelper()

    logger.info('---MPRA design initiated---')
    logger.info('Used config file: ' + str(sys.argv[1]))

    # read in rs input file and obtain respective info
    if config_file.in_variant is not None:
        mongodb_auth = importlib.import_module(config_file.db_auth[0:-3])
        mh = MongoHandler(username=mongodb_auth.username, password=mongodb_auth.password)
        logger.info('-Reading in rsID input file-')
        rs_variants = dh.read_list(config_file.in_variant)
        logger.info('Total rsID variants: ' + f'{len(rs_variants):,}')
        logger.info('-Getting rsID variant info-')
        rs_df = mh.get_variants(rs_variants, database=config_file.db_database, collection=config_file.db_collection_rs)
        logger.info('Obtained rsID variant info for ' + f'{len(rs_df):,}' + ' variants')
        # for now: drop 'RV' column (if it exists; might be excluded while setting up the database)
        if 'RV' in rs_df.columns:
            rs_df = rs_df.drop(labels='RV', axis=1)
    else:
        rs_df = pd.DataFrame()

    # read in vcf input file
    if config_file.in_vcf is not None:
        logger.info('-Reading in VCF-like input file-')
        vcf_variants = pd.read_csv(config_file.in_vcf, sep='\t')
        logger.info('Total VCF-like variants: ' + f'{len(vcf_variants):,}')
        # join rs and vcf for multi check
        vcf_rs_df = pd.concat([rs_df, vcf_variants], sort=False, ignore_index=True)
    else:
        vcf_rs_df = rs_df

    # multi allelic check
    logger.info('-Splitting multi-allelic variants-')
    split_df, bi_count, tri_count, quad_count, unknown = dh.convert_df_to_list(vcf_rs_df, target_col='ALT', sep=',', id_col='ID')
    logger.info('Number of bi-allelic variants: ' + f'{bi_count:,}')
    logger.info('Number of tri-allelic variants: ' + f'{tri_count:,}')
    logger.info('Number of quad-allelic variants: ' + f'{quad_count:,}')
    if unknown > 0:
        logger.info('Number of variants with more than 4 alleles (likely indels): ' + f'{unknown:,}') # but all of them are designed (!) # e.g. rs61087238

    # get relevant enzyme information (either by checking the enzyme cut sites or by using the provided ones)
    logger.info('-Getting enzyme information-')
    if config_file.enz_used is not None:
        if config_file.enz_file_processed is not None:
            enzyme_df = pd.read_csv(config_file.enz_file_processed)
        else:
            enzyme_df = eh.create_enzyme_list(config_file.enz_file)
        enzyme_cut_sites = eh.expanded_cut_site_multi(config_file.enz_used.split(','), enzyme_df)
    else:
        enzyme_cut_sites = config_file.enz_sites

    # read in barcodes
    logger.info('-Reading in barcodes-')
    if config_file.in_barcode_type == 'json':
        import json
        with open(config_file.in_barcode, 'r') as my_file:
            barcodes = json.load(my_file)
        if type(barcodes) == dict: # the pre-generated barcodes are in json files created from dictionaries, therefore they are not read in as lists automatically; however, this is required for downstream barcode stuff
            barcodes = list(barcodes.values())[0]
    else:
        barcodes = dh.read_list(config_file.in_barcode)
    logger.info('Total barcodes: ' + f'{len(barcodes):,}')

    # check barcodes for restriction sites, remove failed barcodes # log: number of removed barcodes; number of remaining barcodes
    logger.info('-Checking barcodes for restriction sites-')
    bc_use, bc_discard = eh.barcode_check(bc=barcodes, order=config_file.de_order, seq_1=config_file.de_seq_1, seq_2=config_file.de_seq_2, seq_3=config_file.de_seq_3, sites=enzyme_cut_sites, cuts=config_file.enz_cumul_cuts_bc)
    logger.info('Suitable barcodes: ' + f'{len(bc_use):,}')
    logger.info('Discarded barcodes: ' + f'{len(bc_discard):,}')

    # read in additional sequencs (e.g. controls)
    if config_file.in_sequence is not None:
        logger.info('-Reading in additional sequences-')
        add_seqs = pd.read_csv(config_file.in_sequence, sep="\t")
        logger.info('Total additional sequences: ' + str(len(add_seqs)))
    else:
        add_seqs = [] # because the length of this variable is used in the calculation below

    # check indel status
    logger.info('-Checking indel status-')
    indel_df, del_count, in_count = gh.indel_check(split_df)
    logger.info('Deletions: ' + f'{int(del_count):,}')
    logger.info('Insertions: ' + f'{int(in_count):,}')

    # remove indels exceeding the maximum size
    logger.info('-Removing indels-')
    indels_before = len(indel_df)
    indel_df = gh.remove_indels(indel_df, size=config_file.set_indel_max_length)
    indels_after = len(indel_df)
    indels_removed = indels_before - indels_after
    if indels_removed > 0:
        logger.info('Total indels removed: ' + f'{indels_removed:,}')

    # check necessary and available number of barcodes
    logger.info('-Checking number of barcodes-')
    bc_needed = (len(indel_df) * (config_file.set_all_features + 1) * (config_file.set_rev_comp + 1) + len(add_seqs)) * config_file.set_barcodes_per_feature
    logger.info('Barcodes required: ' + f'{int(bc_needed):,}')
    if bc_needed > len(bc_use):
        logger.warning('Number of barcodes is insufficient')
        logger.warning('---Exiting---')
        import sys
        sys.exit()
    else:
        logger.info('Number of barcodes is sufficient')

    # get genomic context
    logger.info('-Getting genomic context-')
    genomic_df = gh.get_genomic_context(indel_df, genome=config_file.db_genome, n=config_file.set_feature_size)

    # duplicate ID column for downstream filtering of features failing the restriction check
    genomic_df['original_ID'] = genomic_df['ID'].copy()

    # create features
    if config_file.set_all_features == 1:
        logger.info('-Creating allelic features-')
        feature_df = gh.ref_alt(genomic_df, indels=config_file.set_indel_features, length=config_file.set_feature_size)
    else:
        logger.info('-Creating reference features-')
        feature_df = gh.ref_only(genomic_df, length=config_file.set_feature_size)

    # create reverse complementary features if necessary
    if config_file.set_rev_comp != 0:
        logger.info('-Creating reverse complementary features')
        feature_df = gh.revcomp_features(feature_df)

    # merge genomic_df with additional seqs (if it is present)
    if config_file.in_sequence is not None:
        logger.info('-Merging sequences-')
        all_seqs_df = pd.concat([feature_df, add_seqs], sort=False, ignore_index=True)
    else:
        all_seqs_df = feature_df

    # create intermediate features for restriction check # barcode should not be next to seq because only the first barcode is used for this test(!)
    logger.info('-Checking restriction sites-')
    all_seqs_df = bh.create_intermediate_feature(all_seqs_df, order=config_file.de_order, seq_1=config_file.de_seq_1, seq_2=config_file.de_seq_2, seq_3=config_file.de_seq_3, bc=bc_use[0])

    # check restriction sites and report number of failed features
    cleaned_df, removed_df = eh.feature_check(all_seqs_df, cut_sites=enzyme_cut_sites, cut_count=config_file.enz_cumul_cuts)
    removed_count = len(removed_df)
    if removed_count > 0:
        logger.warning('Features removed due to too many restriction sites: ' + f'{removed_count:,}')

    # Converting to dictionary
    logger.info('-Converting to dictionary-')
    feature_dict = dh.create_seq_dict(cleaned_df)

    # add (shuffled) barcodes and additional sequences
    logger.info('-Creating final features-')
    final_dict = BarcodeHandler().add_shuffled_bc(feature_dict=feature_dict, bc_list=bc_use, order=config_file.de_order,
                                                  five_seq=config_file.de_seq_1, spacer_seq=config_file.de_seq_2,
                                                  three_seq=config_file.de_seq_3, n=config_file.set_barcodes_per_feature)
    logger.info('Total number of created features: ' + f'{len(final_dict):,}')

    # dump (json or tsv)
    if config_file.out_format == 'json':
        logger.info('-Writing final output (JSON)-')
        dh.json_dump(config_file.out_output + '.json', final_dict)
    else:
        logger.info('-Writing final output (TSV)-')
        pd.DataFrame.from_dict(final_dict).to_csv(config_file.out_output + '.tsv', sep='\t')
    removed_df.to_csv(config_file.out_output + '_removed_features.tsv', sep='\t', index=False)

    logger.info('---MPRA design finished---')


if __name__ == '__main__':
    main()
