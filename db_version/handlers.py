import pandas as pd
import re
import json
from pymongo import MongoClient
from pyfaidx import Fasta


# ToDo: Quality Check


class BarcodeHandler:
    def __init__(self):
        pass

    @staticmethod
    def add_shuffled_bc(feature_dict, bc_list, order, five_seq, spacer_seq, three_seq, n=100):
        """add barcodes and additional sequences to design final features"""
        barcodes = pd.DataFrame(bc_list).sample(frac=1).reset_index(drop=True)  # shuffle barcodes and reset index
        barcode_iterator = 0
        dict_out = dict()
        req_zeroes = '{:0' + str(len(str(n))) + '}'  # replaces '"{:03}"' in the number formatting below
        order = order + 'fffff'
        f = 'f'
        a = five_seq
        c = spacer_seq
        e = three_seq
        for key, value in feature_dict.items():
            for i in range(n):
                b = value
                d = barcodes.iat[barcode_iterator, 0]
                dict_out[key + '_' + str(req_zeroes.format(i))] = '{0}{1}{2}{3}{4}'.format(eval(order[0]),
                                                                                           eval(order[1]),
                                                                                           eval(order[2]),
                                                                                           eval(order[3]),
                                                                                           eval(order[4]))
                barcode_iterator += 1
        return dict_out

    @staticmethod
    def create_intermediate_feature(df_in, order, seq_1, seq_2, seq_3, bc):
        """create intermediate feature for restriction check"""
        order = order + 'fffff'  # to prevent out of index when users use e.g. 'abcd' as order
        f = ''  # adds empty string; see line above
        a = seq_1
        c = seq_2
        d = bc
        e = seq_3
        for i in df_in.index:
            b = df_in.at[i, 'FEATURE_SEQ']
            df_in.at[i, 'CHECK_SEQ'] = '{0}{1}{2}{3}{4}'.format(eval(order[0]), eval(order[1]), eval(order[2]),
                                                                eval(order[3]), eval(order[4]))
        return df_in


# ToDo: Quality Check + combine


class EnzymeHandler:
    def __init__(self):
        pass

    @staticmethod
    def cut_site(enzyme_name, enzyme_df):  # not used
        """obtain enzyme-specific restriction sites"""
        try:
            output = enzyme_df.loc[enzyme_df['enzymes'] == enzyme_name]['sites']  # the cut site
            output = output.iloc[0]  # transform cut site to string
        except KeyError:
            raise KeyError("expected column names ['enzymes', 'sites']")
        return output

    @staticmethod
    def expanded_cut_site(enzyme_name, enzyme_df):  # not used
        """obtain enzyme-specific restriction sites (EXPANDED)"""
        try:
            output = enzyme_df.loc[enzyme_df['enzymes'] == enzyme_name]['expanded_sites']  # get expanded cut site
            output = output.iloc[0]  # transform expanded cut site to string
        except KeyError:
            raise KeyError("expected column names ['enzymes', 'expanded_sites']")
        return output

    @staticmethod
    def expanded_cut_site_multi(enzymes_list, enzyme_df):
        """obtain a list of multiple expanded restriction sites"""
        output = []
        try:
            for i in enzymes_list:
                current = enzyme_df.loc[enzyme_df['enzymes'] == i]['expanded_sites']
                current = current.iloc[0]
                output.append(current)
        except KeyError:
            raise KeyError("expected column names ['enzymes', 'expanded_sites']")
        return output

    @staticmethod
    def cut_check(sequence, expanded_cut_site):  # not used
        """check a string for an expanded restriction site and return the number of cut sites"""
        return len(re.findall(expanded_cut_site, sequence))

    @staticmethod
    def cut_check_multi(sequence, sites_list):
        """check a string for multiple restriction enzymes using a list of expanded cut sites"""
        total_cuts = 0
        for i in sites_list:
            total_cuts += len(re.findall(i, sequence))
        return total_cuts

    def barcode_check(self, bc, order, seq_1, seq_2, seq_3, sites, cuts):  # not a good name
        """check barcodes for restriction sites"""
        bc_okay = []
        bc_not_okay = []
        order = order + 'fffff'
        order = order[
                order.find('d') - 1:order.find('d') + 2]  # only take whatever is next to the barcode for the check
        f = ''
        a = seq_1
        c = seq_2
        e = seq_3
        b = ''  # sequence of interest should not be next to barcode + not used
        before = eval(order[0])
        after = eval(order[2])
        whole = list(before + item + after for item in bc)  # list of generator expression
        for i in range(len(whole)):
            d = bc[i]
            hits = self.cut_check_multi(whole[i], sites)
            if hits <= cuts:
                bc_okay.append(d)
            else:
                bc_not_okay.append(d)
        return bc_okay, bc_not_okay

    @staticmethod
    def feature_check(feature_df, cut_sites, cut_count):
        """check created features for restriction sites"""
        for i in feature_df.index:
            feature_df.at[i, 'CUTS'] = 0
            for j in cut_sites:
                feature_df.at[i, 'CUTS'] += len([match for match in re.finditer(j, feature_df.at[i, 'CHECK_SEQ'])])
        failed_IDs = feature_df[feature_df.CUTS > cut_count].original_ID
        cleaned_df = feature_df[~feature_df.original_ID.isin(failed_IDs)]
        removed_df = feature_df[feature_df.original_ID.isin(failed_IDs)]
        # cleaned_df = feature_df[feature_df.CUTS <= cut_count].reset_index()
        # removed_df = feature_df[feature_df.CUTS > cut_count].reset_index()
        return cleaned_df, removed_df

    @staticmethod
    def cut_position(sequence, expanded_cut_site):  # Wrapper around builtin function?
        """get position at which an enzyme cuts"""
        enzyme_hits = re.search(expanded_cut_site, sequence)
        return enzyme_hits.span()

    @staticmethod
    def create_enzyme_list(enzyme_list):
        """tvs list of type II enzymes from 'http://rebase.neb.com/rebase/link_itype2' for python-usable list
            output column names: ['enzymes', 'sites', 'expanded_sites']"""
        with open(enzyme_list, 'r') as f:
            content = f.read()
        content = content.split('\n')  # split into lines
        content = pd.Series(content)  # create series
        content = content.str.split('\t', expand=True)  # split by tabs
        content = pd.DataFrame(content)  # create df
        content = content.drop([1, 3, 4, 5], axis=1)  # drop irrelevant columns
        content = content.ix[content[2].notnull()]  # drop empty rows
        content = content.reset_index(drop=True)  # reset index (because of dropped rows)
        content = content.rename(columns={0: 'enzymes', 2: 'sites'})  # rename columns
        for i in range(len(content.index)):
            content.at[i, 'sites'] = re.sub('\^', '', content.at[i, 'sites'])  # remove ^s from cut site
            content.at[i, 'sites'] = re.sub('\(-?[0-9]+/-?[0-9]+\)', '',
                                            content.at[i, 'sites'])  # remove (int/int) stuff
            content.at[i, 'sites'] = re.sub(',.*', '', content.at[i, 'sites'])  # only keep the first recognition site
        bases = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'W': '[AT]', 'S': '[CG]', 'M': '[AC]', 'K': '[GT]',
                 'R': '[AG]', 'Y': '[CT]', 'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
                 'N': '[ACGT]'}  # ambiguous nucleotides
        for i in range(len(content.index)):
            content.at[i, 'expanded_sites'] = ''.join(
                bases.get(j) for j in content.at[i, 'sites'])  # create expanded sites using bases dict
        return content


# ToDo: Quality Check


class GenomicHandler:
    def __init__(self):
        pass

    @staticmethod
    def replace_chrom(variant_dataframe, column='CHROM'):
        return variant_dataframe.replace({column: {'X': 23, 'Y': 24, 'MT': 25}})

    @staticmethod
    def replacement0(string, pos, char):
        """replace a char in a string defined by its position (e.g. change REF to ALT allele in a sequence) 0-based"""
        string = string[0:pos] + string[pos:(pos + 1)].replace(string[pos], char) + string[(pos + 1):]
        return string

    @staticmethod
    def replacement(string, pos, char):
        """replace a char in a string defined by its position (e.g. change REF to ALT allele in a sequence) 1-based!"""
        pos -= 1
        string = string[0:pos] + string[pos:(pos + 1)].replace(string[pos], char) + string[(pos + 1):]
        return string

    @staticmethod
    def get_genomic_context(df, genome, n=85):
        """get the genomic context of variants
           necessary input columns: ['CHROM', 'POS']; added output columns: ['SEQ', 'SEQ_LONG']"""
        for i in df.index:
            df.at[i, 'SEQ_LONG'] = Fasta(genome)[str(df.at[i, 'CHROM'])][
                                   (df.at[i, 'POS'] - (n + 101)):(df.at[i, 'POS'] + (n + 100))].seq
            df.at[i, 'SEQ'] = df.at[i, 'SEQ_LONG'][100:-100]
        return df

    def ref_alt(self, df, indels=0, length=85):
        """new version to duplicate rows (basically contains rev check)
            duplicate rows, altering SNP_ID to SNP_ID_ref and SNP_ID_alt
            usage: 'variant_data_frame = ref_alt_rows(variant_data_frame, orientation, indels)'
            indels
            (0) create all possible features, for insertions: REF/ALT/REF_SHORT, for deletions: REF/ALT/ALT_SHORT
            (1) only create full-length features, i.e. REF/ALT
            note: additional sequences, e.g. for the creation of a full-length ALT version will be taken from the 3' end of the seq"""
        df2 = pd.DataFrame()
        for i in df.index:
            c_A = df[i:(i + 1)].copy()  # c = current (wrt to for loop)
            c_B = df[i:(i + 1)].copy()
            c_C = df[i:(i + 1)].copy()
            c_A.at[i, 'ID'] = '_'.join(df.at[i, 'ID'].split('_')[:-1]) + '_REF'
            c_B.at[i, 'ID'] = '_'.join(df.at[i, 'ID'].split('_')[:-1]) + '_ALT_' + df.at[i, 'ID'].split('_')[-1]
            if df.at[i, 'INDEL'] == 0:
                if df.at[i, 'REF'] == df.at[i, 'SEQ'][85]:  # ref allele (and, supposedly alt allele) is forward correct
                    c_A.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ']
                    c_B.at[i, 'FEATURE_SEQ'] = self.replacement(df.at[i, 'SEQ'], (length + 1), df.at[i, 'ALT'])
                elif DataHelper.revcomp(df.at[i, 'REF']) == df.at[i, 'SEQ'][85]:  # revcomp correct
                    c_A.at[i, 'FEATURE_SEQ'] = DataHelper.revcomp(df.at[i, 'SEQ'])
                    c_B.at[i, 'FEATURE_SEQ'] = self.replacement(DataHelper.revcomp(df.at[i, 'SEQ']), (length + 1), df.at[i, 'ALT'])
                else:  # neither; assume and use forward
                    c_A.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ']
                    c_B.at[i, 'FEATURE_SEQ'] = self.replacement(df.at[i, 'SEQ'], (length + 1), df.at[i, 'ALT'])
                df2 = df2.append([c_A, c_B], ignore_index=True)
            else:
                c_A.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ']
                if indels == 0:
                    if df.at[i, 'INDEL_TYPE'] == 'insertion':
                        c_B.at[i, 'FEATURE_SEQ'] = self.replacement(df.at[i, 'SEQ'][:-(len(df.at[i, 'ALT']) - 1)],
                                                               (length + 1), df.at[i, 'ALT'])
                        c_C.at[i, 'ID'] = '_'.join(df.at[i, 'ID'].split('_')[:-1]) + '_REF_SHORT'
                        c_C.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ'][:-(len(df.at[i, 'ALT']) - 1)]
                    elif df.at[i, 'INDEL_TYPE'] == 'deletion':
                        c_B.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ'][:length] + df.at[i, 'ALT'] + df.at[i, 'SEQ_LONG'][(100 + length + len(df.at[i, 'REF'])):((100 + length * 2) + len(df.at[i, 'REF']))]
                        c_C.at[i, 'ID'] = '_'.join(df.at[i, 'ID'].split('_')[:-1]) + '_ALT_SHORT_' + \
                                          df.at[i, 'ID'].split('_')[-1]
                        c_C.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ'][:length] + df.at[i, 'ALT'] + df.at[i, 'SEQ'][(length + len(df.at[i, 'REF'])):]
                    df2 = df2.append([c_A, c_B, c_C], ignore_index=True)
                elif indels == 1:
                    if df.at[i, 'INDEL_TYPE'] == 'insertion':
                        c_B.at[i, 'FEATURE_SEQ'] = self.replacement(df.at[i, 'SEQ'][:-(len(df.at[i, 'ALT']) - 1)],
                                                               (length + 1), df.at[i, 'ALT'])
                    elif df.at[i, 'INDEL_TYPE'] == 'deletion':
                        c_B.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ'][:length] + df.at[i, 'ALT'] + df.at[i, 'SEQ_LONG'][(100 + length + len(df.at[i, 'REF'])):((100 + length * 2) + len(df.at[i, 'REF']))]
                    df2 = df2.append([c_A, c_B], ignore_index=True)
        df2 = df2.drop_duplicates(subset='ID')
        df2 = df2.reset_index(drop=True)
        return df2

    # TODO: could be reworked to make it faster
    @staticmethod
    def ref_only(df, length=85):
        """as ref_alt but only creates reference allele features """
        df2 = pd.DataFrame()
        for i in df.index:
            c_A = df[i:(i + 1)].copy()
            c_A.at[i, 'ID'] = '_'.join(df.at[i, 'ID'].split('_')[:-1]) + '_REF'
            if df.at[i, 'INDEL'] == 0:
                if df.at[i, 'REF'] == df.at[i, 'SEQ'][length]:  # ref allele (and, supposedly alt allele) is forward correct
                    c_A.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ']
                elif DataHelper.revcomp(df.at[i, 'REF']) == df.at[i, 'SEQ'][length]:  # revcomp correct
                    c_A.at[i, 'FEATURE_SEQ'] = DataHelper.revcomp(df.at[i, 'SEQ'])
                else:  # neither --> assume and use forward
                    c_A.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ']
            else:
                c_A.at[i, 'FEATURE_SEQ'] = df.at[i, 'SEQ']
            df2 = df2.append([c_A], ignore_index=True)
        df2 = df2.drop_duplicates(subset='ID')
        df2 = df2.reset_index(drop=True)
        return df2

    @staticmethod
    def revcomp_features(df):
        """add revcomp versions of all allelic features"""
        df2 = pd.DataFrame()
        for i in df.index:
            c_A = df[i:(i + 1)].copy()
            c_B = df[i:(i + 1)].copy()
            c_B.at[i, 'ID'] += '_REV'
            c_B.at[i, 'FEATURE_SEQ'] = DataHelper.revcomp(c_B.at[i, 'FEATURE_SEQ'])
            df2 = df2.append([c_A, c_B], ignore_index=True)
        df2.reset_index(drop=True)
        return df2

    @staticmethod
    def indel_check(df):
        """check if a variant is an indel (based on the length of the allele)
        returns the dataframe, number of insertionas and number of deletions"""
        deletion_count, insertion_count = 0, 0
        try:
            for i in range(len(df.index)):
                df.at[i, 'INDEL'] = 0
                df.at[i, 'INDEL_TYPE'] = 'NA'
                ref_len = len(df.at[i, 'REF'])
                alt_len = len(df.at[i, 'ALT'])
                if ref_len != alt_len:  # ref and alt lengths differ
                    df.at[i, 'INDEL'] = 1
                    if ref_len > alt_len:  # ref longer: deletion
                        df.at[i, 'INDEL_TYPE'] = 'deletion'
                        deletion_count += 1
                    elif ref_len < alt_len:  # alt longer: insertion
                        df.at[i, 'INDEL_TYPE'] = 'insertion'
                        insertion_count += 1
        except KeyError:
            raise KeyError("Expected Input Colums ['REF', 'ALT'], added output columns ['INDEL', 'INDEL_TYPE']")
        df = df.reset_index(drop=True)  # reset index
        return df, insertion_count, deletion_count

    @staticmethod
    def remove_indels(df, size=10):
        """remove indels exceeding a certain size limit"""
        try:
            df = df[df['REF'].apply(lambda x: len(x) <= size) & df['ALT'].apply(lambda x: len(x) <= size)]
            df = df.reset_index(drop=True)  # reset index
        except KeyError:
            raise KeyError("Expected Input Columns ['REF', 'ALT']")
        return df


class DataHelper:
    def __init__(self):
        pass

    @staticmethod
    def create_seq_dict(df):
        """create a dictionary containing the ID and the sequence"""
        seq_dict = dict()
        for i in df.index:
            seq_dict[df.at[i, 'ID']] = df.at[i, 'FEATURE_SEQ']  # convert every line into a single dict entry
        return seq_dict

    @staticmethod
    def json_dump(file_name, dictionary):
        with open(file_name, 'w') as f:
            json.dump(dictionary, f, indent=1)

    @staticmethod
    def convert_df_to_list(df, target_col, sep, id_col):
        """split variants containing multiple aternative alleles into one row per alternative allele;
            returns the number of bi-/tri-/quad-allelic and unknown variants"""
        row_accumulator = []
        bi_all, tri_all, quad_all, unknown = 0, 0, 0, 0

        def list_to_rows(row, sep):
            iterator = 1
            split_row = row[target_col].split(sep)
            for s in split_row:
                new_row = row.to_dict()
                new_row[target_col] = s
                new_row[id_col] = new_row[id_col] + '_' + str(iterator)
                row_accumulator.append(new_row)
                iterator += 1

        df.apply(list_to_rows, axis=1, args=(sep,))
        split_df = pd.DataFrame(row_accumulator)
        for i in df.index:
            separator_count = df[target_col][i].count(sep)
            if separator_count == 0:
                bi_all += 1
            elif separator_count == 1:
                tri_all += 1
            elif separator_count == 2:
                quad_all += 1
            else:
                unknown += 1
        return split_df, bi_all, tri_all, quad_all, unknown

    @staticmethod
    def revcomp(seq):
        """reverse complement DNA sequence"""
        comp = str.maketrans('ATCG', 'TAGC')
        return seq.upper().translate(comp)[::-1]  # translate from last to first character

    @staticmethod
    def rev(seq):
        """reverse DNA sequence"""
        return seq.upper()[::-1]  # return from last to first character

    @staticmethod
    def comp(seq):
        """complement DNA sequence"""
        comp = str.maketrans('ATCG', 'TAGC')
        return seq.upper().translate(comp)  # simply translate

    @staticmethod
    def read_list(filename):
        """Read a list of rsIDS"""
        with open(filename, 'r') as f:
            content = f.read().splitlines()  # read in file and split into lines
        return content  # add '.strip()'?


class MongoHandler:
    def __init__(self, username, password, host='localhost', port=27017, timeout=2500):
        self.username = username
        self.password = password
        self.client = MongoClient('mongodb://' + self.username + ':' + self.password + '@' + host + ':' + str(port), serverSelectionTimeoutMS=timeout)

    def get_variants(self, id_list, database='dbsnp_b151', collection='rs'):
        db = self.client[database]  # set database
        coll = db[collection]  # set collection
        df = pd.DataFrame(list((coll.find({'ID': {'$in': id_list}}, {'_id': False}))))  # search for all ids at once
        return df

    def get_single_variant(self, rsID, database='dbsnp_b151', collection='rs'):
        db = self.client[database]
        coll = db[collection]
        return pd.DataFrame(list(coll.find({'ID': rsID}, {'_id': False})))  # return output of single search

    def got_merged(self, id_list, database='dbsnp_b151', collection='merged'):
        db = self.client[database]
        coll = db[collection]
        id_list_mod = [int(i[2:]) for i in id_list]  # create ids suitable for the merged coll
        df = pd.DataFrame(list((coll.find({'OLD': {'$in': id_list_mod}}, {'_id': False}))))
        return df

    def get_single_merged(self, rsID, database='dbsnp_b151', collection='merged'):
        """check a singe variant for 'merged' status"""
        db = self.client[database]
        coll = db[collection]
        simple_ID = int(rsID[2:])
        merged_df = pd.DataFrame(list(coll.find({'OLD': simple_ID}, {'_id': False})))
        return merged_df
