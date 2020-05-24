# PyFFAME (Fast and Flexible Allelic MPRA Designer)

A tool for extremely quick and customizable design of multi-allelic MPRA (massively parallel reporter assays)

## Dependencies

PyFFAME was developed for Unix environments. We suggest running python >= 3.7 on a Ubuntu >= 16.04 machine (or similar). 

### Python packages

- pandas >= 0.24.2
- pyfaidx >= 0.5.5.2

Also see the requirements.txt file.

### Further required data

A reference genome in FASTA format. For example, the GRCh37 (hg19) reference genome can be found here [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz].

A list of enzymes and (expanded) restriction sites, the basis for the one provided in this repository can be found here [ftp://ftp.neb.com/pub/rebase/itype2.txt]. These enzymes and restriction sites are based on REBASE data:
> Richard J Roberts, Tamas Vincze, Janos Posfai, Dana Macelis, REBASE--a database for DNA restriction and modification: enzymes, genes and genomes, _Nucleic Acids Res_, Volume 43, 05 November 2014, Pages D298-D299, https://doi.org/10.1093/nar/gku1046

## Installation

Download or clone this repository. 

If you wish to set up your own data base containing rsID variant info, see the alternative pipeline in the 'db_version' folder of this repository. 

## Usage 

Edit the config_file.py file to customize your design. The options are explained here:

| Parameter | Description | Value |
| --- | --- | --- |
| in_vcf | Path to input file in vcf format (columns: CHROM, POS, ID, REF, ALT) | String (file path) or None |
| in_sequence | Path to input file containing additional sequences (columns: ID, FEATURE_SEQ) | String (file path) or None |
| in_barcode | Path to input file (without header) containing sufficient barcodes for the planned design | String (file path) |
| in_barcode_type | File type of barcode input file (json format or plain text with one barcode per line) | 'json' or 'txt' |
| db_genome | Path to reference genome file | String |
| de_order | Design order, given as a string containing the letters 'a' to 'e'; see below for a more in-depth explanation | String |
| de_seq_1 | First added sequence | String |
| de_seq_2 | Second added sequence | String |
| de_seq_3 | Third added sequence | String |
| set_feature_size | Nucleotides to add to each side of the variant (e.g., 85 will yield a total feature size of 171 nt) | Integer |
| set_all_features | Create all features (ref/alt) or only those for the reference alleles? | 1: all; 0: ref .only |
| set_indel_max_length | Maximum length of indels to include (indels exceeding this size will be removed) | Integer |
| set_indel_features | Only create full-length indel features or create additional features accounting for the size difference? See below for a more in-depth explanation | 0: additional; 1: only full-length |
| set_barcodes_per_feature | Number of barcodes to use per feature | Integer |
| set_rev_comp | Create reverse complementary versions of all features? | 1: yes; 0: no |
| enz_file_processed | Path to processed enzyme file (columns: enzymes, site, expanded_sites) | String (file path) or None |
| enz_file | Path to REBASE enzyme file (ignored if processed file used) | String (file path) |
| enz_used | Enzymes used, comma-delimited (e.g. 'EcoRI, SbfI') | String or None |
| enz_sites | Expanded cut sites, python dictionary (e.g., ['GAATTC', 'CCTGCAGG']) | Dictionary as string |
| enz_cumul_cuts | Total restriction sites expected in final feature | Integer |
| enz_cumul_cuts_bc | Total restriction sites expected in (barcode plus sequences right before and after it; used to screen the barcodes for excess restriction sites) | Integer |
| out_format | Output format, json or tab-separated | 'json' or 'tsv' |
| out_output | Path to designed output | String (file path) |

If no VCF info is supplied (only additional sequences), no genome is necessary. 

Note for single nucleotide variants: the pipeline checks, whether the reference allele for each variants is concordant with the base observed at this position in the genome. If not, the reverse complementary base is checked for concordance. If this reversed base is identical to the genomic position, reference and alternative allele are used reverse complemented. Else the original orientation of the alleles is preserved.

### Further information regarding the design order (parameter de_order):
This parameter facilitates flexible combinations of the designed feature. As mentioned above, the parameter expects a string of length one to five, containing the letters a through e. Each letter corresponds to a certain subsequence of the design as per the following:

- a: first added sequence (parameter de_seq_1)
- b: created feature (i.e. the genomic sequence based on the input)
- c: second added sequence (parameter de_seq_2)
- d: barcode (as provided in the in_barcode file parameter)
- e: third added sequence (parameter de_seq_3)

An example (and the config default) would be 'abcde', resulting in the following final feature: (first added sequence)-(feature)-(second added sequence)-(barcode)-(third added sequence)
If subsequences are not required for a design, the corresponding letter can be omitted.

### Further information regarding features created for insertions and deletions (parameter set_indel_features):
When creating features for insertions and deletions, the resulting genomic sequence for either reference or alternative allele will be of a different length than specified in the config file. For example, in the case of a 3-bp deletion, the reference allele sequence might have a length of 171 bp but the resulting alternative allele sequence will have a length of 168 bp. To remedy this, the default behavior of PyFFAME is to create to create the reference allele sequence, the alternative allele sequence of length 168 and a second alternative allele sequence of length 171 (with padding obtained from the genomic context). Setting this parameter to 1 will only create full-length features, i.e. the first of the two alternative allele sequences will be discarded. When all possible features are created, 'SHORT' / 'LONG' is added to the unique feature ID to indicate this.

### Running the pipeline

Call the script using python3 main.py config.py. Note that the config file has to be in the same folder as your terminal and requires a .py file extension. The provided exemplary files can be modified according to the user's design.

The final designed features each contain a unique ID, created in the following way:

[rsID] + [REF / ALT] + [allele number for variants with more than 1 alternative allele] + [SHORT / LONG] + [consecutive number]

Created output:

File containing the designed features; file containing features which were removed due to additional enzyme restriction sites; log file. 

### Example

Example files (rsID input (for the DB-version), VCF input, additional sequences input) are contained in the example_files folder. The config.py contained in this repository can be used with these files. 

## Barcodes

Barcodes can be supplied by the user. Alternatively, the companion tool PyBarcodes (https://github.com/DNAbased/PyBarcodes) can be used to generate barcodes based on custom parameters.

The following sets of barcodes are already provided in json format:

| Length [nt] | Number of barcodes |
| --- | --- |
| 6 | 36 |
| 7 | 127 |
| 8 | 351 |
| 9 | 1,215 |
| 10 | 4,217 |
| 11 | 10,807 |
| 12 | 38,718 |
| 13 | 99,853 |
| 14 | 370,911 |
| 15 | 954,072 |

The following criteria were applied to obtain these barcodes:

- Removed barcodes with repeats of more than two identical nucleotides
- Removed barcodes without at least one of each nucleotide
- Lower GC threshold: 40 %
- Upper GC threshold: 60 %
- Removed barcodes containing ACA/CAC/GTG/TGT
- Minimum levenshtein distance between barcodes of 3

## Benchmarking

We benchmarked PyFFAME and compared it to two other tools for the design of MPRA assays, MPRAnator and MPRA Design Tools:
> Ilias Georgakopoulos-Soares, Naman Jain, Jesse M Gray, Martin Hemberg, MPRAnator: a web-based tool for the design of massively parallel reporter assay experiments, _Bioinformatics_, Volume 33, Issue 1, 01 January 2017, Pages 137–138, https://doi.org/10.1093/bioinformatics/btw584

> Andrew R Ghazi, Edward S Chen, David M Henke, Namrata Madan, Leonard C Edelstein, Chad A Shaw, Design tools for MPRA experiments, _Bioinformatics_, Volume 34, Issue 15, 01 August 2018, Pages 2682–2683, https://doi.org/10.1093/bioinformatics/bty150

The relevant information can be found in the associated repository (https://github.com/DNAbased/PyFFAME_benchmark).

## Testing (WIP)

```
python3.7 -m unittest
```
