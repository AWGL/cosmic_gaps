import pandas as pd
import csv
import re
import yaml
import argparse


def get_args():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--cosmic_file', action='store', required=True, 
        help='Filepath to input COSMIC file. REQUIRED.'
    )
    parser.add_argument(
        '--config', action='store', required=True, 
        help='Filepath to config file. REQUIRED.'
    )

    return parser.parse_args()


def prefilter_cosmic_data(cosmic_file, config):
    """
    file is too large to go straight into pandas dataframe, prefilter for genes of interest first
    """
    print('Pre-filtering COSMIC file for genes of interest')

    # get a list of genes from all referral types
    all_genes = []
    for referral in config.values():
        all_genes += referral['gene_list']

    # do a crude filter for any lines containing any of the genes (essentially the same as running grep)
    filtered_dict = []
    with open(cosmic_file, 'r') as f:

        for n, line in enumerate(f):
            # take first line as header, reformat to be all caps and only underscores for seperating words
            if n == 0:
                header = [ l.upper().replace(' ', '_').replace('-', '_') for l in line.rstrip().split('\t') ]

            # grep for gene list, save any hits to list
            regex = '|'.join(all_genes)
            if re.search(regex, line):
                filtered_dict.append(line.rstrip().split('\t'))

    # convert list to pandas dataframe
    df = pd.DataFrame.from_records(filtered_dict, columns = header)

    # handy print statements if you need to add a new config
    #print(df.head())
    #print(df.PRIMARY_SITE.unique())

    return df


def filter_by_referral(df, config):
    """
    loop through all referrals and filter for genes/ sites of interest
    """
    for referral in config:

        print(f'Filtering {referral}')

        # extract config for referral type
        gene_list = config[referral]['gene_list']
        site_list = config[referral]['primary_site_list']

        # filter dataframe and add column concatenating gene name and HGVS

        filtered_df = df[ (df.GENE_NAME.isin(gene_list)) & (df.PRIMARY_SITE.isin(site_list)) ]

        filtered_df['GENE_AND_HGVS'] = filtered_df.apply(lambda x: concat_gene(x), axis=1)

        # count based on concatenated gene and HGVS column, turn outputed series back into a dataframe
        mutations_counts = filtered_df.GENE_AND_HGVS.value_counts().to_frame()
        mutations_counts.columns = ['COUNT']
        mutations_counts['DESC'] = mutations_counts.index

        # add gene name and HGVS to new dataframe
        mutations_counts['GENE_NAME'] = mutations_counts.apply(lambda x: extract_gene(x)[0], axis=1)
        mutations_counts['HGVS_C'] = mutations_counts.apply(lambda x: extract_gene(x)[1], axis=1)
        mutations_counts['HGVS_P'] = mutations_counts.apply(lambda x: extract_gene(x)[2], axis=1)
        mutations_counts['HGVS_G'] = mutations_counts.apply(lambda x: extract_gene(x)[3], axis=1)

        # split genomic coords from HGVS_G
        mutations_counts['CHR'] = mutations_counts.apply(lambda x: match_genomic(x)[0], axis=1)
        mutations_counts['START'] = mutations_counts.apply(lambda x: match_genomic(x)[1], axis=1)
        mutations_counts['END'] = mutations_counts.apply(lambda x: match_genomic(x)[2], axis=1)

        # save output to CSV
        mutations_counts.to_csv(
            f'{referral}.csv',
            index=False,
            columns=['CHR', 'START', 'END', 'GENE_NAME', 'HGVS_C', 'HGVS_P', 'HGVS_G', 'COUNT']
        )




def concat_gene(x):
    """
    concatenate gene name and all HGVS records for each row of the dataframe
    """
    return f"{x['GENE_NAME']};{x['HGVSC']};{x['HGVSP']};{x['HGVSG']}"


def extract_gene(x):
    """
    split the concatenated gene name/ HGVS field back into its parts
    change any empty fields to None
    """
    l = x['DESC'].split(';')

    l_new = []
    for i in l:
        if i == '':
            i = 'None'
        l_new.append(i)

    return l_new


def match_genomic(x):
    """
    convert HGVS genomic coordinates into BED format
    """
    pattern = re.compile('([1-9XY]+):g.([0-9]+)[AGTC>deldup]+')
    pattern_2 = re.compile('([1-9XY]+):g.([0-9]+)_([0-9]+)[delinsdupATGC]*')
    search = pattern.match(x['HGVS_G'])
    search_2 = pattern_2.match(x['HGVS_G'])

    if search:
        chr = search.group(1)
        start = str( int(search.group(2)) - 1 )
        end = search.group(2)
    elif search_2:
        chr = search_2.group(1)
        start = str( int(search_2.group(2)) - 1 )
        end = search_2.group(3)
    else:
        chr = '?'
        start = '?'
        end = '?'
    return chr, start, end


if __name__ == '__main__':
    args = get_args()

    # load config
    with open(args.config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # prefilter cosmic data
    cosmic_df = prefilter_cosmic_data(args.cosmic_file, config)

    # generate CSV output per referral type
    filter_by_referral(cosmic_df, config)
