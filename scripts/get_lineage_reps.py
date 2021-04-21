#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2021-01-05
# Last update: 2021-01-21

import pandas as pd
import os
import argparse

pd.set_option('display.max_rows', 800)
pd.set_option('display.max_columns', 50)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Download list of SARS-CoV-2 lineage representatives",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--download", required=True, nargs=1, type=str, default='yes', choices=['yes', 'no'],
                        help="Download pangolin lineage file?")
    parser.add_argument("--keep-only", required=False, nargs='+', type=str,  help="List of country to filter in")
    parser.add_argument("--howmany", required=True, type=int,  help="How many representatives should be sampled per lineage?")
    parser.add_argument("--output", required=True, help="Output file")
    args = parser.parse_args()


    url = 'https://raw.githubusercontent.com/cov-lineages/pangoLEARN/master/pangoLEARN/data/lineages.metadata.csv'
    download = args.download[0]
    location_col = 'country'
    howmany = args.howmany
    keep_only = []
    if args.keep_only != None:
        keep_only = [col.strip() for col in args.keep_only[0].split()]
    output = args.output


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/nextstrain/runX_20210202_b117/quality_assurance/input_files/'
    # os.chdir(path)
    # download = 'yes'
    # location_col = 'country'
    # howmany = 4
    # # keep_only = None
    # keep_only = []
    # output = 'representatives_pangolin.txt'

    path = os.getcwd()
    if download == 'yes':
        if 'pangolin_lineages.csv' not in os.listdir(path):
            os.system('wget %s' % url)
            os.system('mv %s pangolin_lineages.csv' % url.split('/')[-1])
        df = pd.read_csv('pangolin_lineages.csv', encoding='utf-8', sep=',', dtype='str')
    else:
        df = pd.read_csv(url, encoding='utf-8', sep=',', dtype='str')

    if len(keep_only) > 0:
        df = df[df[location_col].isin(keep_only)]

    def fix_country_name(sequence_name):
        parts = sequence_name.split('/')
        country = parts[0]
        country = country.replace('_', '')
        new_name = country + '/' + '/'.join(parts[1:])
        return new_name

    df['sequence_name'] = df['sequence_name'].apply(lambda x: fix_country_name(x))
    lineages = df.groupby('lineage')
    df2 = pd.DataFrame(columns=['strain'])
    for name, dfL in lineages:
        available_samples = dfL['sequence_name'].count()  # genomes in bin

        if available_samples >= howmany:
            random_subset = dfL.sample(n=howmany)
            selected = random_subset['sequence_name'].to_list()
        else:
            random_subset = dfL.sample(n=1)
            selected = random_subset['sequence_name'].to_list()
        for entry in selected:
            dict_row = {}
            dict_row['strain'] = entry
            df2 = df2.append(dict_row, ignore_index=True)

    df2['strain'].sample(n=500, random_state=1).to_csv(output, sep='\t', index=False, header=False)
    print('\nList of genomes representative of pangolin lineages successfully generated\n')
