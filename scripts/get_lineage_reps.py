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
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--download", required=True, nargs=1, type=str, default='yes', choices=['yes', 'no'],
                        help="Download pangolin lineage file?")
    parser.add_argument("--keep-only", required=False, nargs='+', type=str,  help="List of country to filter in")
    parser.add_argument("--output", required=True, help="Output file")
    args = parser.parse_args()


    url = 'https://raw.githubusercontent.com/cov-lineages/pangoLEARN/master/pangoLEARN/data/lineages.metadata.csv'
    download = args.download[0]
    location_col = 'country'
    keep_only = []
    if args.keep_only != None:
        keep_only = [col.strip() for col in args.keep_only[0].split()]
    output = args.output


    # os.chdir(path)
    # download_file = 'yes'
    # keep_only = None
    # keep_only = ['Brazil']
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
    df = df.sort_values(by='sample_date').drop_duplicates(subset=['lineage'])
    # print(df.sort_values(by='lineage'))

    def fix_country_name(sequence_name):
        parts = sequence_name.split('/')
        country = parts[0]
        country = country.replace('_', '')
        new_name = country + '/' + '/'.join(parts[1:])
        return new_name

    print('\nSaved list of genomes representative of pangolin lineages\n')
    df['sequence_name'] = df['sequence_name'].apply(lambda x: fix_country_name(x))
    df['sequence_name'].to_csv(output, sep='\t', index=False, header=False)
