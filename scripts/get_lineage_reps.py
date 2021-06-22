#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2021-01-05
# Last update: 2021-06-22

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
    parser.add_argument("--output", required=True, help="Output file")
    args = parser.parse_args()


    url = 'https://raw.githubusercontent.com/cov-lineages/pangoLEARN/master/pangoLEARN/data/lineages.downsample.csv'
    download = args.download[0]
    output = args.output


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_impacc/nextstrain/run9_20210621_pangoJSON/input_files/'
    # os.chdir(path)
    # download = 'yes'
    # output = 'representatives_pangolin.txt'

    path = os.getcwd()
    if download == 'yes':
        if 'pangolin_lineages.csv' not in os.listdir(path):
            os.system('wget %s' % url)
            os.system('mv %s pangolin_lineages.csv' % url.split('/')[-1])
        df = pd.read_csv('pangolin_lineages.csv', encoding='utf-8', sep=',', dtype='str')
    else:
        df = pd.read_csv(url, encoding='utf-8', sep=',', dtype='str')

    def fix_country_name(sequence_name):
        parts = sequence_name.split('/')
        country = parts[0]
        country = country.replace('_', '').replace('\'', '-')
        new_name = country + '/' + '/'.join(parts[1:])
        return new_name

    df['sequence_name'] = df['sequence_name'].apply(lambda x: fix_country_name(x))
    lineages = df.groupby('lineage')
    df2 = pd.DataFrame(columns=['strain'])
    for name, dfL in lineages:
        available_samples = dfL['sequence_name'].count()  # genomes in bin

        random_subset = dfL.sample(n=1)
        selected = random_subset['sequence_name'].to_list()
        for entry in selected:
            dict_row = {}
            dict_row['strain'] = entry
            df2 = df2.append(dict_row, ignore_index=True)

    df2['strain'].sample(n=500, random_state=1).to_csv(output, sep='\t', index=False, header=False)
    print('\nList of genomes representative of pangolin lineages successfully generated\n')
