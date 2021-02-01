# -*- coding: utf-8 -*-
#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2021-01-25
# Last update: 2021-01-25

from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file with new genomes")
    parser.add_argument("--metadata", required=True, help="Sample metadata file")
    parser.add_argument("--index", required=True, type=str, help="Column with unique identifiers")
    parser.add_argument("--matrix", required=False, help="Quality assurance matrix file")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="Expanded quality assurance matrix file")
    args = parser.parse_args()

    genomes = args.genomes
    metadata = args.metadata
    index = args.index
    matrix = args.matrix
    output1 = args.output1
    output2 = args.output2


    # genomes = path + 'input_files/new_genomes.fasta'
    # metadata = path + 'input_files/querydoc_ISMMS.csv'
    # index = 'SID'
    # matrix = path + 'output_files/qa_matrix1.tsv'
    # output1 = path + 'output_files/metadata_samples.tsv'
    # output2 = path + 'output_files/qa_matrix2.tsv'


    pd.set_option('max_columns', 100)

    # create a dict of existing sequences
    print('\nProcessing sequence file...\n')
    missing = {}
    for fasta in SeqIO.parse(open(genomes), 'fasta'):
        id = fasta.description
        if id not in missing.keys():
            missing[id] = []


    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
            df = df.rename(columns={'Sample-ID': 'id'})
            df = df[~df[index].isin([''])]  # drop row with empty index
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df


    # Load QA matrix
    dfQ = load_table(matrix)

    # Load sample metadata
    dfS = load_table(metadata)
    dfS['region'] = ''
    # dfS = dfS.rename(columns={'Date': 'date', 'Country': 'country', 'State': 'division', 'Site': 'originating_lab'})
    dfS = dfS.rename(columns={'Sample-ID': 'id', 'Collection-date': 'date', 'Country': 'country', 'Division': 'division',
                          'State': 'code', 'Location': 'location', 'Lineage': 'pangolin_lineage', 'Source': 'originating_lab',
                          'Update': 'update'})
    if 'update' in dfS.columns.to_list():
        dfS = dfS[~dfS['update'].isin([''])]  # drop row with empty update information

    dfS['date'] = dfS['date'].replace({'?': ''})
    dfS.fillna('', inplace=True)

    print('Inspecting metadata columns...\n')
    # find sequences with no metadata
    no_metadata = [x for x in missing.keys() if x not in dfS[index].to_list()]
    for id in no_metadata:
        missing[id] = ['all fields']

    # fix place of origin when disagreements between 'place' and 'place_exposure' exists, and flag missing geodata
    def fix_exposure(df):
        geo_columns = ['country', 'division']
        for level in geo_columns:
            exposure_column = level + '_exposure'
            if exposure_column not in df.columns.to_list():
                df[exposure_column] = ''

            for idx, row in df.iterrows():
                # print(df.loc[idx, index])
                if df.loc[idx, level] in ['', 'unknown', np.nan, None]:
                    if df.loc[idx, index] not in missing:
                        missing[df.loc[idx, index]] = [level]
                    else:
                        missing[df.loc[idx, index]] += [level]
                else:
                    if df.loc[idx, index] not in missing:
                        missing[df.loc[idx, index]] = []

                if df.loc[idx, exposure_column].lower() in ['', 'unknown', np.nan, None]:
                    df.loc[idx, exposure_column] = df.loc[idx, level]
        return df
    dfS = fix_exposure(dfS)


    # add state code

    # us_state_abbrev = {
    #     'Alabama': 'AL',
    #     'Alaska': 'AK',
    #     'American Samoa': 'AS',
    #     'Arizona': 'AZ',
    #     'Arkansas': 'AR',
    #     'California': 'CA',
    #     'Colorado': 'CO',
    #     'Connecticut': 'CT',
    #     'Delaware': 'DE',
    #     'District of Columbia': 'DC',
    #     'Washington DC': 'DC',
    #     'Florida': 'FL',
    #     'Georgia': 'GA',
    #     'Guam': 'GU',
    #     'Hawaii': 'HI',
    #     'Idaho': 'ID',
    #     'Illinois': 'IL',
    #     'Indiana': 'IN',
    #     'Iowa': 'IA',
    #     'Kansas': 'KS',
    #     'Kentucky': 'KY',
    #     'Louisiana': 'LA',
    #     'Maine': 'ME',
    #     'Maryland': 'MD',
    #     'Massachusetts': 'MA',
    #     'Michigan': 'MI',
    #     'Minnesota': 'MN',
    #     'Mississippi': 'MS',
    #     'Missouri': 'MO',
    #     'Montana': 'MT',
    #     'Nebraska': 'NE',
    #     'Nevada': 'NV',
    #     'New Hampshire': 'NH',
    #     'New Jersey': 'NJ',
    #     'New Mexico': 'NM',
    #     'New York': 'NY',
    #     'North Carolina': 'NC',
    #     'North Dakota': 'ND',
    #     'Northern Mariana Islands': 'MP',
    #     'Ohio': 'OH',
    #     'Oklahoma': 'OK',
    #     'Oregon': 'OR',
    #     'Pennsylvania': 'PA',
    #     'Puerto Rico': 'PR',
    #     'Rhode Island': 'RI',
    #     'South Carolina': 'SC',
    #     'South Dakota': 'SD',
    #     'Tennessee': 'TN',
    #     'Texas': 'TX',
    #     'Utah': 'UT',
    #     'Vermont': 'VT',
    #     'Virgin Islands': 'VI',
    #     'Virginia': 'VA',
    #     'Washington': 'WA',
    #     'West Virginia': 'WV',
    #     'Wisconsin': 'WI',
    #     'Wyoming': 'WY'
    # }
    # if 'state_code' not in dfS.columns.to_list():
    #     dfS.insert(1, 'state_code', '')
    #     dfS['state_code'] = dfS['division_exposure'].apply(lambda x: us_state_abbrev[x] if x in us_state_abbrev else 'un')


    # find incomplete dates
    for idx, row in dfS.iterrows():
        date = dfS.loc[idx, 'date']
        if date in ['', 'unknown', np.nan, None] or len(date) < 10:
            if dfS.loc[idx, index] not in missing:
                missing[dfS.loc[idx, index]] = ['date']
            else:
                missing[dfS.loc[idx, index]] += ['date']
        else:
            if dfS.loc[idx, index] not in missing:
                missing[dfS.loc[idx, index]] = []

    # fix date format
    dfS['date'] = pd.to_datetime(dfS['date'])
    dfS['date'] = dfS['date'].dt.strftime('%Y-%m-%d')

    if 'strain' not in dfS.columns.to_list():
        dfS.insert(1, 'strain', '')
        dfS['strain'] = dfS[index]
        # dfS['new_strain_name'] = dfS['country_exposure'].astype(str) + '/' + dfS['state_code'].astype(str) + '-IMPACC-' + \
        #                 dfS[index].str.split('-').str[1] + '/' + dfS['date'].str.split('-').str[0]

    # empty matrix dataframe
    dfM = pd.DataFrame()
    missing_data = {x: ', '.join(missing[x]) for x in missing}
    dfM[index] = missing_data.keys()
    dfM['missing_metadata'] = missing_data.values()

    no_genome = [x for x in dfS[index] if x not in missing.keys()]
    dfM['metadata_status'] = dfM['missing_metadata'].apply(lambda x: 'FAIL' if len(x) > 0 else 'PASS')


    # redefine indexes
    dfQ = dfQ.set_index(index)
    dfM = dfM.set_index(index)

    result = pd.concat([dfQ, dfM], axis=1, sort=False)

    result.update(dfM, overwrite=False)
    result.reset_index(inplace=True)
    result = result.rename(columns={"index": index})

    result.fillna('', inplace=True)
    result.loc[result['seq_coverage'] == '', 'seq_coverage_status'] = 'FAIL' # if value in column X is '', column Y = y

    # filter sequences with full metadata
    filter1 = result['seq_coverage_status'] == 'PASS'
    filter2 = result['metadata_status'] == 'PASS'
    # result.where(filter1 & filter2, inplace=False)
    full_metadata = result[(filter1) & (filter2)][index].to_list()

    # export only complete metadata rows
    dfS = dfS[dfS[index].isin(full_metadata)]
    dfS.to_csv(output1, sep='\t', index=False)

    # export QA matrix
    result.to_csv(output2, sep='\t', index=False)

    print('\nDone! TSV files successfully exported.\n')

