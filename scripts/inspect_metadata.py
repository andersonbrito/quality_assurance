# -*- coding: utf-8 -*-
#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2021-01-25
# Last update: 2021-04-02

import pandas as pd
import numpy as np
import argparse
from datetime import datetime


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata1", required=True, help="Core lab metadata file")
    parser.add_argument("--metadata2", required=True, help="Sample metadata file")
    parser.add_argument("--metadata3", required=True, help="Patient metadata file")
    parser.add_argument("--batch", required=True, help="Batch layout file")
    parser.add_argument("--index", required=True, type=str, help="Column with unique identifiers")
    parser.add_argument("--matrix", required=False, help="Quality assurance matrix file")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="Expanded quality assurance matrix file")
    parser.add_argument("--output3", required=True, help="Renaming file 2")
    args = parser.parse_args()

    cl_metadata = args.metadata1
    cs_metadata = args.metadata2
    ci_metadata = args.metadata3
    bl_metadata = args.batch
    index = args.index
    matrix = args.matrix
    output1 = args.output1
    output2 = args.output2
    output3 = args.output3


    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_impacc/nextstrain/run8_20210402_impacc/"
    # cl_metadata = path + 'input_files/metadata_impacc.csv'
    # cs_metadata = path + 'input_files/impacc-virology-clin-sample.csv'
    # ci_metadata = path + 'input_files/impacc-virology-clin-individ.csv'
    # bl_metadata = path + 'input_files/batch_layout.csv'
    # index = 'sample_id'
    # matrix = path + 'output_files/qa/qa_matrix1.tsv'
    # output1 = path + 'output_files/metadata_samples.tsv'
    # output2 = path + 'output_files/qa_matrix2.tsv'
    # output3 = path + 'output_files/rename2.tsv'

    pd.set_option('max_columns', 100)


    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
            df = df.rename(columns={'guspec': index, 'txtpid': 'participant_id'})
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
            df = df.rename(columns={'guspec': index, 'txtpid': 'participant_id'})
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
            df = df[~df[index].isin([''])]  # drop row with empty index
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df


    # Load sample metadata
    dfL = load_table(cl_metadata)
    dfL.fillna('', inplace=True)
    dfL['visit_num'] = dfL['visit_num'].apply(lambda x: x.split(' ')[1] if ' ' in x else x)

    print('\n### Inspecting metadata columns...\n')

    # load clinical sample metadata
    def year_week(y, w):
        return datetime.strptime(f'{y} {w} 3', '%G %V %u')

    dfCS = load_table(cs_metadata)
    dfCS['date'] = dfCS.apply(lambda row: year_week(row['calendar_year'], row['calendar_week']), axis=1)

    # fix date format
    dfCS['date'] = pd.to_datetime(dfCS['date'])
    dfCS['date'] = dfCS['date'].dt.strftime('%Y-%m-%d')


    # load clinical individual metadata
    dfCI = load_table(ci_metadata)
    enrollment_state = {'Arizona': 'Arizona', 'Baylor': 'Texas', 'Boston/BWH': 'Massachusetts', 'Case Western': 'Ohio',
                'Drexel/Tower Health': 'Pennsylvania', 'Emory': 'Georgia', 'Florida': 'Florida', 'Yale': 'Connecticut',
                'ISMMS (Mt Sinai)': 'New York', 'OHSU (Oregon)': 'Oregon', 'OUHSC (Oklahoma)': 'Oklahoma',
                'Stanford': 'California', 'UCLA': 'California', 'UCSF': 'California', 'UT Austin': 'Texas'}
    dfCI['region'] = 'North America'
    dfCI['country'] = 'USA'
    dfCI['division'] = dfCI['enrollment_site'].apply(lambda x: enrollment_state[x] if x in enrollment_state else '')


    # add location
    def variant_category(site):
        state = ''
        if site in enrollment_state:
            state = enrollment_state[site]
        return state

    # melt metadata files
    def melt_metadata(dataframe1, dataframe2, extra_columns):
        if len(extra_columns) > 0:
            print('\t- Adding extra columns to dataframe: ' + ', '.join(extra_columns))

        dataframe2 = dataframe2[extra_columns]
        for column in extra_columns:
            dataframe1[column] = ''
            for idx, row in dataframe2.iterrows():
                if idx in dataframe1.index:
                    dataframe1.loc[idx, column] = dataframe2.loc[idx, column]
        dataframe1.reset_index(inplace=True)
        return dataframe1


    # redefine indexes, round 1
    dfL = dfL.set_index(index)
    dfCS = dfCS.set_index(index)
    dfL = melt_metadata(dfL, dfCS, ['date', 'trajgroup'])

    # redefine indexes, round 2
    dfL = dfL.set_index('participant_id')
    dfCI = dfCI.set_index('participant_id')
    dfL = melt_metadata(dfL, dfCI, ['region', 'country', 'division', 'enrollment_site', 'admit_age', 'sex', 'race', 'ethnicity', 'death'])


    # Load QA matrix
    dfQ = load_table(matrix)
    # create a dict of existing sequences
    print('\n\n### Processing QA matrix...\n')
    missing = {}
    list_id = dfQ[index].tolist()
    for id in list_id:
        missing[id] = []


    # add missing metadata into QA matrix
    for idx, row in dfL.iterrows():
        sample_id = dfL.loc[idx, index]
        dict_row = {}
        if sample_id not in dfQ[index].tolist():# and sample_id not in dict_row[index]:# and sample_id not in found:
            for column in dfQ.columns.tolist():
                if column in dfL.columns.tolist():
                    dict_row[column] = dfL.loc[dfL[index] == sample_id, column].iloc[0]
            dfQ = dfQ.append(dict_row, ignore_index=True)


    print('\n### Processing batch layout...\n')

    # load batch layout
    dfB = load_table(bl_metadata)
    dfB['visit_num'] = dfB['visit_num'].apply(lambda x: x.split(' ')[1] if ' ' in x else x)

    dfL = dfL.set_index(index)
    dfB = dfB.set_index(index)


    # detect batch layout inconsistencies
    batch_issues = {}
    for idx, row in dfL.iterrows():
        if idx not in dfB.index.tolist():
            if idx not in batch_issues:
                batch_issues[idx] = ['all_fields']
            else:
                batch_issues[idx] += ['all_fields']
        else:
            issues = []
            if dfL.loc[idx, 'participant_id'] != dfB.loc[idx, 'participant_id']:
                issues.append('participant_id')
            if dfL.loc[idx, 'visit_num'] != dfB.loc[idx, 'visit_num']:
                issues.append('visit_num')

            if idx not in batch_issues:
                batch_issues[idx] = issues
            else:
                batch_issues[idx] += issues

    dfL.reset_index(inplace=True)
    dfQ.reset_index(inplace=True)


    # find sequences with no metadata
    no_metadata = [x for x in missing.keys() if x not in dfL[index].tolist()]
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
    dfL = fix_exposure(dfL)


    # add state code
    us_state_abbrev = {
        'Alabama': 'AL',
        'Alaska': 'AK',
        'American Samoa': 'AS',
        'Arizona': 'AZ',
        'Arkansas': 'AR',
        'California': 'CA',
        'Colorado': 'CO',
        'Connecticut': 'CT',
        'Delaware': 'DE',
        'District of Columbia': 'DC',
        'Washington DC': 'DC',
        'Florida': 'FL',
        'Georgia': 'GA',
        'Guam': 'GU',
        'Hawaii': 'HI',
        'Idaho': 'ID',
        'Illinois': 'IL',
        'Indiana': 'IN',
        'Iowa': 'IA',
        'Kansas': 'KS',
        'Kentucky': 'KY',
        'Louisiana': 'LA',
        'Maine': 'ME',
        'Maryland': 'MD',
        'Massachusetts': 'MA',
        'Michigan': 'MI',
        'Minnesota': 'MN',
        'Mississippi': 'MS',
        'Missouri': 'MO',
        'Montana': 'MT',
        'Nebraska': 'NE',
        'Nevada': 'NV',
        'New Hampshire': 'NH',
        'New Jersey': 'NJ',
        'New Mexico': 'NM',
        'New York': 'NY',
        'North Carolina': 'NC',
        'North Dakota': 'ND',
        'Northern Mariana Islands': 'MP',
        'Ohio': 'OH',
        'Oklahoma': 'OK',
        'Oregon': 'OR',
        'Pennsylvania': 'PA',
        'Puerto Rico': 'PR',
        'Rhode Island': 'RI',
        'South Carolina': 'SC',
        'South Dakota': 'SD',
        'Tennessee': 'TN',
        'Texas': 'TX',
        'Utah': 'UT',
        'Vermont': 'VT',
        'Virgin Islands': 'VI',
        'Virginia': 'VA',
        'Washington': 'WA',
        'West Virginia': 'WV',
        'Wisconsin': 'WI',
        'Wyoming': 'WY'
    }
    if 'state_code' not in dfL.columns.to_list():
        dfL['state_code'] = dfL['division_exposure'].apply(lambda x: us_state_abbrev[x] if x in us_state_abbrev else 'un')

    # find incomplete dates
    for idx, row in dfL.iterrows():
        date = dfL.loc[idx, 'date']
        if date in ['', 'unknown', np.nan, None] or len(date) < 10:
            if dfL.loc[idx, index] not in missing:
                missing[dfL.loc[idx, index]] = ['date']
            else:
                missing[dfL.loc[idx, index]] += ['date']
        else:
            if dfL.loc[idx, index] not in missing:
                missing[dfL.loc[idx, index]] = []


    if 'strain' not in dfL.columns.to_list():
        dfL.insert(0, 'strain', '')
        # dfL['strain'] = dfL[index]
        dfL['strain'] = dfL['country_exposure'].astype(str) + '/' + dfL['state_code'].astype(str) + '-' + \
                        dfL['core_specific_ID'] + '/' + dfL['date'].str.split('-').str[0]

    dfR = dfL[['sample_id', 'strain']]
    dfR.to_csv(output3, sep='\t', index=False, header=False)


    # empty matrix dataframe
    def qa_results(dict_records, missing_col, status_col):
        df = pd.DataFrame()
        missing_data = {x: ', '.join(dict_records[x]) for x in dict_records}
        df[index] = missing_data.keys()
        df[missing_col] = missing_data.values()
        df[status_col] = df[missing_col].apply(lambda x: 'FAIL' if len(x) > 0 else 'PASS')
        return df


    # add metadata QA status
    dfMM = qa_results(missing, 'missing_metadata', 'metadata_status')
    dfQ = dfQ.set_index(index)
    dfMM = dfMM.set_index(index)
    dfQ = melt_metadata(dfQ, dfMM, ['missing_metadata', 'metadata_status'])


    # add batch layout QA status
    dfBL = qa_results(batch_issues, 'batch_issues', 'batch_status')
    dfQ = dfQ.set_index(index)
    dfBL = dfBL.set_index(index)
    dfQ = melt_metadata(dfQ, dfBL, ['batch_issues', 'batch_status'])


    result = dfQ
    result.fillna('', inplace=True)
    result.loc[result['seq_coverage'] == '', 'seq_coverage_status'] = 'FAIL' # if value in column X is '', column Y = y

    # filter sequences with full metadata
    filter1 = result['seq_coverage_status'] == 'PASS'
    filter2 = result['metadata_status'] == 'PASS'
    filter3 = result['batch_status'] == 'PASS'


    # result.where(filter1 & filter2, inplace=False)
    full_metadata = result[(filter1) & (filter2) & (filter3)][index].to_list()

    # export only complete metadata rows
    dfL = dfL[dfL[index].isin(full_metadata)]
    dfL.to_csv(output1, sep='\t', index=False)

    # export QA matrix
    result.to_csv(output2, sep='\t', index=False)



    print('\nDone! TSV files successfully exported.\n')
