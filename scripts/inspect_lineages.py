# -*- coding: utf-8 -*-

import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate updated QA matrix, with mutation QA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--lineages", required=True, help="Pangolin lineages output file")
    parser.add_argument("--metadata", required=True, help="Latest metadata file")
    parser.add_argument("--matrix", required=True, help="Latest version of QA matrix")
    parser.add_argument("--index", required=True, type=str, help="Sample ID column name")
    parser.add_argument("--output1", required=True, help="Updated QA matrix")
    parser.add_argument("--output2", required=True, help="Final metadata file�")
    args = parser.parse_args()
    
    lineages = args.lineages
    metadata = args.metadata
    matrix = args.matrix
    index = args.index
    output1 = args.output1
    output2 = args.output2

#     path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_impacc/nextstrain/run9_20210621_pangoJSON/'
#     metadata = path + 'output_files/metadata/metadata1.tsv'
#     lineages = path + 'output_files/metadata/lineage_report.csv'
#     matrix = path + 'output_files/qa/qa_matrix3.tsv'
#     index = 'sample_id'
#     output1 = path + 'output_files/qa/qa_matrix4.tsv'
#     output2 = path + 'output_files/assured_data/metadata.tsv'


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
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df


    # load QA matrix
    df = load_table(matrix)
    df.fillna('', inplace=True)
    qa_entries = df[index].to_list()
    df = df.set_index(index)


    # load lineages file
    df2 = load_table(lineages)
    df2 = df2[df2['taxon'].isin(qa_entries)]

    # search lineage assignments
    dic_lineages = {}
    for idx, row in df2.iterrows():
        seqid = df2.loc[idx, 'taxon']
        lineage = df2.loc[idx, 'lineage']
        if seqid in qa_entries:
            dic_lineages[seqid] = lineage

    # store mutations in QA matrix
    df.reset_index(inplace=True)
    df['lineage'] = df[index].map(dic_lineages)
    df.fillna('', inplace=True)

    # assign PASS FAIL to mutations
    def failed(id):
        if id in dic_lineages:
            status = 'PASS'
        else:
            if df.loc[id, 'seq_coverage_status'] == 'FAIL':  # sequences that failed at metadata inspection
                status = 'NA'
            else:
                status = 'FAIL'
        return status

    df = df.set_index(index)
    df['lineage_status'] = df.index.map(failed)
    df.reset_index(inplace=True)

    # load metadata file
    df3 = load_table(metadata)
    df3.fillna('', inplace=True)
    df3['lineage'] = df3[index].map(dic_lineages)

    # export QA matrix
    df.to_csv(output1, sep='\t', index=False)

    # export metadata file
    df3.to_csv(output2, sep='\t', index=False)

    print('\nDone! TSV files successfully exported.\n')
