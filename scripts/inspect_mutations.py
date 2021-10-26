# -*- coding: utf-8 -*-

import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate updated QA matrix, with mutation QA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--matrix", required=True, help="Latest version of QA matrix")
    parser.add_argument("--index", required=True, type=str,  help="Sample ID column name")
    parser.add_argument("--mutations", required=True, help="Mutation TSV file")
    parser.add_argument("--insertions", required=True, help="Insertions CSV file from nextstrain")
    parser.add_argument("--output", required=True, help="Updated QA matrix")
    args = parser.parse_args()

    matrix = args.matrix
    index = args.index
    mutations = args.mutations
    insertions = args.insertions
    output = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_impacc/nextstrain/2021-05-17_update_pipeline/'
    # matrix = path + 'output_files/qa/qa_matrix2.tsv'
    # index = 'sample_id'
    # mutations = path + 'output_files/sequences/mutations.tsv'
    # insertions = path + 'output_files/sequences/nextalign.insertions.csv'
    # output = path + 'output_files/genomes.tsv'


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

    # load mutation file
    df2 = load_table(mutations)
    df2 = df2[df2['sequence'].isin(qa_entries)]
    seq_errors = df2[df2['mut_type'].isin(['nonstop', 'nonsense', 'frameshift'])]

    # redefine indexes
    df = df.set_index(index)
    # seq_errors = seq_errors.set_index('sequence')
    dic_errors = {}
    for idx, row in seq_errors.iterrows():
        seqid = seq_errors.loc[idx, 'sequence']
        if seqid in qa_entries:
            comment = seq_errors.loc[idx, 'mut_type'] + '(' + seq_errors.loc[idx, 'position'].split('..')[0] + ':' + seq_errors.loc[idx, 'var_codon'] + ')'
            if seqid not in dic_errors:
                dic_errors[seqid] = comment
            else:
                dic_errors[seqid] += ";"+comment


    # load insertion file
    df3 = load_table(insertions)
    df3 = df3[df3['seqName'].isin(qa_entries)]
    df3.fillna('', inplace=True)
    insert_errors = df3[~df3['insertions'].isin([''])]

    for idx, row in insert_errors.iterrows():
        seqid = insert_errors.loc[idx, 'seqName']
        if seqid in qa_entries:
            comment = 'insertion(' + insert_errors.loc[idx, 'insertions'] + ')'
            if seqid not in dic_errors:
                dic_errors[seqid] = comment
            else:
                dic_errors[seqid] += ";"+comment


    # store mutations in QA matrix
    df.reset_index(inplace=True)
    df['mutations'] = df[index].map(dic_errors)
    df.fillna('', inplace=True)

    # assign PASS FAIL to mutations
    def failed(id):
        if id in dic_errors:
            status = 'FAIL'
        else:
            if df.loc[id, 'seq_coverage_status'] == 'FAIL':  # sequences that failed at metadata inspection
                status = 'NA'
            else:
                status = 'PASS'
        return status
    
    dic_errors = {}
    df = df.set_index(index)
    df['mutation_status'] = df.index.map(failed)
    df.reset_index(inplace=True)

    # export QA matrix
    df.to_csv(output, sep='\t', index=False)
