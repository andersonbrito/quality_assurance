# -*- coding: utf-8 -*-
#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2021-01-24
# Last update: 2021-01-24

from Bio import SeqIO
import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate final QA matrix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--matrix", required=True, help="Latest version of QA matrix")
    parser.add_argument("--genomes", required=True, help="Complete genomes")
    parser.add_argument("--outliers", required=True, help="List of molecular clock outliers")
    parser.add_argument("--index", required=True, type=str,  help="Sample ID column name")
    parser.add_argument("--output1", required=True, help="Final version of QA matrix")
    parser.add_argument("--output2", required=True, help="Quality genomes")
    args = parser.parse_args()

    matrix = args.matrix
    genomes = args.genomes
    outliers = args.outliers
    index = args.index
    output1 = args.output1
    output2 = args.output2

    # matrix = path + 'qa/qa_matrix2.tsv'
    # outliers = path + 'root2tip/outliers_list.txt'
    # genomes = path + 'sequences/filtered_seqs.fasta'
    # index = 'SID'
    # output1 = path + 'assured_data/qa_matrix_final.tsv'
    # output2 = path + 'assured_data/genomes.tsv'

    # dict_rename = {}
    # for line in open(rename_list).readlines():
    #     sid, strain = line.strip().split('\t')
    #     dict_rename[strain] = sid


    # create a dict of genomes
    print('\nProcessing sequence file...\n')
    sequences = {}
    for fasta in SeqIO.parse(open(genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        sequences[id] = str(seq)


    # qa matrix
    df1 = pd.read_csv(matrix, encoding='utf-8', sep='\t', dtype='str')
    df1.fillna('', inplace=True)
    qa_entries = df1[index].to_list()

    # molecular clock outliers metadata
    values = [x.strip() for x in open(outliers).readlines()]
    data = {index: values}
    df2 = pd.DataFrame.from_dict(data)
    # df2 = df2.replace({index: dict_rename})
    df2['clock_status'] = ''
    df2 = df2[df2[index].isin(qa_entries)]

    # redefine indexes
    df1 = df1.set_index(index)
    df2 = df2.set_index(index)

    def failed(id):
        if id in df2.index.to_list(): # sequences that failed at root-to-tip analysis
            status = 'FAIL'
        # elif id not in dict_rename.values():
        # elif id not in sequences.keys():
        elif df1.loc[id, 'metadata_status'] == 'FAIL': # sequences that failed at metadata inspection
            status = 'FAIL'
        elif df1.loc[id, 'seq_coverage_status'] == 'FAIL': # sequences that failed at metadata inspection
            status = 'FAIL'
        else:
            status = 'PASS'
        return status

    # merge frames
    df3 = pd.concat([df1, df2[~df2.index.isin(df1.index)]]) # add new rows, if applicable
    df3.update(df2, overwrite=False)
    df3.reset_index(inplace=True)
    df3['clock_status'] = df3[index].map(failed)

    # export QA matrix
    df3.to_csv(output1, sep='\t', index=False)

    df3 = df3.set_index(index)
    with open(output2, 'w') as outfile:
        for id in sequences.keys():
            if df3.loc[id, 'clock_status'] == 'PASS': # filter out quality sequences
                entry = ">" + id + "\n" + sequences[id].upper() + "\n"
                outfile.write(entry)
