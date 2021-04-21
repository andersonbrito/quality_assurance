#!/usr/bin/python
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append newly sequenced genomes to current genome dataset, and export metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file with newly sequenced genomes")
    parser.add_argument("--index", required=True, type=str, help="Column with unique identifiers")
    parser.add_argument("--refgenome-size", required=True, type=int,  help="Reference genome size")
    parser.add_argument("--max-missing", required=True, type=int,  help="Maximum percentage of Ns or gaps (int = 1-100)")
    parser.add_argument("--output", required=True, help="Quality assurance TSV file")
    parser.add_argument("--output2", required=True, help="Renaming list")
    args = parser.parse_args()

    genomes = args.genomes
    index = args.index
    max_gaps = args.max_missing
    genome_size = args.refgenome_size
    output = args.output
    output2 = args.output2

    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_impacc/nextstrain/run8_20210402_impacc/"
    # genomes = path + "input_files/mock_sequences.fasta"
    # index = 'sample_id'
    # genome_size = 29903
    # max_gaps = 30
    # output = path + "output_files/qa_matrix1.tsv"
    # output2 = path + "output_files/rename.tsv"

    min_size = genome_size - int(genome_size * (max_gaps / 100))

    # store only new sequences in a dictionary, ignoring existing ones
    list_id = []
    list_visitnum = []
    list_txtpid = []
    list_seqid = []
    list_coverage = []
    list_status = []
    old_newname = {}
    for fasta in SeqIO.parse(open(genomes),'fasta'):
        id, seq = fasta.description, fasta.seq

        fullid = id
        visitnum  = id.split('|')[1].split(':')[1].split('_')[2]
        txtpid = id.split('|')[1].split(':')[1].split('_')[1]
        seqid = id.split('|')[0].split(':')[1].split('.')[0]
        id = id.split('|')[1].split(':')[1].split('_')[0]
        if fullid not in old_newname:
            old_newname[fullid] = id

        list_txtpid.append(txtpid)
        list_visitnum.append(visitnum)
        list_seqid.append(seqid)
        list_id.append(id)

        size = len(str(seq).replace('N', '').replace('-', ''))
        coverage = ''
        status = ''
        if size > min_size:
            coverage = str(round(size / genome_size, 3))
            status = 'PASS'
        else:
            coverage = str(round(size / genome_size, 3))
            status = 'FAIL'
        list_coverage.append(coverage)
        list_status.append(status)

    # create and export dataframe
    df = pd.DataFrame()
    df['sample_id'] = list_id
    df['participant_id'] = list_txtpid
    df['visit_num'] = list_visitnum
    df['core_specific_ID'] = list_seqid
    df['seq_coverage'] = list_coverage
    df['seq_coverage_status'] = list_status
    df.to_csv(output, sep='\t', index=False)

    # export renaming file
    with open(output2, 'w') as outfile:
        for old, new in old_newname.items():
            outfile.write(old + '\t' + new + '\n')

    # sequences filtered out due to low coverage
    print('\n### Filtering low coverage sequences ###\n')
    l = 1
    for idx, row in df.iterrows():
        id = df.loc[idx, index]
        status = df.loc[idx, 'seq_coverage_status']

        if status == 'FAIL':
            print('\t' + str(l) + '. ' + id)
            l += 1

    print('\n### Final result\n')
    print(str(len(list_id)) + ' sequences found in the original sequence file')
    print(str(l-1) + ' genomes were FILTERED OUT due low coverage (>' + str(max_gaps) + '% missing data)\n')
    print(str(len(list_id)-l+1) + ' genomes included in FINAL dataset\n')
