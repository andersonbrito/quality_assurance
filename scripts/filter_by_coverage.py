#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2021-01-25
# Last update: 2021-01-25

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
    args = parser.parse_args()

    genomes = args.genomes
    index = args.index
    max_gaps = args.max_missing
    genome_size = args.refgenome_size
    output = args.output

    # genomes = path + "input_files/new_genomes.fasta"
    # genome_size = 29903
    # max_gaps = 10
    # output = path + "output_files/qa_matrix1.tsv"

    min_size = genome_size - int(genome_size * (max_gaps / 100))

    # store only new sequences in a dictionary, ignoring existing ones
    high_coverage = {}
    low_coverage = {}
    for fasta in SeqIO.parse(open(genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        size = len(str(seq).replace('N', '').replace('-', ''))
        if size > min_size:
            if id not in high_coverage: # avoid potential duplicates
                high_coverage[id] = str(round(size/genome_size, 3))
        else:
            if id not in low_coverage: # avoid potential duplicates
                low_coverage[id] = str(round(size/genome_size, 3))

    # create and export dataframe
    df = pd.DataFrame()
    df[index] = list(high_coverage.keys()) + list(low_coverage.keys())
    df['seq_coverage'] = list(high_coverage.values()) + list(low_coverage.values())
    df['seq_coverage_status'] = ['PASS' for x in high_coverage] + ['FAIL' for x in low_coverage]
    df.to_csv(output, sep='\t', index=False)


    # sequences filtered out due to low coverage
    print('\n### Filtering low coverage sequences out ###\n')
    l = 1
    for id in low_coverage:
        print('\t' + str(l) + '. ' + id)
        l += 1

    print('\n### Final result\n')
    print(str(len(high_coverage)+len(low_coverage)) + ' sequences found in the original sequence file')
    print(str(len(low_coverage)) + ' genomes were FILTERED OUT due low coverage (>' + str(max_gaps) + '% missing data)\n')
    print(str(len(high_coverage)) + ' genomes included in FINAL dataset\n')

