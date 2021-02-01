#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2021-01-06
# Last update: 2021-01-06


import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append newly sequenced genomes to current genome dataset, and export metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file with latest genomes from GISAID")
    parser.add_argument("--new-genomes", required=True, help="FASTA file with newly sequenced genomes")
    parser.add_argument("--keep", required=True, help="TXT file with accession number of genomes to be included")
    parser.add_argument("--remove", required=True, help="TXT file with accession number of genomes to be removed")
    parser.add_argument("--refgenome-size", required=True, type=int,  help="Reference genome size")
    parser.add_argument("--max-missing", required=True, type=int,  help="Maximum percentage of Ns or gaps (int = 1-100)")
    parser.add_argument("--output1", required=True, help="FASTA file containing lineage representative and new genomes")
    parser.add_argument("--output2", required=True, help="FASTA file containing low coverage genomes")
    args = parser.parse_args()

    genomes = args.genomes
    new_genomes = args.new_genomes
    keep = args.keep
    remove = args.remove
    max_gaps = args.max_missing
    genome_size = args.refgenome_size
    outfile1 = args.output1
    outfile2 = args.output2


    # genomes = path + "gisaid_hcov-19.fasta"
    # new_genomes = path + "new_genomes.fasta"
    # keep = path + 'keep.txt'
    # remove = path + "remove.txt"
    # genome_size = 29420
    # max_gaps = 5
    # outfile1 = path + "sequences_temp.fasta"
    # outfile2 = path + 'low_coverage.fasta'

    min_size = genome_size - int(genome_size * max_gaps / 100)

    # store only new sequences in a dictionary, ignoring existing ones
    newly_sequenced = {}
    short_sequences = {}
    for fasta in SeqIO.parse(open(new_genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        if '/' in id:
            id = id.split('/')[0]
        size = len(str(seq).replace('N', '').replace('-', ''))
        if size > min_size:
            if id not in newly_sequenced.keys(): # avoid potential duplicates
                newly_sequenced[id] = str(seq)
        else:
            if id not in short_sequences.keys(): # avoid potential duplicates
                short_sequences[id] = str(seq)

    # create a list of the existing sequences
    all_sequences = {}
    for fasta in SeqIO.parse(open(genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        # print(id)
        id = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        seq = str(seq).replace('N', '').replace('-', '')
        size = len(seq)
        if size > min_size:
            all_sequences[id] = str(seq)

    # create a list of sequences to be added in all instances
    keep_sequences = {}
    mismatch = []
    for id in sorted(open(keep, "r").readlines()):
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            if id not in keep_sequences.keys():
                try:
                    keep_sequences[id] = all_sequences[id]
                except:
                    mismatch.append(id)

    # create a list of sequences to be ignored in all instances
    remove_sequences = []
    for id in open(remove, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            remove_sequences.append(id)


    # export only sequences to be used in the nextstrain build
    c = 1
    removed = 0
    sequences = {**keep_sequences, **newly_sequenced}
    print('\n### Exporting sequences\n')
    exported = []
    with open(outfile1, 'w') as output1:
        for id in sequences.keys():
            if id not in remove_sequences: # filter out unwanted sequences
                entry = ">" + id + "\n" + sequences[id].upper() + "\n"
                exported.append(id)
                output1.write(entry)
                if '/' not in id: # search for newly sequenced genomes, named "F999"
                    print('* ' + str(c) + '. ' + id)
                else:
                    print(str(c) + '. ' + id)
            else:
                removed += 1
                c -= 1
            c += 1

    # export low coverage sequences
    filtered_out = []
    with open(outfile2, 'w') as output2:
        for id in short_sequences.keys():
            entry = ">" + id + "\n" + short_sequences[id].upper() + "\n"
            filtered_out.append(id)
            output2.write(entry)

    # mismatched sequence headers
    m = 1
    if len(mismatch) > 0:
        print('\n### List of genomes not found in input sequence file\n')
        for id in mismatch:
            print(str(m) + '. ' + id)
            m += 1
    else:
        print('\tAll requested genomes were successfully incorporated in the analysis...')


    # excluding sequences
    print('\n### Excluding sequences ###\n')
    e = 1
    for id in remove_sequences:
        if id in sequences:
            print(str(e) + '. ' + id)
            e += 1

    # sequences filtered out due to low coverage
    print('\n### Filtering low coverage sequences out ###\n')
    l = 1
    for id in filtered_out:
        if id in short_sequences:
            print(str(l) + '. ' + id)
            l += 1


    print('\n### Final result\n')

    print('Lab file contains ' + str(len(newly_sequenced)) + ' sequences')
    print('GISAID file contains ' + str(len(all_sequences)) + ' sequences\n')

    print(str(len(mismatch)) + ' genomes in keep.txt were NOT FOUND on GISAID database')
    print(str(len(keep_sequences)) + ' genomes ADDED from GISAID dataset')
    print(str(len(newly_sequenced)) + ' newly sequenced genomes were added')
    print(str(removed) + ' genomes were REMOVED according to remove.txt')
    print(str(len(filtered_out)) + ' genomes were FILTERED OUT due low coverage (>' + str(max_gaps) + '% missing data)\n')
    print(str(len(exported)) + ' genomes included in FINAL dataset\n')

