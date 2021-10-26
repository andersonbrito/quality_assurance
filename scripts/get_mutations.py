#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2020-12-16
# Last update: 2021-06-22

from Bio import SeqIO
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
from Bio.Alphabet import Gapped
from Bio.Alphabet import IUPAC
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import sys
import argparse

pd.set_option('display.max_columns', 500)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Search for genetic changes and potential sequencing errors.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="FASTA file with genome alignment")
    parser.add_argument("--annotation", required=True, help="Reference genome annotation file")
    parser.add_argument("--ref-name", required=True, help="Name of the reference sequence, as shown in the alignment")
    parser.add_argument("--export-only", required=False, nargs='+', type=str,
                        choices=['synonymous', 'nonsynonymous','insertion', 'deletion', 'nonstop', 'nonsense', 'frameshift', 'ambiguity'],
                        help="Categories to be exclusively exported")
    parser.add_argument("--ignore", required=False, nargs='+', type=str,
                        choices=['synonymous', 'nonsynonymous','insertion', 'deletion', 'nonstop', 'nonsense', 'frameshift', 'ambiguity'],
                        help="Categories to be ignored (it only works if '--export-only' is OFF)")
    parser.add_argument("--output", required=True, help="TSV file with variants")
    args = parser.parse_args()


    alignment = args.input
    reference_file = args.annotation
    ref_genome = args.ref_name
    export_list = args.export_only
    ignore_list = args.ignore
    output = args.output

    export_only = []
    ignore = []
    if export_list != None:
        export_only = [col.strip() for col in export_list[0].split()]
    if ignore_list != None:
        ignore = [col.strip() for col in ignore_list[0].split()]


    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_impacc/nextstrain/2021-05-17_update_pipeline/"
    # alignment = path + 'output_files/sequences/aligned.fasta'
    # reference_file = path + 'input_files/reference.gb'
    # annotation = SeqIO.read(reference_file, "genbank") # annotation file
    # ref_genome = 'Wuhan/Hu-1/2019'
    # export_only = ['nonstop', 'nonsense', 'frameshift']
    # ignore = []
    # output = path + 'variants.tsv'


    """ ANNOTATIONS """
    # pull annotation coordinates
    list_traits = []
    for seq_record in SeqIO.parse(reference_file, "gb"):
        for feature in seq_record.features:
            # print(feature.qualifiers)
            if feature.type == 'CDS':# or feature.type == 'probe':
                name = str(feature.qualifiers['gene'][0])
                positions = ''
                # print(str(feature.location))

                if 'join' in str(feature.location):
                    rgx = re.search(r'(\[.+\))', str(feature.location))
                    if rgx:
                        # print(rgx.group(0))
                        parts = rgx.group(0).split(",")
                        for part in parts:
                            start = int(part.strip().split(']')[0][1:].split(':')[0])
                            end = int(part.strip().split(']')[0][1:].split(':')[1])
                            sense = part.split('(')[1][:-1]
                            # print(name, sense, positions)
                            if sense == '-':
                                start = start + 1
                            positions = str(start + 1) + ':' + str(end)

                            id = name + '|' + sense + '|' + positions
                            if id not in list_traits:
                                list_traits.append(id)
                else:
                    sense = str(feature.location).split('(')[1][:-1]
                    start = int(str(feature.location).split(']')[0][1:].split(':')[0])
                    end = int(str(feature.location).split(']')[0][1:].split(':')[1])
                    # if sense == '-':
                    #     start = start + 1
                    positions = str(start + 1) + ':' + str(end)
                    id = name + '|' + sense + '|' + positions
                    if id not in list_traits:
                        list_traits.append(id)


    if len(list_traits) > 0:
        print('\nGenomic traits being scanned = ' + str(len(list_traits)) + ':\n')
        for t in list_traits:
            print('\t- ' + t)
    else:
        print('No genomic annotation was found. Stopping mutation detection...\n')
        sys.exit()



    """ CODONS """
    # Read original sequence file, and create list of codons
    dfCodons = pd.DataFrame(columns=['genomes'] + list_traits)

    # search for insertion in reference genome
    seq_index = SeqIO.index(alignment, 'fasta')
    ref_seq = seq_index[ref_genome].seq
    insert_pos = [pos for pos, nuc in enumerate(ref_seq, 1) if nuc == '-']

    insertions = {}
    if len(insert_pos) > 0:
        print('\nInsertions detected in the following positions in \"' + ref_genome + '\":\n')
        cur_idx = 0
        cur_ins = 0
        lag = 0
        print(str(insert_pos))

        for i in insert_pos:
            scanned = [site for list_sites in insertions.values() for site in list_sites]
            if i > cur_ins + 1: # not consecutive
                lag = len(scanned)
                print('gap openning - ' + str(i-lag))
                insertions[i] = [i-lag]
                cur_idx = i
            else: # is consecutive
                insertions[cur_idx].append(i-lag)
                print('\textension - ', str(i-lag))
            # print('')
            cur_ins = i
            print(insertions)

        print('\nThe sites above were recorded as insertions, and will be ignored in the next steps of mutation detection.')
    else:
        print('\nNo insertions were detected...\n')


    lines = []
    list_columns = ['sequence', 'gene', 'sense', 'position', 'mut_name', 'mut_type', 'ref_codon', 'var_codon',
                    'ref_aa', 'aa_pos', 'var_aa', 'grantham_score']
    sequences = SeqIO.parse(open(alignment), 'fasta')
    for entry in sequences:
        id, seq = entry.id, str(entry.seq)

        mut_name = ''
        # search for insertions
        list_seq = []
        lag = 0
        for ins in insert_pos:
            if seq[ins-1] != '-':
                if ins in insertions:
                    # insert_seq = seq[insertions[ins][0]-1:insertions[ins][-1]]
                    insert_seq = seq[ins-1: ins + len(insertions[ins])-1]
                    mut_name = 'c.' + str(insertions[ins][0]) + '_' + str(insertions[ins][0]+1) + 'ins' + insert_seq
                    # print(ins, mut_name, insertions[ins], insert_seq)
                    line = [id, '.', '.', str(ins), mut_name, 'insertion', '.', '.', '.', '.', '.', '']
                    print('\t'.join(line))
                    lines.append(line)

        for num, base in enumerate(seq, 1):
            if num not in insert_pos:
                # print(num, base)
                list_seq.append(base)
        seq = ''.join(list_seq) # genomes, with insertions removed

        dict_row = {}
        for trait in list_traits:
            pos = trait.split('|')[-1]
            start = int(pos.split(":")[0]) - 1
            end = int(pos.split(":")[1])
            sense = trait.split('|')[1]
            orf = str(seq)[start:end]

            if sense == '-': # reverse complement
                start, end = start - 1, end - 1
                orf = str(Seq(orf, generic_dna).reverse_complement())
            codons = [orf[i:i + 3] for i in range(0, len(orf), 3)]

            dict_row['genomes'] = id
            dict_row[trait] = codons
        dfCodons = dfCodons.append(dict_row, ignore_index=True)

    # put reference genome as first row of dataframe
    dfCodons["new"] = range(1, len(dfCodons) + 1)
    idx_refgenome = dfCodons[dfCodons['genomes'] == ref_genome].index.values.astype(int)[0]
    dfCodons.loc[idx_refgenome, 'new'] = 0
    dfCodons = dfCodons.sort_values("new")
    dfCodons = dfCodons.drop('new', axis=1)
    dfCodons.reset_index(drop=True, inplace=True)


    """ VARIANTS """
    def get_aa(seq):
        seq = Seq(seq, Gapped(IUPAC.unambiguous_dna))
        seq = seq.translate(table=1, to_stop=False)
        return str(seq)

    # generate dictionary from Grantham matrix
    matrix = matlist.grant
    grantham_matrix = {}
    for subs, score in matrix.items():
        grantham_score = 215 - score # covert from Vogt's score to Grantham's
        aa1, aa2 = subs
        if (aa1, aa2) not in grantham_matrix:
            grantham_matrix[str(aa1 + '_' + aa2)] = grantham_score
        if (aa2, aa1) not in grantham_matrix:
            grantham_matrix[str(aa2 + '_' + aa1)] = grantham_score


    # find variants
    exon_end = {}
    for column in list_traits:
        orf, sense, coord = column.split('|')
        start, end = [int(num) for num in coord.split(":")]
        if orf not in exon_end:
            exon_end[orf] = 0
        print('\n*', orf, sense, coord)
        l = dfCodons[column].tolist()
        list_codons = list(map(list, zip(*l)))
        for pos, codon_set in enumerate(list_codons):
            # print(pos, codon_set)
            aa_pos = pos + 1 + exon_end[orf]
            if sense == '-':
                nuc_pos = end - (pos * 3)
            else:
                nuc_pos = (((pos + 1) * 3) - 3) + start
            # print(orf, sense, aa_pos, nuc_pos, codon_set)

            start_nucpos, end_nucpos = nuc_pos, nuc_pos + 2
            if sense == '-':
                end_nucpos = nuc_pos - 2
            start_nucpos, end_nucpos = min([start_nucpos, end_nucpos]), max([start_nucpos, end_nucpos])

            variant_count = set(codon_set)
            if len(variant_count) > 1: # check if codons is strictly conserved
                # print(variant_count)
                ref_codon = codon_set[0] # reference codon
                for index, codon in enumerate(codon_set[1:], 1):
                    # print(index, codon)
                    genome = dfCodons.loc[index, 'genomes']
                    mut_type = ''
                    g_score = ''
                    atcg = set(codon).issubset('ATCG') # count numbers of ATCGs
                    if atcg:
                        ref_aa = get_aa(ref_codon)
                        var_aa = get_aa(codon)
                        if ref_codon != codon: # check if codons are identical
                            if '*' in var_aa: # stop codon
                                mut_name = orf + ':' + ref_aa + str(aa_pos) + var_aa
                                mut_type = 'nonsense'
                            else:
                                if ref_aa != var_aa:
                                    mut_name = orf + ':' + ref_aa + str(aa_pos) + var_aa
                                    if '*' in ref_aa:
                                        mut_type = 'nonstop'
                                    else:
                                        mut_type = 'nonsynonymous'
                                        g_score = grantham_matrix[str(ref_aa + '_' + var_aa)]
                                else:
                                    mut_name = orf + ':' + str(start_nucpos) + ref_codon + '>' + codon
                                    mut_type = 'synonymous'
                                    g_score = 0
                            line = [genome, orf, sense, str(start_nucpos) + '..' + str(end_nucpos), mut_name,
                                    mut_type, ref_codon, codon, ref_aa, str(aa_pos), var_aa, str(g_score)]
                            if len(export_only) > 0:
                                if mut_type in export_only:
                                    lines.append(line)
                                    print('\t'.join(line))
                            else:
                                if mut_type not in ignore:
                                    lines.append(line)
                                    print('\t'.join(line))

                    else:
                        ref_aa = get_aa(ref_codon)
                        var_aa = '.'
                        def gap_counter(string):
                            return string.count('-')
                        # print(ref_codon, gap_counter(ref_codon), codon, gap_counter(codon))
                        if gap_counter(ref_codon) == 0 and gap_counter(codon) == 3:
                            var_aa = '-'
                            mut_name = 'p.' + ref_aa + str(aa_pos) + 'del'
                            mut_type = 'deletion'
                        elif gap_counter(ref_codon) == 0 and 0 < gap_counter(codon) < 3:
                            mut_name = 'p.' + ref_aa + str(aa_pos) + 'fsX'
                            mut_type = 'frameshift'
                        elif gap_counter(codon) == 0 and atcg == False:
                            var_aa = get_aa(codon)
                            mut_name = 'p.' + ref_aa + str(aa_pos) + 'X'
                            mut_type = 'ambiguity'
                        else:
                            mut_name = 'x.NA'
                            mut_type = 'unknown'
                        line = [genome, orf, sense, str(start_nucpos) + '..' + str(end_nucpos), mut_name,
                                mut_type, ref_codon, codon, ref_aa, str(aa_pos), var_aa, str(g_score)]

                        if len(export_only) > 0:
                            if mut_type in export_only:
                                lines.append(line)
                                print('\t'.join(line))
                        else:
                            if mut_type not in ignore:
                                lines.append(line)
                                print('\t'.join(line))

        print('')
        exon_end[orf] += len(list_codons)
        # print(orf, exon_end[orf])

    table = pd.DataFrame(lines, columns=list_columns)
    table = table.sort_values(by=['sequence', 'mut_type'])
    # print(table)
    table.to_csv(output, sep='\t', index=False)
