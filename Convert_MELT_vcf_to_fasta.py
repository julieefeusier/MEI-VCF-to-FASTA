#!/usr/bin/python
import re

import sys
from argparse import ArgumentParser

parser = ArgumentParser()
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-i', '--in',
                    metavar='INPUT VCF',
                    dest="i",
                    required=True,
                    help='name/path to VCF')

parser.add_argument('-o', '--out',
                    metavar='OUTPUT FASTA',
                    dest="o",
                    help='output FASTA file. Default will have the same prefix as vcf')

parser.add_argument('-m',
                    metavar='MELT HG19 TRANSPOSON',
                    dest="m",
                    type=str,
                    help='MELT Hg19/Hg38 human transposon families (Alu, LINE1, or SVA)')

parser.add_argument('-f', '--fasta',
                    metavar='TRANSPOSON FASTA',
                    dest="f",
                    help='transposon fasta used for MELT analysis')

args = parser.parse_args()

if args.i is None:
        raise NameError('Must include name/path to MELT VCF file with option -i')
else:
        input_file = args.i

if args.o is None:
        output_tmp = args.i
        output_file = output_tmp.rsplit(".",1)[0] + ".fasta"
else:
        output_file= args.o
fasta=open(output_file, 'w+')

if args.m is None and args.f is None:
        raise NameError('Must include either hg19 transposon name (Alu, LINE1, or SVA) with option -m or other MELT transposon file with option -f')
elif args.m is None and args.f is not None:
        if (args.f == "Alu") or (args.f == 'LINE1') or (args.f == 'SVA'):
                raise NameError('Please input a fasta file that is not named Alu, LINE1, or SVA. Or specify the transposon using option -m')
        else:
                MEI = args.f
elif args.m is not None and args.f is None:
        if (args.m == "Alu") or (args.m == 'LINE1') or (args.m == 'SVA'):
                MEI = args.m
        else:
                raise NameError('Must include hg19 transposon name (Alu, LINE1, or SVA) with option -m')

###These are the hg19 transposon FASTA sequences as providied by MELTv2.1.4
AluY= list('ggccgggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggcggatcacgaggtcaggagatcgagaccatcctggctaacacggtgaaaccccgtctctactaaaaatacaaaaaattagccgggcgtggtggcgggcgcctgtagtcccagctactcgggaggctgaggcaggagaatggcgtgaacccgggaggcggagcttgcagtgagccgagatcgcgccactgcactccagcctgggcgacagagcgagactccgtctca')

SVA = list('ctccctctccctcaccctctccccatggtctccctctccctctctttagtctcgttcactcagtgctcaatgatagcagcctgccttggcctcccaaagtgccgagattgcagcctctgcccggccgccaccccgtctgggaagtgaggagtgtctccgcctggccacccatcgtctgggatgtgaggagcgtctctgccctgccgcccatcgtctgagatgtggggagcacctctgcccggccgccccgtccgggatgtgaggagcgtcgctgcccggccgccccgtctgagaagtgaggagaccctctgcctggcaaccgctccatctgagaagtgaggagcccctccgcccggcagccgccctgtctgagaagtgaggagcccctccgcccagcagccacctggtccgggagggaggtgggggggtcagccccccgcccggccagccgccccgtccgggagggaggtgggggggtcagcccccagcccggccagccgccccgtccgggaagtgaggggcgcctctgcccggccgcccctactgggaagtgaggagccactttgcccggccagccactctgtccgggagggaggtgggggggtcagccccccgcccggccagccgccccgtccgggagggaggtggggggatcagccccccgcccagccagccgccccgtccgggagggaggtgggggggtcagccccccgcccggccagccgccctgtccgggaggtgaggggcgcctctgcccggccgcgcctactggaaagtgaggagcccctctgcccggccaccaccccgtctgggaggtgtgcccaacagctcattgagaaggggccatgatgacaatggcggttttgtggaatagaaaggggggaaaggtggggaaaagattgagaaatcggatggttgccgtgtctgtgtagaaagaggtagacctgggagacttttcattttgttctgtactaagaaaaattcttctgccttgggatcctgttgatcggtgaccttacccccaaccctgtgctctctgaaacatgtgctgtatccactcagggttgaatggattaagagcggtgcaagatgtgctttgttaaacagatgcttgaaggcagcatgctccttaagagtcatcaccactccctaatctcaagtacccagggacacaaacactgcggaaggccgcagggtcctctgcctaggaaaaccagagacctttgttcacttgtttatctgctgaccttccctccactattgtcctgtgaccctgccaaatccccctctgtgagaaacacccaagaatgatcaat')

LINE1 = list('gggggaggagccaagatggccgaataggaacagctccggtctacagctcccagcgtgagcgacgcagaagacggtgatttctgcatttccatctgaggtaccgggttcatctcactagggagtgccagacagtgggcgcaggccagtgtgtgtgcgcaccgtgcgcgagccgaagcagggcgaggcattgcctcacctgggaagcgcaaggggtcagggagttccctttctgagtcaaagaaaggggtgacggtcgcacctggaaaatcgggtcactcccacccgaatattgcgcttttcagaccggcttaagaaacggcgcaccacgagactatatcccacacctggctcggagggtcctacgcccacggaatctcgctgattgctagcacagcagtctgagatcaaactgcaaggcggcaacgaggctgggggaggggcgcccgccattgcccaggcttgcttaggtaaacaaagcagccgggaagctcgaactgggtggagcccaccacagctcaaggaggcctgcctgcctctgtaggctccacctctgggggcagggcacagacaaacaaaaagacagcagtaacctctgcagacttaagtgtccctgtctgacagctttgaagagagcagtggttctcccagcacgcagctggagatctgagaacgggcagacagactgcctcctcaagtgggtccctgactcctgacccccgagcagcctaactgggaggcaccccccagcaggggcacactgacacctcacacggcagggtattccaacagacctgcagctgagggtcctgtctgttagaaggaaaactaacaaccagaaaggacatctacaccgaaaacccatctgtacatcaccatcatcaaagaccaaaagtagataaaaccacaaagatggggaaaaaacagaacagaaaaactggaaactctaaaacgcagagcgcctctcctcctccaaaggaacgcagttcctcaccagcaacggaacaaagctggatggagaatgattttgacgagctgagagaagaaggcttcagacgatcaaattactctgagctacgggaggacattcaaaccaaaggcaaagaagttgaaaactttgaaaaaaatttagaagaatgtataactagaataaccaatacagagaagtgcttaaaggagctgatggagctgaaaaccaaggctcgagaactacgtgaagaatgcagaagcctcaggagccgatgcgatcaactggaagaaagggtatcagcaatggaagatgaaatgaatgaaatgaagcgagaagggaagtttagagaaaaaagaataaaaagaaatgagcaaagcctccaagaaatatgggactatgtgaaaagaccaaatctacgtctgattggtgtacctgaaagtgatgtggagaatggaaccaagttggaaaacactctgcaggatattatccaggagaacttccccaatctagcaaggcaggccaacgttcagattcaggaaatacagagaacgccacaaagatactcctcgagaagagcaactccaagacacataattgtcagattcaccaaagttgaaatgaaggaaaaaatgttaagggcagccagagagaaaggtcgggttaccctcaaaggaaagcccatcagactaacagtggatctctcggcagaaaccctacaagccagaagagagtgggggccaatattcaacattcttaaagaaaagaattttcaacccagaatttcatatccagccaaactaagcttcataagtgaaggagaaataaaatactttatagacaagcaaatgttgagagattttgtcaccaccaggcctgccctaaaagagctcctgaaggaagcgctaaacatggaaaggaacaaccggtaccagccgctgcaaaatcatgccaaaatgtaaagaccatcgagactaggaagaaactgcatcaactaatgagcaaaatcaccagctaacatcataatgacaggatcaaattcacacataacaatattaactttaaatataaatggactaaattctgcaattaaaagacacagactggcaagttggataaagagtcaagacccatcagtgtgctgtattcaggaaacccatctcacgtgcagagacacacataggctcaaaataaaaggatggaggaagatctaccaagccaatggaaaacaaaaaaaggcaggggttgcaatcctagtctctgataaaacagactttaaaccaacaaagatcaaaagagacaaagaaggccattacataatggtaaagggatcaattcaacaagaggagctaactatcctaaatatttatgcacccaatacaggagcacccagattcataaagcaagtcctcagtgacctacaaagagacttagactcccacacattaataatgggagactttaacaccccactgtcaacattagacagatcaacgagacagaaagtcaacaaggatacccaggaattgaactcagctctgcaccaagcagacctaatagacatctacagaactctccaccccaaatcaacagaatatacctttttttcagcaccacaccacacctattccaaaattgaccacatagttggaagtaaagctctcctcagcaaatgtaaaagaacagaaattataacaaactatctctcagaccacagtgcaatcaaactagaactcaggattaagaatctcactcaaagccgctcaactacatggaaactgaacaacctgctcctgaatgactactgggtacataacgaaatgaaggcagaaataaagatgttctttgaaaccaacgagaacaaagacaccacataccagaatctctgggacgcattcaaagcagtgtgtagagggaaatttatagcactaaatgcctacaagagaaagcaggaaagatccaaaattgacaccctaacatcacaattaaaagaactagaaaagcaagagcaaacacattcaaaagctagcagaaggcaagaaataactaaaatcagagcagaactgaaggaaatagagacacaaaaaacccttcaaaaaatcaatgaatccaggagctggttttttgaaaggatcaacaaaattgatagaccgctagcaagactaataaagaaaaaaagagagaagaatcaaatagacacaataaaaaatgataaaggggatatcaccaccgatcccacagaaatacaaactaccatcagagaatactacaaacacctctacgcaaataaactagaaaatctagaagaaatggatacattcctcgacacatacactctcccaagactaaaccaggaagaagttgaatctctgaatagaccaataacaggctctgaaattgtggcaataatcaatagtttaccaaccaaaaagagtccaggaccagatggattcacagccgaattctaccagaggtacatggaggaactggtaccattccttctgaaactattccaatcaatagaaaaagagggaatcctccctaactcattttatgaggccagcatcattctgataccaaagccgggcagagacacaaccaaaaaagagaattttagaccaatatccttgatgaacattgatgcaaaaatcctcaataaaatactggcaaaccgaatccagcagcacatcaaaaagcttatccaccatgatcaagtgggcttcatccctgggatgcaaggctggttcaatatacgcaaatcaataaatgtaatccagcatataaacagagccaaagacaaaaaccacatgattatctcaatagatgcagaaaaagcctttgacaaaattcaacaacccttcatgctaaaaactctcaataaattaggtattgatgggacgtatttcaaaataataagagctatctatgacaaacccacagccaatatcatactgaatgggcaaaaactggaagcattccctttgaaaaccggcacaagacagggatgccctctctcaccgctcctattcaacatagtgttggaagttctggccagggcaatcaggcaggagaaggaaataaagggtattcaattaggaaaagaggaagtcaaattgtccctgtttgcagacgacatgattgtatatctagaaaaccccatcgtctcagcccaaaatctccttaagctgataagcaacttcagcaaagtctcaggatacaaaatcaatgtacaaaaatcacaagcattcttatacaccaacaacagacaaacagagagccaaatcatgggtgaactcccattcgtaattgcttcaaagagaataaaatacctaggaatccaacttacaagggatgtgaaggacctcttcaaggagaactacaaaccactgctcaaggaaataaaagaggacacaaacaaatggaagaacattccatgctcatgggtaggaagaatcaatatcgtgaaaatggccatactgcccaaggtaatttacagattcaatgccatccccatcaagctaccaatgactttcttcacagaattggaaaaaactactttaaagttcatatggaaccaaaaaagagcccgcattgccaagtcaatcctaagccaaaagaacaaagctggaggcatcacactacctgacttcaaactatactacaaggctacagtaaccaaaacagcatggtactggtaccaaaacagagatatagatcaatggaacagaacagagccctcagaaataatgccgcatatctacaactatctgatctttgacaaacctgagaaaaacaagcaatggggaaaggattccctatttaataaatggtgctgggaaaactggctagccatatgtagaaagctgaaactggatcccttccttacaccttatacaaaaatcaattcaagatggattaaagatttaaacgttaaacctaaaaccataaaaaccctagaagaaaacctaggcattaccattcaggacataggcgtgggcaaggacttcatgtccaaaacaccaaaagcaatggcaacaaaagacaaaattgacaaatgggatctaattaaactaaagagcttctgcacagcaaaagaaactaccatcagagtgaacaggcaacctacaacatgggagaaaattttcgcaacctactcatctgacaaagggctaatatccagaatctacaatgaacttaaacaaatttacaagaaaaaaacaaacaaccccatcaaaaagtgggcgaaggacatgaacagacacttctcaaaagaagacatttatgcagccaaaaaacacatgaagaaatgctcatcatcactggccatcagagaaatgcaaatcaaaaccactatgagatatcatctcacaccagttagaatggcaatcattaaaaagtcaggaaacaacaggtgctggagaggatgcggagaaataggaacacttttacactgttggtgggactgtaaactagttcaaccattgtggaagtcagtgtggcgattcctcagggatctagaactagaaataccatttgacccagccatcccattactgggtatatacccaaatgagtataaatcatgctgctataaagacacatgcacacgtatgtttattgcggcactattcacaatagcaaagacttggaaccaacccaaatgtccaacaatgatagactggattaagaaaatgtggcacatatacaccatggaatactatgcagccataaaaaatgatgagttcatatcctttgtagggacatggatgaaattggaaaccatcattctcagtaaactatcgcaagaacaaaaaaccaaacaccgcatattctcactcataggtgggaattgaacaatgagatcacatggacacaggaaggggaatatcacactctggggactgtggtggggtcgggggaggggggagggatagcattgggagatatacctaatgctagatgacacattagtgggtgcagcgcaccagcatggcacatgtatacatatgtaactaacctgcacaatgtgcacatgtaccctaaaacttagagtat')

Alu_subfamilies = {}
Alu_subfamilies = {
    'AluYa1_1': ['g145a'],
    'AluYa1_2': ['g237c'],
    'AluYa2': ['g145a', 'c174t'],
    'AluYa3_1': ['t89c', 'c96a', 'g145a'],
    'AluYa3_2': ['t89c', 'c174t', 'g237c'],
    'AluYa3_3': ['t89c', 'c96a', 'g237c'],
    'AluYa3_4': ['g145a', 'c174t', 'g237c'],
    'AluYa3_5': ['c96a', 'g145a', 'c174t'],
    'AluYa4_1': ['t89c', 'g145a', 'c174t', 'g237c'],
    'AluYa4_2': ['t89c', 'c96a', 'g145a', 'g237c'],
    'AluYa4_3': ['t89c', 'c96a', 'c174t', 'g237c'],
    'AluYa4_4': ['c96a', 'g145a', 'c174t', 'g237c'],
    'AluYa4_5': ['t89c', 'c96a', 'g145a', 'c174t'],
    'AluYa5a1': ['a82g', 't89c', 'c96a', 'g145a', 'c174t', 'g237c'],
    'AluYa5a2': ['g72a', 't89c', 'c96a', 'i133a', 'g145a', 'c174t', 'g237c'],
    'AluYa5b1': ['t89c', 'c96a', 'g145a', 'a161g', 'c174t', 'g237c'],
    'AluYa5b2': ['t86a', 't89c', 'g90c', 'c96a', 'g145a', 'c174t', 'g237c'],
    'AluYa5c1': ['t89c', 'c96a', 'g145a', 'c174t', 'a234g', 'g237c'],
    'AluYa5': ['t89c', 'c96a', 'g145a', 'c174t', 'g237c'],
    'AluYa8': ['t89c', 'c96a', 'a123c', 'd133', 'g145a', 'c166t', 'c174t', 'g237c'],
    'AluYb3a1': ['c57t', 'c236t', 'c248g'],
    'AluYb3a2': ['t144c', 'c236a', 'c248g', 'a252g', 'c263a'],
    'AluYb5': ['c57t', 'c64t', 'c98a', 't144c', 'c236t'],
    'AluYb6_1': ['c57t', 'c64t', 'c98a', 'c236t', 'c248g', 'i252gcagtcc', 'a252g'],
    'AluYb6_2': ['c57t', 'c64t', 'c98a', 't144c', 'g211a', 'c236t'],
    'AluYb7_1': ['c57t', 'c64t', 'c98a', 't144c', 'g211a', 'c248g', 'i252gcagtcc', 'a252g'],
    'AluYb7_2': ['c57t', 'c64t', 'c98a', 'g211a', 'c236t', 'c248g', 'i252gcagtcc', 'a252g'],
    'AluYb7_3': ['c57t', 'c98a', 't144c', 'g211a', 'c236t', 'c248g', 'i252gcagtcc', 'a252g'],
    'AluYb7_4': ['c57t', 'c64t', 'c98a', 't144c', 'c236t', 'c248g', 'i252gcagtcc', 'a252g'],
    'AluYb8': ['c57t', 'c64t', 'c98a', 't144c', 'g211a', 'c236t', 'c248g', 'i252gcagtcc', 'a252g'],
    'AluYb9': ['c57t', 'c64t', 'c98a', 't144c', 'c174g', 'g211a', 'c236t', 'c248g', 'i252gcagtcc', 'a252g'],
    'AluYb10': ['c57t', 'c64t', 'c98a', 't144c', 'c174g', 'g211a', 'c236t', 'c248g', 'i252gcagtcc'],
    'AluYb11': ['c57t', 'c64t', 'c98a', 't144c', 'c174g', 'i201t', 'g211a', 'c236t', 'c248g', 'i252gcagtcc'],
    'AluYbc3a': ['i133a', 'g148a', 'g231a', 'c236a', 'c248g', 'c263a'],
    'AluYc1': ['g148a'],
    'AluYc2': ['c98a', 'g148a'],
    'AluYc5': ['d86-97', 't144c', 'g145a', 'g209a', 'c230g', 't244a', 'a252c'],
    'AluYd': ['d86-97'],
    'AluYd2': ['d86-97', 't144c'],
    'AluYd3a1': ['c23t', 'd86-97', 't144c', 'i177a'],
    'AluYd3': ['c23t', 'd86-97', 't144c'],
    'AluYd8': ['d86-97', 't144c', 'g145a', 'g209a', 'c230g', 't244a', 'a252c', 'g267a'],
    'AluYe2': ['t144a', 'd265-266'],
    'AluYe4': ['t144a', 'i207c', 'a210g', 't219c'],
    'AluYe5': ['t144a', 'i207c', 'a210g', 't219c', 'd265-266'],
    'AluYe6': ['t144a', 'i207c', 'a210g', 't219c', 'g253a', 'd265-266'],
    'AluYf1': ['t173g'],
    'AluYf2': ['g99a', 'i124a', 't124c', 'i133a', 't173g'],
    'AluYg5b3': ['g52a', 'a94g', 'g143a', 'g152c', 'c229t', 'a247g', 'g271a'],
    'AluYg6a2': ['g52a', 'g143a', 'g152c', 'c154t', 't173a', 'g175a', 'c229t', 'g271a'],
    'AluYg6': ['g52a', 'g143a', 'g152c', 't173a', 'c229t', 'g271a'],
    'AluYh3': ['c48a', 'c64t', 't111c'],
    'AluYh7': ['c48a', 'c64t', 'a97g', 't111c', 'a161g', 'a167g', 'a234g'],
    'AluYh9': ['c48a', 'c64t', 'a97g', 't111c', 'a161g', 'a167g', 'c230t', 'a234g', 't249c'],
    'AluYi6': ['c23t', 'i133a', 'g146a', 'c236t', 'a252c', 'g261c'],
    'AluYj3': ['a125g', 'a171g', 'c238t'],
    'AluYj4': ['a19g', 'a125g', 'a171g', 'c238t'],
    'AluYk13': ['c10t', 'c57t', 'c64t', 'g99a', 'c142a', 't144c', 'g155t', 'c174t', 'g176a', 'g207a', 'g269c', 't274g', 'c276t']
}


###This function finds the correct nucleotide position to fix
def sub_location(str_i):
    value_i = ''
    length_i = ''
    if len(str_i) == 3:
        value_i = str_i[2]
        length_i = str_i[1]
    elif len(str_i) == 4:
        value_i = str_i[3]
        length_i = str_i[1:3]
    elif len(str_i) == 5:
        value_i = str_i[4]
        length_i = str_i[1:4]
    elif len(str_i) == 6:
        value_i = str_i[5]
        length_i = str_i[1:5]
    return value_i, int(length_i)

###This function replaces the consensus region with a dash for a deleted region
def deletion_range(x):
    x = x[1:]
    result = []
    for part in x.split(','):
        if '-' in part:
            a, b = part.split('-')
            a, b = int(a), int(b)
            result.extend(range(a, b + 1))
        else:
            a = int(part)
            result.append(a)
    return result

###This function adds inserted sequences into the consensus region    
def insertion_location(str_i):
    insert_name = re.split('(\d+)',str_i)
    return insert_name

###This function sorts insertions reverse numerically
def insertion_sort(list_a):
    b = {}
    b_sorted = []
    for i in list_a:
        if (i[1].isnumeric() == False):
            b[i[0]] = i[1:]
        elif (i[2].isnumeric() == False):
            b[i[0:2]] = i[2:]
        elif (i[3].isnumeric() == False):
            b[i[0:3]] = i[3:]
        elif (i[4].isnumeric() == False):
            b[i[0:4]] = i[4:]
    for key, value in sorted(b.items(), key=lambda item: int(item[0]), reverse=True):
        b_sorted.append(key + value)
    return(b_sorted)

###Takes the TE input for the consensus sequence
MEI_contig = []
if MEI == 'Alu':
        MEI_contig = AluY[:]
        MEI_end = 281
elif MEI == 'LINE1':
        MEI_contig = LINE1[:]
        MEI_end = 6019
elif MEI == 'SVA':
        MEI_contig = SVA[:]
        MEI_end = 1316
elif MEI == args.f:
        transposon_line = ""
        ###This extracts the single- or multi-line FASTA
        with open(MEI,'r') as transposon_file:
                for t_line in transposon_file:
                        if re.search(">",t_line):
                                continue
                        transposon_line += t_line
        transposon_file.close()
        transposon_line1 = transposon_line.replace('\n','')
        MEI_contig = list(transposon_line1.lower())
        MEI_end = len(transposon_line1)-1

MEI_contig_original = MEI_contig[:]

with open(input_file,'r') as vcf_file:
    for line in vcf_file:
        MEI_contig = MEI_contig_original[:]
                
        ###This skips all header lines in the vcf
        if re.search("#", line):
            continue

        ###This skips Alu elements that are classified as SVA elements
        if MEI == 'Alu':
            if re.search("MEINFO=SVA", line):
                continue

        ###This parses the vcf file
        vcf_split = line.split('\t')
        info_split = vcf_split[7].split(';')
        diff_info2 = info_split[5]
        MEI_info = diff_info2.split(",")
        MEI_sub = MEI_info[0].split('=') ##
        Alu_sub = MEI_sub[1]

        if re.search("DIFF=null", line):
            diff_info = info_split[5].split(':')
            diff_subs = []

        else:
            diff_info = info_split[6].split(':')
            remove_dif = diff_info[1].split(':')
            diff_subs = remove_dif[0].split(',')

        ###Add Subfamily changes for MELTv2.2.0
        if Alu_sub in Alu_subfamilies:
            for i in Alu_subfamilies[Alu_sub]:
                diff_subs.append(i)    
        #remove duplicated mutations
        diff_subs = list( dict.fromkeys(diff_subs) )
        diff_ind = []
        diff_del = []
        for i in diff_subs:
            if i[0] =='i':
                diff_ind.append(i)
            if i[0] == 'd':
                diff_del.append(i)
        
        ###Replacement of SNVs and deletions
        for i in diff_subs:
            if (i[0] == 'g') or (i[0] == 'a') or (i[0] == 'c') or (i[0] == 't'):
                sub, location = sub_location(i)
                new_location = location - 1
                MEI_contig[new_location] = sub
        for i in diff_del:
            drange = deletion_range(i)
            for n in drange:
                m = n-1
                MEI_contig[m] = '-'
                        
        ###This section adds dashes for 5' truncated loci    
        start = MEI_info[1]
        if str(start) > str(MEI_end):
            start = MEI_end
        if start > '1':
            truncation = "d1-" + str(int(start)-1)
            drange = deletion_range(truncation)
            for n in drange:
                m = n-1
                MEI_contig[m] = '-'
                        
        ###This section adds dashes for 3' truncated loci
        stop = MEI_info[2]
        if str(stop) > str(MEI_end):
            stop = MEI_end
        if str(stop) != str(MEI_end):
            truncation = "d" + str(int(stop)+1) + "-" + str(MEI_end)
            drange = deletion_range(truncation)
            for n in drange:
                m = n-1
                MEI_contig[m] = '-'            
                
        ###Replacement of inserted loci
        diff_ind = [e[1:] for e in diff_ind]
        diff_ind = insertion_sort(diff_ind)
        for i in diff_ind:
            in_name = insertion_location(i)
            MEI_contig.insert(int(in_name[1])-1,in_name[2])

        ###Remove dashes from sequence
        while '-' in MEI_contig: MEI_contig.remove('-')
                
        ###Print the fasta sequence
        fasta.write(str('>' + vcf_split[0] + '_' + vcf_split[1] + '_' + info_split[5] + '_' + info_split[6] + '\n'))
        fasta.write(str(''.join(MEI_contig) + '\n'))
vcf_file.close()
