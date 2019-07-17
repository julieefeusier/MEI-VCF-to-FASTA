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
                    help='MELT HG19 human transposon names (Alu, LINE1, or SVA)')

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

                ###This parses the vcf file
                vcf_split = line.split('\t')
		info_split = vcf_split[7].split(';')
                diff_info2 = info_split[5]
                MEI_info = diff_info2.split(",")

                ###The INFO column is different when DIFF=null. This corrects for these instances
                if re.search("DIFF=null", line):
                        diff_info = info_split[5].split(':')
                        remove_dif = diff_info[1].split(':')
                        diff_subs = remove_dif[0].split(',')
                        MEI_info[0] = diff_info[0].split(':')
                        MEI_info[1] = diff_subs[-3]
                        MEI_info[2] = diff_subs[-2]
                else:
                        diff_info = info_split[6].split(':')
                        remove_dif = diff_info[1].split(':')
                        diff_subs = remove_dif[0].split(',')
		
                ###Replacement of SNVs and deletions
		for i in diff_subs:
                        if (i[0] == 'g') or (i[0] == 'a') or (i[0] == 'c') or (i[0] == 't'):
				sub, location = sub_location(i)
				new_location = location - 1
						
				MEI_contig[new_location] = sub
			elif (i[0] == 'd'):
		                drange = deletion_range(i)
				for n in drange:
					m = n-1
					MEI_contig[m] = '-'
                        
		###This section adds dashes for 5' truncated loci	
                start = MEI_info[1]
                if start > '1':
                        truncation = "d1-" + str(int(start)-1)
                        drange = deletion_range(truncation)
                        for n in drange:
                                m = n-1
                                MEI_contig[m] = '-'
                        
                ###This section adds dashes for 3' truncated loci
                stop = MEI_info[2]
                if str(stop) != str(MEI_end):
                        truncation = "d" + str(int(stop)+1) + "-" + str(MEI_end)
                        drange = deletion_range(truncation)
                        for n in drange:
                                m = n-1
                                MEI_contig[m] = '-'			
                
                ###Replacement of inserted loci
                for i in diff_subs:
                        if i[0] == 'i':
                                 in_name = insertion_location(i)
                                 MEI_contig.insert(int(in_name[1])-1,in_name[2])

                ###Remove dashes from sequence
                while '-' in MEI_contig: MEI_contig.remove('-')
                
                ###Print the fasta sequence
                fasta.write(str('>' + vcf_split[0] + '_' + vcf_split[1] + '_' + info_split[5] + '_' + info_split[6] + '\n'))
		fasta.write(str(''.join(MEI_contig) + '\n'))
vcf_file.close()
