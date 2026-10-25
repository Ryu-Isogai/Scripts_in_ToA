#!/usr/bin/env python3
# coding:utf_8

import argparse
import random
import math
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description='Lower the coverage of the Rhodelphis sequence in a phylogenomic alignment by random selection.')
    parser.add_argument('-f', '--fasta_file', required=True, type=str, help='input FASTA file')
    parser.add_argument('-r', '--select_rate', required=True, type=int, help='percentage of valid positions to retain (0-100)')
    parser.add_argument('-n', '--trial_num', required=True, type=int, help='number of replicates to generate (>= 1)')
    parser.add_argument('-s', '--seed', required=False, default=100, type=int, help='seed value for the random sampling')
    parser.add_argument('-o', '--output_prefix', required=True, type=str, help='prefix of the output files')
    return parser.parse_args()

args=get_args()

# read input fasta file
input_fasta = {}

with open(args.fasta_file, mode='r') as fi:
    for record in SeqIO.parse(fi, 'fasta'):
        input_fasta[record.id] = record.seq
        
# name of target
Rhodelphis1 = "RhodelphisLimneticus"
Rhodelphis2 = "RhodelphisMarinus"

Rhodel1_seq = input_fasta[Rhodelphis1]
Rhodel2_seq = input_fasta[Rhodelphis2]


# get index list
def valid_aa_index(seq):
    index_list = []
    # valid amino acids list
    amino_acids_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    for i in range(len(seq)):
        if seq[i].upper() in amino_acids_list:
            index_list.append(i)
    return index_list

Rhodel1_valid_aa_index = valid_aa_index(Rhodel1_seq)
Rhodel2_valid_aa_index = valid_aa_index(Rhodel2_seq)


def extract_seq(seq, select_index_list):
    # make new seq corresponding to `select_index_list`
    new_seq = ''
    for i in range(len(seq)):
        if i in select_index_list:
            new_seq += seq[i]
        else:
            new_seq += '-'
    
    return new_seq


# number of selecting positions
R1_select_num = math.ceil(len(Rhodel1_valid_aa_index) * args.select_rate /100)
R2_select_num = math.ceil(len(Rhodel2_valid_aa_index) * args.select_rate /100)

# make low coverage sequences
for i in range(args.trial_num):
    # set seed
    current_seed = args.seed + i
    random.seed(current_seed)
    print(f'Trial {i+1} (seed:{current_seed})')
    
    # ramdom sampling
    R1_select_index = sorted(random.sample(Rhodel1_valid_aa_index,R1_select_num))
    R1_reduced_seq = extract_seq(Rhodel1_seq, R1_select_index)
    
    R2_select_index = sorted(random.sample(Rhodel2_valid_aa_index,R2_select_num))
    R2_reduced_seq = extract_seq(Rhodel2_seq, R2_select_index)
    
    # output files
    output_file_prefix = f'{args.output_prefix}_RLowCov{args.select_rate}_s{current_seed}'
    with open(output_file_prefix + '.fasta', mode='w') as fo:
        fo.write(f'>{Rhodelphis1}_select{args.select_rate}\n' + R1_reduced_seq + '\n')
        fo.write(f'>{Rhodelphis2}_select{args.select_rate}\n' + R2_reduced_seq)
    with open(output_file_prefix + '_R1.log', mode='w') as fo:
        fo.write('\n'.join([str(x+1) for x in R1_select_index]))
    with open(output_file_prefix + '_R2.log', mode='w') as fo:
        fo.write('\n'.join([str(x+1) for x in R2_select_index]))

print('\nComplete jobs!!\n')
