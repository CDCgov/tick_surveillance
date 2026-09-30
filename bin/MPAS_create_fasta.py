#!/usr/bin/env python

import os
# from tkinter import TRUE
import pandas as pd
import csv
import sys


# setup command two line arguments:
# 1. sequencing_report.xlsx
# 2. targets.tsv

sequencing_report: sys.argv[1]
targets: sys.argv[2]


# 1. read in all_data tab and metadata tabs from sequening_report.xlsx file as two separate dataframes and read in the targets.tsv file

val = pd.read_excel(sys.argv[1], sheet_name='all_data')
meta = pd.read_excel(sys.argv[1], sheet_name='metadata')
test_result = pd.read_excel(sys.argv[1], sheet_name='TestingResults')
targets = pd.read_table(sys.argv[2])

# 2. clean up each dataframe to only include specified fields

#all_data DF: only indluce Index, abundance, sequence number, subject, percent identity, and sequence
seqs = val[['Index', 'abundance', 'primer_name', 'sequence_number', 'subject','percent_identity', 'observed_sequence']]

# remove control sequences
seqs = seqs[seqs['subject'].str.contains('control', case=False)==False]

#metadata DF: only include CSID, Index, Tick_Genus_species, Lifestage, state
meta_sub = meta[['Index', 'CSID', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'State']]

#results sub
if 'Molecular_Tick_ID' in test_result.columns:
    test_result_sub = test_result[['CSID', 'Molecular_Tick_ID']]
else:
    print('Molecular_Tick_ID not present, morphological ID will be used in tree')
    test_result_sub = test_result[['CSID']]

#targets DF:
targets_sub = targets[['primer_name', 'ref_sequence_name', 'sequence']]
targets_sub = targets_sub[targets_sub['ref_sequence_name'].str.contains('control', case=False)==False]

# 3. Join seq and meta_sub into one DF based on Index. Replace all spaces with underscore.
meta_sub_seqs = pd.merge(seqs, meta_sub, on='Index').replace(' ', '_', regex=True)

# 4. Keep all sequences with abundance greter than 50
seq_50 = meta_sub_seqs[meta_sub_seqs['abundance'] >=50]

# merge seq_50 and test_result_sub
seq_50 = seq_50.merge(test_result_sub, on='CSID', how='outer')

#5. Add new column that contains only the first 7 characters of the sequence_number

seq_50['new_seq_number'] = seq_50['sequence_number'].str[:7]
#seq_50 = seq_50.assign(new_seq_number = seq_50['sequence_number'].astype(str).str[:7])

# 6. Group representative sequences found in each state and tick species and create CSV - use molecular tick ID column if available.
#seq_50_group = seq_50.groupby(['State','sequence_number', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'primer_name', 'observed_sequence', 'new_seq_number'], dropna=False).size().reset_index(name='n')
#seq_50_group_reorder = seq_50_group[['State', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'new_seq_number', 'primer_name', 'n', 'observed_sequence']]
if 'Molecular_Tick_ID' in seq_50.columns:
    seq_50_group = seq_50.groupby(['State','sequence_number', 'Molecular_Tick_ID', 'Lifestage', 'primer_name', 'observed_sequence', 'new_seq_number'], dropna=False).size().reset_index(name='n')
    seq_50_group_reorder = seq_50_group[['State', 'Molecular_Tick_ID', 'Lifestage', 'new_seq_number', 'primer_name', 'n', 'observed_sequence']]
    seq_50_new = seq_50_group_reorder.rename(columns={'Molecular_Tick_ID': 'Tick_ID'})
else:
    seq_50_group = seq_50.groupby(['State','sequence_number', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'primer_name', 'observed_sequence', 'new_seq_number'], dropna=False).size().reset_index(name='n')
    seq_50_group_reorder = seq_50_group[['State', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'new_seq_number', 'primer_name', 'n', 'observed_sequence']]
    seq_50_new = seq_50_group_reorder.rename(columns={'Morphological_Ectoparasite_Genus_Species': 'Tick_ID'})

# 7. Funtions to Convert sample CSV and ref seq CSV files to fasta files

def write_fasta(in_seq, out_seq):

    try:
        fin = open(in_seq, 'r')
        fout = open(out_seq,'a')
    except:
        return -1

    with fin:
        reader = csv.DictReader(fin, skipinitialspace=True)
        for line in reader:
            state = line['State']
            genus = line['Tick_ID']
            lifestage = line['Lifestage']
            number = line['new_seq_number']
            tot = line['n']
            ob_seq = line['observed_sequence']
            print(f'>{state}_{genus}_{lifestage}_seq{number}_n={tot}\n{ob_seq}', file=fout)
    return

def write_ref_fasta(infile, outfile):
    try:
        f_in = open(infile, 'r')
        f_out = open(outfile,'a')
    except:
        return -1

    with f_in:
        read = csv.DictReader(f_in, skipinitialspace=True)
        for line in read:
            name = line['ref_sequence_name']
            ref_seq = line['sequence']
            print(f'>{name}\n{ref_seq}', file=f_out) 
    return

# 8.Create csv files and fasta files for samples and ref sequences based on primers

for primers, seqs in seq_50_new.groupby(['primer_name']):
    seqs.to_csv(f'{primers}.csv', index=False)
    write_fasta(f'{primers}.csv', f'{primers}_all.fasta')

for primer, seq in targets_sub.groupby(['primer_name']):
    seq.to_csv(f'{primer}_ref.csv', index=False)
    write_ref_fasta(f'{primer}_ref.csv', f'{primer}_all.fasta')