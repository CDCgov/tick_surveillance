#!/usr/bin/env python

'''
Python script that creates new target file based on iPM number and study type provided by user.
Created by: LMO
Date: 02/10/25

'''

# 1. load python modules
import pandas as pd
import numpy as npy
import csv
import os
import sys


# 2. load target excel tab
targets = pd.read_excel(sys.argv[1], sheet_name='targets')

# 3. load primers tab
primers = pd.read_excel(sys.argv[1], sheet_name='primers')

# 4. primer mix user input - This will be a param
pm_no = sys.argv[2]
split_pm = pm_no.split(",")

# 5. surveillance or other?
study_type_in = sys.argv[3]
#make all study_type input lower case
study_type = study_type_in.lower()
#check that study_type is surveillance or research
study_type_options = ('surveillance', 'research')

if study_type in study_type_options:
    pass
else:
    exit("Error: Incorrect study_type provided. The input --study_type must be surveillance or research. Exit process")

# 6. create df of primes names based on input primer_mix
select_primers = primers[primers['primer_mix'].isin(split_pm)]
select_primers_new = select_primers.drop(columns=['primer_mix']) #remove primer_mix column
uniq_pfile = select_primers_new.drop_duplicates() # primers for primer.tsv file, duplicates removed
#checks if input primer_mix name is found in the primers tab of refseq_all.xlsx
if len(uniq_pfile) > 0:
    pass
else:
    exit("Error: The provided --primer_mix cannnot be found in the the primers tab of reafseq_all.xlsx. Exit process.")


# 7. capture only the primer name
p_name = select_primers['primer_name']
unique_select_primers = list(set(p_name)) #remove duplicate primers

# 8. select only ref sequences that are for the selected iPM primer name
new_targets = targets[targets['primer_name'].isin(unique_select_primers)]
#checks if primer_name from targets tab is in primer_name of primers tab
if len(new_targets) > 0: 
    new_targets = pd.DataFrame(new_targets)
else: 
    exit('Error: The targets tab primer_name column cannot be found in primer_name column of primers tab in refseq_all.xlsx. Exit process.')

# 9. rename  columns requried in targets file
# if Genus_species_strain_Genbank_Tick_Host_Location_Primers, rename to 'reference_sequence_name'
rename_col = 'Genus_species_strain_Genbank_Tick_Host_Location_Primers'
new_refcol_name = 'ref_sequence_name'

if rename_col in new_targets.columns:
    new_targets.rename(columns={rename_col: new_refcol_name}, inplace=True)
else:
    pass

# 10. conditional based on user input if desired target file is for surveillance or 'other'. Update reporting column name
# if for surviellance, the target df should contain 'surveillance_report_col_format'
# if not surviellance, target df should countain 'reseach_report_col_format'

if study_type == 'surveillance':
    filter_targets = new_targets[['ref_sequence_name', 'species', 'primer_name', 'SURVEILLANCE_reporting_columns', 'SURVEILLANCE_min_percent_identity', 'SURVEILLANCE_min_percent_aligned', 'SURVEILLANCE_max_percent_gaps', 'SURVEILLANCE_internal_control', 'sequence']]
    filter_targets = filter_targets.rename(columns={'SURVEILLANCE_reporting_columns':'reporting_columns', 'SURVEILLANCE_min_percent_identity': 'min_percent_identity', 'SURVEILLANCE_min_percent_aligned':'min_percent_aligned', 'SURVEILLANCE_max_percent_gaps':'max_percent_gaps', 'SURVEILLANCE_internal_control':'internal_control'})

else:
    filter_targets = new_targets[['ref_sequence_name', 'species', 'primer_name', 'RESEARCH_reporting_columns', 'RESEARCH_min_percent_identity', 'RESEARCH_min_percent_aligned', 'RESEARCH_max_percent_gaps', 'RESEARCH_internal_control', 'sequence']]
    filter_targets = filter_targets.rename(columns={'RESEARCH_reporting_columns':'reporting_columns', 'RESEARCH_min_percent_identity': 'min_percent_identity', 'RESEARCH_min_percent_aligned':'min_percent_aligned', 'RESEARCH_max_percent_gaps':'max_percent_gaps', 'RESEARCH_internal_control':'internal_control'})
    

# 11. filter this df to only caputre columns required in the targets file
more_filter_target = filter_targets[['ref_sequence_name', 'species', 'primer_name', 'reporting_columns', 'min_percent_identity', 'min_percent_aligned', 'max_percent_gaps', 'internal_control', 'sequence']]

# 12. save tsv file
more_filter_target.to_csv(f"targets.tsv", sep='\t', index=False)

# 13. save primer file
primer_file = uniq_pfile[['primer_name', 'primer_f_name', 'primer_f_seq', 'primer_r_name', 'primer_r_seq']]
primer_file.to_csv(f"primers.tsv", sep='\t', index=False)
