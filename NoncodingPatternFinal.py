'''
Created on Jul 27, 2021

@author: daily
'''
import sys
import subprocess

# install dependency (pandas)
reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
if 'pandas' not in installed_packages:
    print('Pandas dependency not found. Attempting to install...')
    print()
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas'])
    reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
    installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
    if 'pandas' not in installed_packages:
        print('Installation failed. Please manually install the pandas library before continuing.')
        print()
        input("Press Enter to close...")
        exit()
    else:
        print('Installation successful. Continuing...')
        print()

import pandas as pd
import os.path
pd.set_option('display.max_columns', None)
pd.set_option('display.min_rows', 50)

# check if annotation file exists
if os.path.isfile('Annotations.gff3'):
    print("Parsing Annotation File...")
    print()
else:
    print("Could not find the 'Annotations.gff3' file. Has it been decompressed?")
    print()
    input("Press Enter to close...")
    exit()

# import annotation data and split based on scaffold
raw_data = pd.read_csv('Annotations.gff3', sep='\t', header=None, prefix='Column', skiprows=2)
split_data = dict()
grCols = ['Column0']
split_data = {i: y for i, (d, y) in enumerate(raw_data.groupby(grCols))}

# setup output files
stats_out = open('overall_statistics.txt', 'w')
print('Scaffold','ncRNA/mRNA Match', 'Total ncRNA', 'ncRNA Match/Total', 'Total mRNA', 'mRNA Match/Total', sep='\t', file=stats_out)
annotations_out = open('matching_annotations.gff3', 'w')

# parse annotation data and print to output files
for r in range(len(split_data)):
    match_count = 0
    split_data[r].reset_index(drop=True,inplace=True)
    
    if split_data[r].at[0, 'Column0'].find('#') < 0:
        all_ncRNA_data = split_data[r].query('Column2 == "ncRNA" | Column2 == "lnc_RNA"')
        all_ncRNA_data.reset_index(drop=True,inplace=True)
        all_ncRNA_data = all_ncRNA_data.astype({'Column3' : 'int64', 'Column4' : 'int64'})
        
        all_mRNA_data = split_data[r].query('Column2 == "mRNA"')
        all_mRNA_data.reset_index(drop=True,inplace=True)
        all_mRNA_data = all_mRNA_data.astype({'Column3' : 'int64', 'Column4' : 'int64'})
        
        filtered_data = pd.concat([split_data[r].query('Column2 == "mRNA"'), split_data[r].query('Column2 == "ncRNA" | Column2 == "lnc_RNA"')], axis=0)
        sorted_data = filtered_data.sort_values(by='Column3')
        sorted_data.reset_index(drop=True,inplace=True)
        sorted_data = sorted_data.astype({'Column3' : 'int64', 'Column4' : 'int64'})
            
        for i in range(len(sorted_data)):
            if (sorted_data.at[i, 'Column2'] == 'ncRNA' or sorted_data.at[i, 'Column2'] == 'lnc_RNA') and sorted_data.at[i, 'Column6'] == '-' and i < (len(sorted_data)-1):
                if sorted_data.at[i+1, 'Column2'] == 'mRNA' and sorted_data.at[i+1, 'Column6'] == '+':
                    match_count += 1
                    print(sorted_data.loc[i])
                    print(sorted_data.loc[i+1])
                    print(sorted_data.at[i, 'Column0'], sorted_data.at[i, 'Column1'], sorted_data.at[i, 'Column2'], sorted_data.at[i, 'Column3'], sorted_data.at[i, 'Column4'], sorted_data.at[i, 'Column5'], sorted_data.at[i, 'Column6'], sorted_data.at[i, 'Column7'], sorted_data.at[i, 'Column8'], sep='\t', file=annotations_out)
                    print(sorted_data.at[i+1, 'Column0'], sorted_data.at[i+1, 'Column1'], sorted_data.at[i+1, 'Column2'], sorted_data.at[i+1, 'Column3'], sorted_data.at[i+1, 'Column4'], sorted_data.at[i+1, 'Column5'], sorted_data.at[i+1, 'Column6'], sorted_data.at[i+1, 'Column7'], sorted_data.at[i+1, 'Column8'], sep='\t', file=annotations_out)
                    print()
            elif (sorted_data.at[i, 'Column2'] == 'ncRNA' or sorted_data.at[i, 'Column2'] == 'lnc_RNA') and sorted_data.at[i, 'Column6'] == '+' and i > 0:
                if sorted_data.at[i-1, 'Column2'] == 'mRNA' and sorted_data.at[i-1, 'Column6'] == '-':
                    match_count += 1
                    print(sorted_data.loc[i])
                    print(sorted_data.loc[i-1])
                    print(sorted_data.at[i, 'Column0'], sorted_data.at[i, 'Column1'], sorted_data.at[i, 'Column2'], sorted_data.at[i, 'Column3'], sorted_data.at[i, 'Column4'], sorted_data.at[i, 'Column5'], sorted_data.at[i, 'Column6'], sorted_data.at[i, 'Column7'], sorted_data.at[i, 'Column8'], sep='\t', file=annotations_out)
                    print(sorted_data.at[i-1, 'Column0'], sorted_data.at[i-1, 'Column1'], sorted_data.at[i-1, 'Column2'], sorted_data.at[i-1, 'Column3'], sorted_data.at[i-1, 'Column4'], sorted_data.at[i-1, 'Column5'], sorted_data.at[i-1, 'Column6'], sorted_data.at[i-1, 'Column7'], sorted_data.at[i-1, 'Column8'], sep='\t', file=annotations_out)
                    print()
                    
        if match_count > 0:
            print(sorted_data.at[0, 'Column0'], match_count, len(all_ncRNA_data), round(match_count/len(all_ncRNA_data), 2), len(all_mRNA_data), round(match_count/len(all_mRNA_data), 2), sep='\t', file=stats_out)
        else:
            print(sorted_data.at[0, 'Column0'], match_count, len(all_ncRNA_data), 0, len(all_mRNA_data), 0, sep='\t', file=stats_out)
            
stats_out.close()
annotations_out.close()

# calculate genome totals and append to stats output file
scaffold_stats = pd.read_csv('overall_statistics.txt', sep='\t')

total_stats_out = open('overall_statistics.txt', 'a')

print('Genome Totals', scaffold_stats['ncRNA/mRNA Match'].sum(), scaffold_stats['Total ncRNA'].sum(), round(scaffold_stats['ncRNA/mRNA Match'].sum()/scaffold_stats['Total ncRNA'].sum(), 2), scaffold_stats['Total mRNA'].sum(), round(scaffold_stats['ncRNA/mRNA Match'].sum()/scaffold_stats['Total mRNA'].sum(), 2), sep='\t', file=total_stats_out)

total_stats_out.close()

print('Done.')
print()
input("Press Enter to close...")