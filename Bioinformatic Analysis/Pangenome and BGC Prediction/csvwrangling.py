import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import scipy.stats as stats


def parse_fasta_headers(fasta_file):
    '''
    Parse fasta file to get record name and record id for each entry in a multifasta

    Parameters
    ----------
    fasta_file : (string) Fasta file name.

    Returns
    -------
    ids_output : (list) String of id for each record
    descript_output : (list) String of name for each record

    '''
    ids_output = []
    descript_output = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        ids_output.append(record.id.split('.')[0])
        descript_output.append(record.description)
    
    return ids_output, descript_output

## Add consensus protein column for each protein

df = pd.read_csv('04-roary/gene_presence_absence.csv')
if df.astype(str).apply(lambda x: x.str.contains('.p01')).any().any():
    raise Exception("Using a spreadsheet parser or a text parser, please remove all substrings .p01 from the gene_presence_absence.csv file")


fasta_file = '04-roary/pan_genome_reference.fa'
[ids_output,descript_output] = parse_fasta_headers(fasta_file)

df['all_locus_tags'] = df.iloc[:,14:92].astype(str).apply(lambda row: ';'.join(row), axis=1)
df_new = df.drop(df.columns[14:92], axis=1)
consensus_locus_tag = []
for index, row in df_new.iterrows():
    all_locus_tags = row['all_locus_tags']
    all_locus_tags_temp = all_locus_tags.split('\t')
    all_locus_tags_new = [entry.split(';') for entry in all_locus_tags_temp]
    all_locus_tags_list = []
    [all_locus_tags_list.extend(sublist) for sublist in all_locus_tags_new]


    intersection = set(all_locus_tags_list) & set(ids_output)
    consensus_locus_tag.append(', '.join(intersection))
    
df_new['consensus_locus_tag'] = consensus_locus_tag

#%% import database hits

# amrfinder
df1 = pd.read_csv('database_hits/amrfinder_results.csv')
df1 = df1.iloc[:, [0,1]]
df1.columns = ['locus_tag','amrfinder_hit']

# VFDB
df2 = pd.read_csv('database_hits/VFDB_results.csv')
df2 = df2.iloc[:,[3,2]]
df2.columns = ['locus_tag','virulence_hit']
df2 = df2.assign(locus_tag=df2['locus_tag'].str.split('; ')).explode('locus_tag')
df2 = df2[df2['locus_tag'] != '-']

# DefenseFinder
df3 = pd.read_csv('database_hits/defensefinder_results.csv')
df3 = df3.iloc[:,[1,2]]
df3.columns = ['locus_tag','defensefinder_hit']

# DefenseFinder
df4 = pd.read_csv('database_hits/eggNOG_results.csv')
df4 = df4.iloc[:,[0,6]]
df4.columns = ['locus_tag','COG_category']
df4['COG_category'].replace('-',np.nan, inplace=True)


df_data = pd.merge(df1, df2, on='locus_tag', how='outer')
df_data = pd.merge(df_data, df3, on='locus_tag', how='outer')
df_data = pd.merge(df_data, df4, on='locus_tag', how='outer')

df_data = df_data.rename(columns={"locus_tag": "consensus_locus_tag"})

df_final = pd.merge(df_new,df_data, on='consensus_locus_tag', how='outer')


df_filtered = df_final.dropna(subset=['Gene'])
#filter specific lineages in case of contamination
# df_filtered = df_final[~(df_final['all_locus_tags'].str.contains('BJEEGK') & (df_final['No. isolates'] == 1))]

df_filtered.to_csv('pangenome_output_table.csv')




