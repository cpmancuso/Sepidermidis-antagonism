import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import scipy
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



def parse_proteins_from_genbank(genbank_file, protein_list, descript_list):
    '''
    Check a genbank file for a list of proteins, extract their sequences and return
    
    Parameters
    ----------
    genbank_file : (string) Genbank file name.
    protein_list : (list of strings) protein names to be extracted.
    descript_list : (list of strings) descriptions paired with protein list

    Returns
    -------
    proteins_output : (list) SeqRecords for each protein found in the genbank

    '''
    proteins_output = []
    for contig in SeqIO.parse(genbank_file,'genbank'):
        for record in contig.features:
            if record.type == 'CDS':
                locus_tag = record.qualifiers["locus_tag"][0]
                if locus_tag in protein_list:
                    index = protein_list.index(locus_tag)
                    if "translation" in record.qualifiers.keys():
                        prot_seq = record.qualifiers["translation"][0]
                        descript = descript_list[index]
                        proteins_output.append(SeqRecord(Seq(prot_seq),id=locus_tag,description=descript))
                    else:
                        print(f"{locus_tag} has no translation in genbank file.")
                        prot_seq = 'XX'
                        descript = descript_list[index]
                        proteins_output.append(SeqRecord(Seq(prot_seq),id=locus_tag,description=descript))
                    
    return proteins_output

def recursive_flatten(array):
    result = []
    for item in array:
        if isinstance(item, np.ndarray):
            result.extend(recursive_flatten(item))
        else:
            result.append(item)
    return result

def merge_protein_databases(roary_folder,database_hits_folder):
    '''
    Combines pangenome, consensus proteins, and database hits for amrfinder, 
    VFDB, defensefinder, and eggNOG

    Parameters
    ----------
    roary_folder : folder containing roary outputs, containing 
    'gene_presence_absence.csv' and 'pan_genome_reference.fa'.
    
    database_hits_folder : folder containing 'amrfinder_results.csv', 
    'VFDB_results.csv', 'defensefinder_results.csv', and 'eggNOG_results.csv'. 
    See sample CSVs for formatting of these tables.

    Returns
    -------
    pangenome_table : dataframe containing pangenome table supplemented with 
    database hits and consensus proteins.

    '''    

    ## Add consensus protein column for each protein
    df_temp = pd.read_csv(roary_folder+'/gene_presence_absence.csv')
    if df_temp.astype(str).apply(lambda x: x.str.contains('.p01')).any().any():
        raise Exception("Using a spreadsheet parser or a text parser, please remove all substrings .p01 from the gene_presence_absence.csv file")
    
    
    fasta_file = roary_folder+'/pan_genome_reference.fa'
    [ids_output,descript_output] = parse_fasta_headers(fasta_file)
    
    df_temp['all_locus_tags'] = df_temp.iloc[:,14:92].astype(str).apply(lambda row: ';'.join(row), axis=1)
    df_new = df_temp.drop(df_temp.columns[14:92], axis=1)
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
    
    # import database hits
    
    # amrfinder
    df1 = pd.read_csv(database_hits_folder+'/amrfinder_results.csv')
    df1 = df1.iloc[:, [0,1]]
    df1.columns = ['locus_tag','amrfinder_hit']
    
    # VFDB
    df2 = pd.read_csv(database_hits_folder+'/VFDB_results.csv')
    df2 = df2.iloc[:,[3,2]]
    df2.columns = ['locus_tag','virulence_hit']
    df2 = df2.assign(locus_tag=df2['locus_tag'].str.split('; ')).explode('locus_tag')
    df2 = df2[df2['locus_tag'] != '-']
    
    # DefenseFinder
    df3 = pd.read_csv(database_hits_folder+'/defensefinder_results.csv')
    df3 = df3.iloc[:,[1,2]]
    df3.columns = ['locus_tag','defensefinder_hit']
    
    # DefenseFinder
    df4 = pd.read_csv(database_hits_folder+'/eggNOG_results.csv')
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
    pangenome_table = df_filtered
    return pangenome_table

def append_gainloss_results(pangenome_table,genome_folder,gainloss_folder):
    '''
    Extracts relevant data from gainloss folders and appends to pangenome table.

    Parameters
    ----------
    pangenome_table : dataframe containing pangenome table supplemented with 
    database hits and consensus proteins.
    
    genome_folder : folder containing annotated genomes as '{lineage}.gbff files.
    
    gainloss_folder : folder containing gainloss outputs as '{lineage}_regions.mat files.

    Returns
    -------
    pangenome_table : dataframe containing pangenome table supplemented with 
    database hits and consensus proteins.

    '''
    #%% Get all locus tags from gainloss regions
    all_gainloss_locus_tags = []
    gbk_files = glob.glob(genome_folder+"/*.gbff")
    # gbk_files = ['input_genomes\sepi_clade_1.gbff']
    for genbank_file in gbk_files:
        filename = genbank_file.split('\\')[-1]
        lineage = filename.split('.')[0]
        
        for contig in SeqIO.parse(genbank_file,'genbank'):
            for record in contig.features:
                if record.type == 'CDS':
                    locus_tag = record.qualifiers["locus_tag"][0]
                    locus_tag_header = locus_tag.split('_')[0]
                    break
        
        print(lineage)
        mat_data = scipy.io.loadmat(f'gainloss_results/{lineage}_regions.mat')
        table = mat_data['all_regions_table']
        
        if table.size:
            for n in range(table.shape[1]):
                region = table[:,n]
                if region[0]:
                    locusnums = region['Gene_locustags']
                    locusnums_list = recursive_flatten(locusnums)
                    locus_tag_list = [locus_tag_header + '_' + s for s in locusnums_list]
                    all_gainloss_locus_tags.extend(locus_tag_list[:])
    
    #%% Add data to csv table
    pangenome_table['num_gainlosses'] = 0
    for index, row in pangenome_table.iterrows():
        gene_tags = row['all_locus_tags'].split(';')
        intersect = set(gene_tags) & set(all_gainloss_locus_tags)
        pangenome_table['num_gainlosses'][index] = len(intersect)
    
    pangenome_table.to_csv('pangenome_output_table.csv')
    return pangenome_table


def main():
    roary_folder = '04-roary'
    database_hits_folder = 'database_hits'
    genome_folder = 'Input Genomes'
    gainloss_folder = 'gainloss_results'
    
    ## Create pangenome table
    pangenome_table = merge_protein_databases(roary_folder,database_hits_folder)
    updated_pangenome_table = append_gainloss_results(pangenome_table,genome_folder,gainloss_folder)
    
    
    ## Calculate stats using pangenome table
    df = updated_pangenome_table
    num_lineages = max(df.loc[:, 'No. isolates'])
    core_thresh = 0.95*num_lineages
    shell_thresh = 0.15*num_lineages


    #%% calculate enrichment of gene group
    total_core = sum(df.loc[:, 'No. isolates']>core_thresh)
    total_shell = sum(df.loc[:, 'No. isolates']>shell_thresh)-total_core
    total_cloud = sum(df.loc[:, 'No. isolates']>0)-total_core-total_shell
    group = ['Total']


    columns = ['amrfinder_hit','virulence_hit','defensefinder_hit']
    COG_dictionary = pd.read_csv('COG_dictionary.csv')

    #%% Calculate enrichment in recent gains/losses
    total_gainloss = sum(df.loc[:, 'num_gainlosses']>0)
    group = ['Total']
    core_num = [total_core]
    gainloss_num = [total_gainloss]
    odds = [0]
    p = [0]

    for col in columns:
        non_nan_indices = np.where(df[col].notna())
        core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
        pseudo=False
        if core==0:
            pseudo = True
            core=1 #add pseudo count
        gainloss = sum((df.loc[non_nan_indices, 'num_gainlosses']>0))
        odds_ratio, p_value = stats.fisher_exact([[gainloss, core],[total_gainloss-gainloss,total_core-core]])
        group.append(col)
        if pseudo:
            core=0        
        core_num.append(core)
        gainloss_num.append(gainloss)
        odds.append(odds_ratio)
        p.append(p_value)

    for code in COG_dictionary['Code']:
        non_nan_indices = np.where(df['COG_category'].fillna('').str.contains(code))
        core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
        gainloss = sum((df.loc[non_nan_indices, 'num_gainlosses']>0))
        odds_ratio, p_value = stats.fisher_exact([[gainloss, core],[total_gainloss-gainloss,total_core-core]])
        group.append(code)
        core_num.append(core)
        gainloss_num.append(gainloss)
        odds.append(odds_ratio)
        p.append(p_value)
        
        
    df_stats = pd.DataFrame({'Group':group,'Core Number':core_num,'Gainloss Number':gainloss_num,'Odds Ratio':odds,'p_value':p})
    df_stats.to_csv('gainloss_stats_output.csv')


    #%% Calculate enrichment in non-core
    total_noncore = total_shell+total_cloud
    group = ['Total']
    core_num = [total_core]
    noncore_num = [total_noncore]
    odds = [0]
    p = [0]

    for col in columns:
        non_nan_indices = np.where(df[col].notna())
        core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
        pseudo=False
        if core==0:
            pseudo = True
            core=1 #add pseudo count
        noncore = sum(((df.loc[non_nan_indices, 'No. isolates']>0)&(df.loc[non_nan_indices, 'No. isolates']<=core_thresh))|(df.loc[non_nan_indices, 'num_gainlosses']>0))
        odds_ratio, p_value = stats.fisher_exact([[noncore, core],[total_noncore-noncore,total_core-core]])
        group.append(col)
        if pseudo:
            core=0  
        core_num.append(core)
        noncore_num.append(noncore)
        odds.append(odds_ratio)
        p.append(p_value)

    for code in COG_dictionary['Code']:
        non_nan_indices = np.where(df['COG_category'].fillna('').str.contains(code))
        core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
        noncore = sum(((df.loc[non_nan_indices, 'No. isolates']>0)&(df.loc[non_nan_indices, 'No. isolates']<=core_thresh))|(df.loc[non_nan_indices, 'num_gainlosses']>0))
        odds_ratio, p_value = stats.fisher_exact([[noncore, core],[total_noncore-noncore,total_core-core]])
        group.append(code)
        core_num.append(core)
        noncore_num.append(noncore)
        odds.append(odds_ratio)
        p.append(p_value)
        
        
    df_stats = pd.DataFrame({'Group':group,'Core Number':core_num,'Noncore Number':noncore_num,'Odds Ratio':odds,'p_value':p})
    df_stats.to_csv('accessory_stats_output.csv')

if __name__ == "__main__":
    main()






