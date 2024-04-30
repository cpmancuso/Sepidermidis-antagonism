from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import scipy
import pandas as pd
import numpy as np


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

#%% Get all locus tags from gainloss regions
all_gainloss_locus_tags = []
gbk_files = glob.glob("Input Genomes/*.gbff")
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
pangenome_table = pd.read_csv('pangenome_output_table.csv')
pangenome_table['num_gainlosses'] = 0
for index, row in pangenome_table.iterrows():
    gene_tags = row['all_locus_tags'].split(';')
    intersect = set(gene_tags) & set(all_gainloss_locus_tags)
    pangenome_table['num_gainlosses'][index] = len(intersect)

pangenome_table.to_csv('pangenome_output_table.csv')