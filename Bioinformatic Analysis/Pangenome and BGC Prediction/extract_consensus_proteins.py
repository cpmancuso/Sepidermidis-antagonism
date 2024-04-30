from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob


def extract_proteins_from_genbank(genbank_file, protein_list, descript_list):
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



# Make list of consensus proteins to look for
fasta_file = '04-roary/pan_genome_reference.fa'
[ids_output,descript_output] = parse_fasta_headers(fasta_file)
consensus_proteins = []

# append proteins to list
gbk_files = glob.glob("Input Genomes/*.gbff")
for genbank_file in gbk_files:
    print(genbank_file)
    proteins_output = extract_proteins_from_genbank(genbank_file,ids_output,descript_output)
    consensus_proteins.extend(proteins_output)

# write to file
SeqIO.write(consensus_proteins,'consensus_proteins.faa','fasta')