######################
# GAINLOSS SNAKEMAKE #
######################

# Version History:
# # 2023.1 Chris: Merged Jacob Case snakemake with Arolyn Gain Loss scripts


''' PRE-SNAKEMAKE '''
import sys, subprocess

MY_SCRIPTS_DIRECTORY = "./scripts"
sys.path.insert(0, MY_SCRIPTS_DIRECTORY)
SCRIPTS_DIRECTORY = "/scratch/mit_lieberman/scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)

spls = "samples_gainloss.csv"

import basic_snakemake_functions as bsf # generic functions to read samples.csv etc.
from itertools import compress
import sys

# Define couple of lists from samples_case.csv
# # Path = points to snakemake mapping diversity folder
# # Sample = sample name; must match mapping step
# # ReferenceGenome = Not used
# # Outgroup = Not used
# # Clade = the number of the lineage to which the sample belongs, name of folder for ReferenceGenome
[PATH_ls, SAMPLE_ls, REF_Genome_ls, OUTGROUP_ls, CLADEID_ls] = bsf.read_samplesCSV(spls)
# # The unique clades
CLADES_ls = set(CLADEID_ls)
REFGENOMES_ls = set(REF_Genome_ls)

''' SNAKEMAKE '''
rule all:
    input:
        # # Data links only # #
        expand("data/diversity/{sampleID}_ref_{cladeID}.diversity.mat",zip,sampleID=SAMPLE_ls, cladeID=CLADES_ls),
        expand("2-clade_coverage_matrix/{cladeID}/contigs.fasta",cladeID=CLADES_ls),
        expand("2-clade_coverage_matrix/{cladeID}/{cladeID}.gbff",cladeID=CLADES_ls),
        # # Coverage Matrix # #
        expand("2-clade_coverage_matrix/{cladeID}/coverage_matrix.mat",cladeID=CLADES_ls),
        # # Gain Loss # #
        expand("3-gain_loss_output/{cladeID}/regions/{cladeID}_regions.mat",cladeID=CLADES_ls),

''' HELPER FUNCTIONS for switching from sample-wise rules to clade-wise rules'''

def get_clade_wildcards(cladeID):
    is_clade = [int(i == cladeID) for i in CLADEID_ls]
    sampleID_clade = list(compress(SAMPLE_ls,is_clade))
    reference_clade = list(compress(REF_Genome_ls,is_clade))
    outgroup_clade = list(compress(OUTGROUP_ls,is_clade))
    return sampleID_clade,reference_clade,outgroup_clade
    
def get_sampleID_names(wildcards):  
    sampleID_clade,_,_ = get_clade_wildcards(wildcards.cladeID)
    return sampleID_clade

def get_ref_genome(wildcards):
    _,reference_clade,_ = get_clade_wildcards(wildcards.cladeID)
    return reference_clade

def get_diversity(wildcards):
    sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
    diversity_mat = expand("data/diversity/{sampleID}_ref_{cladeID}.diversity.mat",sampleID=sampleID_clade, cladeID=wildcards.cladeID)
    return diversity_mat   

# make data links

#single process for all diversity, since need to make on per-sample basis
rule make_all_sample_data_links:
    output:
        div_mat_links = expand("data/diversity/{sampleID}_ref_{cladeID}.diversity.mat",zip,sampleID=SAMPLE_ls, cladeID=CLADEID_ls),
    run:
        # subprocess.run( "rm -fr data/ " ,shell=True) # clean it up prior run
        # subprocess.run( "mkdir -p data/ " ,shell=True)
        for i in range(len(SAMPLE_ls)):
            sourcefile = PATH_ls[i] + SAMPLE_ls[i] + "_ref_" + CLADEID_ls[i] + "_aligned.sorted.strain.variant.diversity.mat"
            linkfile = " data/diversity/" + SAMPLE_ls[i] + "_ref_" + CLADEID_ls[i] + ".diversity.mat"
            subprocess.run( "ln -fs -T " + sourcefile + linkfile,shell=True)

rule make_assembly_links:
    output:
        clade_fasta = "2-clade_coverage_matrix/{cladeID}/{cladeID}.fasta",
        clade_genbank = "2-clade_coverage_matrix/{cladeID}/{cladeID}.gbff",
    params:
        ref_genome_folder = "/scratch/mit_lieberman/projects/cm_staph/2302_acera_gainloss/clade_assemblies/",
    run:
        subprocess.run("ln -fs -T " + params.ref_genome_folder + wildcards.cladeID + "/contigs.fasta "+output.clade_fasta ,shell=True)
        subprocess.run("ln -fs -T " + params.ref_genome_folder + wildcards.cladeID +"/"+ wildcards.cladeID + ".gbff " +output.clade_genbank ,shell=True)


# Build input for coverage matrix
rule string_sampleID_names:
    params:
        sampleID_names = get_sampleID_names,
    output:
        string_sampleID_names = "1-temp_pos/clade_{cladeID}_string_sampleID_names.txt",
    run:
        with open( output.string_sampleID_names, "w") as f: 
            print(*params.sampleID_names, sep=" ", file=f)

rule string_diversity_mat:
    input:
        diversity_mat = get_diversity,
    output:
        string_diversity_mat = "1-temp_pos/clade_{cladeID}_string_diversity_mat.txt",
    run:
        with open( output.string_diversity_mat ,"w") as f: 
            print(*input.diversity_mat, sep=" ", file=f)

rule coverage_matrix_mat:
    input: 
        string_sampleID_names = "1-temp_pos/clade_{cladeID}_string_sampleID_names.txt",
        string_diversity_mat = "1-temp_pos/clade_{cladeID}_string_diversity_mat.txt",
    output:
        coverage_mat = "2-clade_coverage_matrix/{cladeID}/coverage_matrix.mat",
    log:
        'logs/{cladeID}_coverage_matrix.log',
    shell:
        """
        module add mit/matlab/2015b; 
        matlab -r "path('{SCRIPTS_DIRECTORY}',path); get_all_coverage_snakemake( '{input.string_sampleID_names}', '{input.string_diversity_mat}', '{output.coverage_mat}' )" > {log} 2>&1
        """

# Run gain loss on each clade

rule gain_loss:
    input:
        coverage_mat = rules.coverage_matrix_mat.output.coverage_mat,
        clade_genbank = rules.make_assembly_links.output.clade_genbank,
        clade_fasta = rules.make_assembly_links.output.clade_fasta,
    output:
        regions = "3-gain_loss_output/{cladeID}/regions/{cladeID}_regions.mat",
    params:
        data_dir = "2-clade_coverage_matrix/{cladeID}/",
        output_dir = "3-gain_loss_output/", #matlab appends cladeID
    log:
        'logs/{cladeID}_gainloss.log',
    shell:
        """
        module add mit/matlab/2018a; 
        matlab -r "gainloss_snakemake( '{wildcards.cladeID}' , '{params.data_dir}' , '{params.output_dir}' )" > {log} 2>&1
        """