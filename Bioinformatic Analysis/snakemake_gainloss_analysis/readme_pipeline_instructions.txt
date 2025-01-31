These steps should be performed in order to run the Staphylococcus epidermidis gainloss analysis descriped in this work:

Inputs:
Lineage coassemblies, Sequencing fastqs (see Bioproject #PRJNA1052084) not included here

Steps:
- Run Snakemake pipeline to produce gainloss results


Results:
*masked.fasta: masked coassembly
*.gbff: masked coassembly
*reg*.fasta: sequence of gained or lost region
*regions.mat: table of regions with metadata about contents 
NOTE: relevant results have been precomputed were stored as inputs in the "Pangenome and BGC Prediction" folder

To-dos:
- update for new computing clusters
- build options to reduce what gainloss outputs are included