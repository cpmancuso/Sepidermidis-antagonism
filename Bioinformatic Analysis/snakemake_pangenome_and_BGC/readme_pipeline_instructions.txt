These steps should be performed in order to run the pangenome and BGC analysis descriped in this work:

Inputs:
fastq files (see Bioproject #PRJNA1052084, not included here), genome assembly (*.gbff files)

Steps:
- Edit samples.csv to point to files
- Run Snakemake pipeline to produce antiSMASH, BiGScape, and Roary results and BGC summary tables

Results:
02-regions: antismash region contigs folder
03-bigscape: bigscape BGC clustering results
04-roary: pangenome data folder
05-summarize-results/BGCs_collated_by_family.csv: summary table of all BGCs clustered by BigScape
05-summarize-results/BGCs_collated_by_source.csv: summary table of all BGCs clustered by knownclusterblast results, even if not clustered by BiGScape

NOTE: pre-computed results folders 02/03/04 are too large for GitHub, see data download instructions

To-dos:
Change snakemake to keep only desired output files, since analyses produce many