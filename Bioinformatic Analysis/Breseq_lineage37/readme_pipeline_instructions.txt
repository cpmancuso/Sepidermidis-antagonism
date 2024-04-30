These steps should be performed in order to run the breseq analysis descriped in this work:

Inputs:
fastq files (see Bioproject #SUB14054091, not included here), genome assembly (*.gbff files)

Steps:
- Run Snakemake pipeline to produce breseq results

Results:
Output/index.html, in depth analysis of isolate 37.3, aka stock 313
Output/breseq_comparison.html, comparison of all lineage 37 isolates and colonies selected on polymyxin B


To-dos:
Change snakemake to keep only desired output files, since breseq produces many