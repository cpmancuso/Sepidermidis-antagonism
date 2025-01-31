These steps should be performed in order to run the non-Staphylococcus epidermidis assembly and annotation as descriped in this work:

Inputs:
Raw sequencing fastq files (see Bioproject #PRJNA1052084 and #PRJNA1215987) not included here

Steps:
- Edit samples.csv to point to files
- Run Snakemake pipeline to produce gainloss results


Results:
02-annotation/*/*.gbff: annotated assembly


To-dos:
- merge with other pipelines