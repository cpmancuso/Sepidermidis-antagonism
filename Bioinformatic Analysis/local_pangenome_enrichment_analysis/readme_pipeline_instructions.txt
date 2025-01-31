These steps should be performed in order to run the Staphylococcus epidermidis pangenome and BGC enrichment statistical analysis descriped in this work:

Inputs:
Results of genome assembly, gainloss, and annotation pipeline
NOTE: input_genomes folder is too large for GitHub, see raw data download instructions

Steps:
- Download Antismash/BiGScape/Roary and gain/loss results from respective snakemake folders
- Run extract_consensus_proteins.py
- Submit new_consensus_protein-faa to AMRfinder server
- Submit new_consensus_protein-faa to VFDB server
- Submit new_consensus_protein-faa to DefenseFinder server
- Submit new_consensus_protein-faa to eggnog-mapper server
- Run calculate_pangenome_enrichment.py
- convert COG category letters to COG category names, using COG dictionary

Results:
gainloss_stats_output.csv: enrichment of categories in gainloss regions
accessory_stats_output.csv: enrichment of categories in noncore regions

To-dos:
automated annotation of AMR/VFDB/DefenseFinder/eggNOG-mapper (or wait for improvement to Bakta annotation overwrite system)
