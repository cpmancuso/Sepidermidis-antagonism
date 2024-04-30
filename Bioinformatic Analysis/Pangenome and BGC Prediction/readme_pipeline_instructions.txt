These steps should be performed in order to run the Staphylococcus epidermidis pangenome and BGC analysis descriped in this work:

Inputs:
Results of genome assembly, gainloss, and annotation pipeline
NOTE: input_genomes is too large for GitHub, see raw data download instructions

Steps:
- Run Snakemake pipeline to produce antismash, BiGScape, and Roary results
- Run extract_consensus_proteins.py
- Submit new_consensus_protein-faa to AMRfinder server
- Submit new_consensus_protein-faa to VFDB server
- Submit new_consensus_protein-faa to DefenseFinder server
- Submit new_consensus_protein-faa to eggnog-mapper server
- Run csvwrangling.py
- Run gainloss.py
- Run stats_enrichment.py
- convert COG category letters to COG category names, using COG dictionary
- collate bigscape families and antismash regions of interest into table

Results:
gainloss_stats_output.csv: enrichment of categories in gainloss regions
combined_stats_output.csv: enrichment of categories in noncore regions

To-dos:
automated local analysis with all .py scripts
automated collation of bigscape and antismash results (pick representative seq, find known_cluster_blast results)
automated annotation of AMR/VFDB/DefenseFinder/eggNOG-mapper (or wait for improvement to Bakta annotation overwrite system)
