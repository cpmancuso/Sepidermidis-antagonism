**Intraspecies warfare restricts strain coexistence in human skin microbiomes** \
The manuscript corresponding to this repository has not been peer reviewed, and this repository is likely to undergo additional ammendments before publication.

This repository contains all code necessary to reproduce the analyses described in:

> **Intraspecies warfare restricts strain coexistence in human skin microbiomes** \
> Christopher P. Mancuso, Jacob S. Baker, Evan Qu, A. Delphine Tripp, Ishaq O. Balogun, Tami D. Lieberman

This project is part of my postdoctoral research in the Lieberman Lab.

**Overview** \
This repository is grouped into three major parts:
- **Bioinformatic Analysis**: contains the snakemake pipelines (to be run on a SLURM computing cluster) and processed data outputs from the bioinformatic analyses. 
- **Image Analysis**: contains the MATLAB pipeline and processed image outputs used in the aforementioned manuscript. 
- **Interaction Analysis**: contains the MATLAB pipeline and data sturctures needed to generate all the plots in the manuscript. 

**Instructions** \
See "instructions_to_run" document in this folder

**Dependencies** \
**Bioinformatic Analysis**: All conda environments used for SLURM cluster operation are included in the relevant snakemake folders. Installation of all required conda environments may take up to 1-2h depending on your system and the download speed of external databases. As these scripts are written for execution on the MIT Engaging computing cluster (a SLURM system), they will need significant modifications before running on your system. The v.1.0 (original GitHub commit) pipeline was tested on Snakemake version 7.32.3. \
**Image & Interaction Analysis**: All MATLAB scripts are included, and were tested on a MATLAB 2023a environment with the Bioinformatics, Image Processing, and Statistics & Machine Learning toolboxes.

**Data Availability** \
Intermediate and processed data files for each set of pipelines are included in each folder. \
Original image files for the antagonism screens are saved at: https://doi.org/10.6084/m9.figshare.25726482.v1

Sequencing and genome assemblies were performed as part of "Highly-resolved within-species dynamics in the human facial skin microbiome", see Jacob's repo: https://github.com/jsbaker760/highres_dynamics 

Raw sequencing data for running Snakemakes is currently being made available on the NCBI Sequence Read Archive under Bioproject #PRJNA1052084 and #PRJNA1215987.

**Expected Run Times** \
See readme in each folder for additiona instructions.
- **Bioinformatic Analysis**: the pipeline takes up to 12h to complete on the MIT Engaging Cluster. Final results (but not intermediate data files) are included in this repository.
- **Image Analysis**: the code takes up to 1h per set of 96 images, including hands-on guided correction of mis-identified colonies. Pre-analyzed data files are included in this repository.
- **Interaction Analysis**: the code takes <5 minutes to generate all the figures in the manuscript.
