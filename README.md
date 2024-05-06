**Intraspecies warfare restricts strain coexistence in human skin microbiomes** \
The manuscript corresponding to this repository has not been peer reviewed, and this repository is likely to undergo additional ammendments before publication.

This repository contains all code necessary to reproduce the analyses described in:

> **Intraspecies warfare restricts strain coexistence in human skin microbiomes** \
> Christopher P. Mancuso, Jacob S. Baker, Evan Qu, A. Delphine Tripp, Ishaq O. Balogun, Tami D. Lieberman

This project is part of my postdoctoral research in the Lieberman Lab.

**Overview** \
This repository is grouped into three major parts:
- **Bioinformatic Analysis**: contains the snakemake pipelines and processed data outputs from the bioinformatic analyses. As these scripts are written for execution on a particular computing cluster, they will need significant modifications before running on your system.
- **Image Analysis**: contains the MATLAB pipeline and processed image outputs used in the aforementioned manuscript. This code should work on your system, provided you have the requisite MATLAB toolboxes.
- **Interaction Analysis**: contains the MATLAB pipeline and data sturctures needed to generate all the plots in the manuscript. This code should work on your system, provided you have the requisite MATLAB toolboxes.

**Data Availability** \
Intermediate and processed data files for each set of pipelines are included in each folder. \
Original image files for the antagonism screens are saved at: https://doi.org/10.6084/m9.figshare.25726482.v1

Sequencing and genome assemblies were performed as part of "Highly-resolved within-species dynamics in the human facial skin microbiome", see Jacob's repo: https://github.com/jsbaker760/highres_dynamics 

Raw sequencing data for running Snakemakes is currently being made available on the NCBI Sequence Read Archive under Bioproject #SUB14054091.
