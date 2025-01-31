%% MOBILE ELEMENT GAIN/LOSS ANALYSIS
function gainloss_snakemake(this_clade_name,dir_data_assemblies,dir_output)

%% Summary
% This script uses a coverage matrix over a lineage pan-genome assembly in
% order to to identify candidate gain/loss regions.

%% Directory setup
% Main directory:
dir_main = char(pwd);
path( dir_main,path );

% Directory for my scripts:
dir_myscripts = [dir_main '/miniscripts' ];
path(dir_myscripts,path);

%% Where to find data on lineage pan-genome assemblies
% {this_clade_name} e.g. sepi_clade_1

% dir_data_assemblies contains
% 1. {this_clade_name}.fasta: fasta nucleotide sequence, tested on contigs.fasta from spades
% 2. {this_clade_name}.gbff: genbank nulceotide annotation file, tested on output from Bakta v1.6.1
% 3. {this_clade_name}_coverage_matrix.mat: coverage matrix output from get_all_coverage_snakemake() containing SampleNames, all_coverage_per_bp, all_maf_per_bp, cov_modes

%% Output directory
dir_save_coverage = [dir_output '/' this_clade_name '/coverage'];
if ~exist( dir_save_coverage, 'dir' )
    mkdir( dir_save_coverage )
end

dir_save_regions = [dir_output '/' this_clade_name '/regions'];
if ~exist( dir_save_regions, 'dir' )
    mkdir( dir_save_regions )
end

%% Identify Gain/Loss Elements
[all_regions_table, num_regions_final] = gainloss_lineage( this_clade_name,dir_data_assemblies, dir_save_coverage, dir_save_regions);
% exit(0)