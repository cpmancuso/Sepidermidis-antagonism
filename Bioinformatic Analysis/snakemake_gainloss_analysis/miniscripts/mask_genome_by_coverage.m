function mask_genome_by_coverage(copynum_cutoff, spades_cutoff, mask_custom, this_clade_name, dir_data_assemblies, path_genome,path_genbank,path_covmat,dir_save_coverage)

save_masked_files = 1;
length_cutoff = 200;

% Get genome stats
[ChrStarts, GenomeLength, ~, genome_contig_names, genome_contig_seqs] = get_genome_stats(path_genome);
contig_num_by_index = p2chrpos(1:1:GenomeLength,ChrStarts); contig_num_by_index = contig_num_by_index(:,1);
contig_start_indices = ChrStarts+1; % since ChrStarts starts with 0
contig_end_indices = [ ChrStarts(2:end), GenomeLength ];
num_contigs = numel(ChrStarts);

contig_lengths = cellfun(@(x) numel(x), genome_contig_seqs );
for j=1:num_contigs
    name_split = strsplit(genome_contig_names{j},'_');
    %spades saves the kmer coverage of each contig in the contig name
    spades_contig_coverage(j) = str2num(name_split{end});
end

% Annotations
if ~exist([dir_data_assemblies '/cds_sorted.mat'], 'file')
    genbank_to_CDS(dir_data_assemblies,path_genome,path_genbank)
end % makes cds_sorted.mat if it doesn't already exist
load([dir_data_assemblies '/cds_sorted.mat']) 

% Load coverage matrix
load( path_covmat ) % all_coverage_per_bp, SampleNames
clear all_maf_per_bp % not used in this version
clear cov_modes % not used in this version
num_samples = numel(SampleNames);

%% Compute coverage statistics and normalized coverage matrix
fprintf(1,['Generating coverage matrices...' '\n'])

% Calculate mean / stdev based on top 2Mb of contigs
% stats_contigs = find(contig_start_indices<2000000);

% Calculate mean / stdev based on top 10% of long contigs by coverage
long_contigs = contig_lengths>20000;
stats_contigs = find(spades_contig_coverage(long_contigs)>prctile(spades_contig_coverage(long_contigs),90)); %assumes contigs are ordered longest to shortest
stats_positions = [];
fprintf(1,['Calculating based on: ' num2str(numel(stats_contigs)) ' contigs' '\n'])
for n=1:numel(stats_contigs)
    stats_positions = [stats_positions contig_start_indices(stats_contigs(n)):contig_end_indices(stats_contigs(n))];
end

avg_cov_by_sample = mean(all_coverage_per_bp(:,stats_positions),2); 
stdev_cov_by_sample = std(single(all_coverage_per_bp(:,stats_positions)),0,2); 
median_cov_by_sample = median(all_coverage_per_bp(:,stats_positions),2); 
difference_mean_in_sdunits = (single(all_coverage_per_bp) - repmat(avg_cov_by_sample,1,size(all_coverage_per_bp,2))) ./ repmat(stdev_cov_by_sample,1,size(all_coverage_per_bp,2)) ;

%% Produce coverage statistics for each contig 
fprintf(1,['Finding low coverage contigs ...' '\n'])
all_coverage_per_bp_copynum = single(all_coverage_per_bp)./single(median_cov_by_sample);
median_copy_by_contig=single(zeros(num_samples,num_contigs));
norm_median_cov_by_contig=single(zeros(num_samples,num_contigs));

for i=1:num_samples
    for j=1:numel(contig_start_indices)   
        median_copy_by_contig(i,j) = median(all_coverage_per_bp_copynum(i,contig_start_indices(j):contig_end_indices(j)));
        norm_median_cov_by_contig(i,j) = median(difference_mean_in_sdunits(i,contig_start_indices(j):contig_end_indices(j)));
    end
end

%mask samples with suspiciously high number of contigs
contigs_present_per_sample = sum(median_copy_by_contig > copynum_cutoff,2);
contaminated_samples = contigs_present_per_sample > mean(contigs_present_per_sample)+2*std(contigs_present_per_sample);

median_copy_by_contig(contaminated_samples,:)=NaN;
norm_median_cov_by_contig(contaminated_samples,:)=NaN;

for i=1:num_samples
    if contaminated_samples(i)
        fprintf(1,['Sample ' SampleNames{i} ' has ' num2str(contigs_present_per_sample(i)) ' high cov contigs and may be contaminated...' '\n'])
    else
        fprintf(1,['Sample ' SampleNames{i} ' has ' num2str(contigs_present_per_sample(i)) ' high cov contigs...' '\n'])
    end
end

%% Mask genome
% calculate a lineage specific cutoff, but note that spades filtering is
% stored, not applied
spades_cutoff = spades_cutoff.*median(spades_contig_coverage(stats_contigs));

% keep contigs with high spades coverage OR high coverage in at least one sample
mask_copynum = (max(median_copy_by_contig)>=copynum_cutoff);
mask_spades = spades_contig_coverage>=spades_cutoff;
% also mask short contigs
mask_length = contig_lengths>=length_cutoff;
% mask = (mask_copynum | mask_spades); %if want to trust spades, risks
% contamination
mask = mask_copynum;
mask = mask & mask_length;
mask = mask & mask_custom;
mask_list=find(mask);
pos_mask=ismember(contig_num_by_index,mask_list)';
fprintf(1,['Plotting masked contigs...' '\n'])
make_clade_coverage_plot(SampleNames, this_clade_name, contig_start_indices, contig_lengths,...
    median_copy_by_contig, norm_median_cov_by_contig,...
    mask_copynum,mask_spades,mask_length,mask_custom,mask,...
    copynum_cutoff,spades_cutoff,length_cutoff, spades_contig_coverage, dir_save_coverage)

%% Produce masked genome and annotations

%NOTE: this creates masked versions of all files which can be loaded later for reanalysis
if save_masked_files
    fprintf(1,['Saving masked data...' '\n'])
    % masks
    save([dir_save_coverage '/' this_clade_name '_masks.mat'],'mask_copynum','mask_spades','mask_length','mask_custom','mask','copynum_cutoff','spades_cutoff','length_cutoff')
    
    % genome
    masked_fasta = struct('Sequence',genome_contig_seqs(mask),'Header',genome_contig_names(mask));
    output_fasta = [ dir_save_coverage '/' this_clade_name '_masked.fasta' ];
    if exist( output_fasta, 'file' ) % remove if already exists
        delete( output_fasta )
    end
    fastawrite(output_fasta,masked_fasta)
    close all

    % annotations
    CDS=CDS(mask);
    save([dir_save_coverage '/' this_clade_name '_cds_sorted_masked.mat'],'CDS')

    % coverage
    all_coverage_per_bp=all_coverage_per_bp(:,pos_mask);
    save([dir_save_coverage '/' this_clade_name '_coverage_matrix_masked.mat'],'all_coverage_per_bp','SampleNames','-v7.3');
end