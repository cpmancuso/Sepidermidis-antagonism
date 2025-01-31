function [ all_regions_table, num_regions_final ] = ...
    gainloss_lineage( this_clade_name, ...
    dir_data_assemblies, dir_save_coverage, dir_save_regions )

%% Summary

% This function detects gain/loss regions in a lineage

%% Where to find data on lineage pan-genome assemblies
% {this_clade_name} e.g. sepi_clade_1

% dir_data_assemblies contains
% 1. {this_clade_name}.fasta: fasta nucleotide sequence, tested on contigs.fasta from spades
% 2. {this_clade_name}.gbff: genbank nulceotide annotation file, tested on output from Bakta v1.6.1
% 3. {this_clade_name}_coverage_matrix.mat: coverage matrix output from get_all_coverage_snakemake() containing SampleNames, all_coverage_per_bp, all_maf_per_bp, cov_modes

% Paths to downloaded files
path_genome = [ dir_data_assemblies '/' this_clade_name '.fasta'];
path_genbank = [ dir_data_assemblies '/' this_clade_name '.gbff'];
path_covmat =  [ dir_data_assemblies '/coverage_matrix.mat' ];

%% Load data for this cluster

% Mask low coverage contigs
if ~exist([dir_save_coverage '/' this_clade_name '_coverage_matrix_masked.mat'],'file')
    fprintf(1,['Masking genome for ' this_clade_name '...' '\n'])
    load(path_covmat,'SampleNames')
    % masking parameters
    copynum_cutoff = 0.5;
    spades_cutoff = 0.05; %fraction of max contigs
    % input custom mask, placeholder
    
    mask_custom = true(1,numel(fastaread(path_genome)));
    path_annotations = [dir_save_coverage '/' this_clade_name '_cds_sorted_masked.mat'];

    mask_genome_by_coverage(copynum_cutoff, spades_cutoff, mask_custom, this_clade_name, dir_data_assemblies, path_genome,path_genbank,path_covmat,dir_save_coverage)
end

% point to masked files instead
path_genome = [dir_save_coverage '/' this_clade_name '_masked.fasta'];
path_annotations = [dir_save_coverage '/' this_clade_name '_cds_sorted_masked.mat'];
path_covmat = [dir_save_coverage '/' this_clade_name '_coverage_matrix_masked.mat'];

% Get genome stats from fasta file
[ChrStarts, GenomeLength, ~, ~, genome_contig_seqs] = get_genome_stats(path_genome);
contig_num_by_index = p2chrpos(1:1:GenomeLength,ChrStarts); contig_num_by_index = contig_num_by_index(:,1);
contig_start_indices = ChrStarts+1; % since ChrStarts starts with 0
contig_end_indices = [ ChrStarts(2:end), GenomeLength ];
contig_lengths = cellfun(@(x) numel(x), genome_contig_seqs );

% Load annotations and coverage matrix 
fprintf(1,['Loading masked data for ' this_clade_name '...' '\n'])
load(path_annotations, 'CDS') 
load(path_covmat,'all_coverage_per_bp','SampleNames')
num_samples = numel(SampleNames);

%% Compute coverage statistics and normalized coverage matrix
fprintf(1,['Generating coverage matrices...' '\n'])

% Calculate mean / stdev based on masked contigs
avg_cov_by_sample = mean(single(all_coverage_per_bp),2);
stdev_cov_by_sample = std(single(all_coverage_per_bp),0,2);
median_cov_by_sample = median(single(all_coverage_per_bp),2);
all_coverage_per_bp_copynum = single(all_coverage_per_bp)./single(median_cov_by_sample);

% Double-normalized coverage matrix
start_pos=1; %previously used for buffering
% First normalize by sample
difference_mean_in_sdunits = (single(all_coverage_per_bp) - repmat(avg_cov_by_sample,1,size(all_coverage_per_bp,2))) ./ repmat(stdev_cov_by_sample,1,size(all_coverage_per_bp,2)) ;
% Then normalize by position normalize by sample
avg_cov_bypos = mean(difference_mean_in_sdunits);
cov_mat_doublenorm = (single(difference_mean_in_sdunits) - repmat(avg_cov_bypos,size(all_coverage_per_bp,1),1)) ;

%% Find candidate gain/loss events
fprintf(1,['Finding candidate gain/loss events...' '\n'])

% Parameters: regions absent
min_size_loss = 2000; % bp
loose_threshold_for_loss = 0.25; % copy number
max_cutoff_for_loss = 0.15; % copy number
max_avg_cov_for_loss = 2.5; % reads
min_avg_copynum_for_loss_high_control = 0.75; % copy number

% Parameters: regions present
% min_size_gain = min_size_loss; % bp
loose_threshold_for_gain = 0.5; % copy number % 0.5
min_cutoff_for_gain = min_avg_copynum_for_loss_high_control;  % copy number
% min_avg_cov_for_gain = 10; % reads
% max_avg_copynum_for_gain_low_control = max_cutoff_for_loss; % copy number
% max_avg_cov_for_gain_low_control = max_avg_cov_for_loss; % copy number

% Parameters: testing contigs
max_contig_len_to_test = GenomeLength; % test all contigs (even though obviously long ones won't show up as gain/losses)

% Sample filter
min_cov_to_eval_sample = 20; % reads

% Find positions/samples with high/low coverage
has_high_coverage = all_coverage_per_bp_copynum > loose_threshold_for_gain; 
has_low_coverage = all_coverage_per_bp_copynum < loose_threshold_for_loss;
has_high_coverage(:,1:2)=0;  has_high_coverage(:,end-1:end)=0; %force the first and last positions to be 0 so that there are the same number of starts and ends
has_low_coverage(:,1:2)=0;  has_low_coverage(:,end-1:end)=0; %force the first and last positions to be 0 so that there are the same number of starts and ends

% Initialize
copy_number_variants_strains=[];
copy_number_variants_ends=[];
copy_number_variants_starts=[];
isdel=[];
outputfile = fopen('output.txt','w');
% Find candidate gains and losses
for i=1:num_samples
    fprintf(outputfile,['Sample ' num2str(i) '/' num2str(num_samples) ': ' SampleNames{i} '\n']);
    fprintf(1,['Sample ' num2str(i) '/' num2str(num_samples) ': ' SampleNames{i} '\n'])
    
    if avg_cov_by_sample(i)>min_cov_to_eval_sample % only looking at strains with high enough coverage

        % Candidate gain/loss regions in this sample: 
        del_starts_0=find(diff(has_low_coverage(i,:))>0)+1; % diff is one shorter 
        del_ends_0=find(diff(has_low_coverage(i,:))<0);
        % Impose contig boundaries
        [ del_starts, del_ends ] = break_at_contig_boundaries( del_starts_0, del_ends_0, ...
            contig_num_by_index, contig_start_indices, contig_end_indices, start_pos );
        % Also test all contigs with length under max_contig_len_to_test
        del_starts = [ del_starts, contig_start_indices( contig_lengths <= max_contig_len_to_test ) ];
        del_ends = [ del_ends, contig_end_indices( contig_lengths <= max_contig_len_to_test ) ];
        % Filter candidate losss
        for j=1:numel(del_starts)
            if (del_ends(j) - del_starts(j) + 1) >= min_size_loss ... %length
                    && mean(all_coverage_per_bp_copynum(i,del_starts(j):del_ends(j))) <= max_cutoff_for_loss ... %normalized 
                    && mean(all_coverage_per_bp(i,del_starts(j):del_ends(j))) <= max_avg_cov_for_loss ...%raw data
                    && max( mean( all_coverage_per_bp_copynum(:,del_starts(j):del_ends(j)),2 ) ) >= min_avg_copynum_for_loss_high_control % require another sample to have the region              
                % Record candidate loss
                copy_number_variants_strains(end+1)=i;
                copy_number_variants_starts(end+1)=del_starts(j);
                copy_number_variants_ends(end+1)=del_ends(j);
                isdel(end+1)=1;
            end
        end
        
    end
    
end

fprintf(1,['Number of candidate gain/loss regions: ' num2str(numel(copy_number_variants_strains)) '.\n'] )

%% Generate clickable plot

% Convert back to genome position
if numel(copy_number_variants_strains)>0
    div_clickable_scatter_coverage( ...
        copy_number_variants_strains, copy_number_variants_starts+start_pos-1, copy_number_variants_ends+start_pos-1, isdel, ...
        SampleNames, all_coverage_per_bp, all_coverage_per_bp_copynum, cov_mat_doublenorm, ...
        this_clade_name, GenomeLength, ChrStarts, genome_contig_seqs, CDS, ...
        dir_save_regions )
end


%% Collapse overlapping regions across samples
fprintf(1,['Collapsing gains/loss regions...' '\n'] )

% Boolean matrix of all regions
all_regions_bool = zeros( num_samples, GenomeLength, 'logical' ); % initialize
start_pos_genome = copy_number_variants_starts+start_pos-1;
end_pos_genome = copy_number_variants_ends+start_pos-1; 
for r0=1:numel(start_pos_genome)
    all_regions_bool( copy_number_variants_strains(r0), start_pos_genome(r0):end_pos_genome(r0) ) = 1;
end
%check for samples that remain contaminated after contig filtering
all_regions_per_sample = sum(all_regions_bool,2);
potential_contaminants = find(all_regions_per_sample>500000); %do not expect gain loss regions longer than 0.5Mb

% Detect all collapsed regions
all_regions_bool_combined = sum(all_regions_bool);
all_regions_bool_combined( all_regions_bool_combined>1 ) = 1;
all_regions_starts_0 = find( diff(all_regions_bool_combined) > 0 ) + 1;
all_regions_ends_0 = find( diff(all_regions_bool_combined) < 0 );
if all_regions_bool_combined(end)==1 % case where region includes last position on genome
    all_regions_ends_0 = [ all_regions_ends_0, GenomeLength ];
end
% Split across contig boundaries if regions on adjacent contigs happened to merge
[ all_regions_starts, all_regions_ends ] = break_at_contig_boundaries( all_regions_starts_0, all_regions_ends_0, ...
    contig_num_by_index, contig_start_indices, contig_end_indices, 1 );

% Print number of total regions
num_regions_final = numel(all_regions_starts);
fprintf(1,['Number of collapsed gain/loss regions: ' num2str(num_regions_final) '.\n'] )


%% Save info for all regions found
fprintf(1,['Recording info on gains/loss regions...' '\n'] )


% Make a table make a info on all regions
% Make a fasta with all region sequences
% Make a summary plot for each region

% Initialize
all_regions_table = struct;

% Loop through regions
for r=1:numel(all_regions_starts)
    
    % Basic info on region
    region_name = [ this_clade_name '-reg-' num2str(r) ];
    start_index = all_regions_starts(r);
    end_index = all_regions_ends(r);
    contig_index = contig_num_by_index(start_index);
    if contig_index > numel(CDS)
        continue %no prokka annotations in region
    end
    start_contig_pos = p2chrpos(start_index,ChrStarts);
    end_contig_pos = p2chrpos(end_index,ChrStarts);
    
    % Make a table
    % Names
    all_regions_table(r).Name = region_name;
    all_regions_table(r).Cluster = this_clade_name;
    % Compute how many colonies have this region present
    bool_samples_pos = ( mean( all_coverage_per_bp_copynum(:,start_index:end_index),2 ) >= min_cutoff_for_gain );
    bool_samples_neg = ( mean( all_coverage_per_bp_copynum(:,start_index:end_index),2 ) <= max_cutoff_for_loss );
    bool_samples_ambig = ~bool_samples_pos & ~bool_samples_neg; 
    num_samples_pos = sum( bool_samples_pos );
    num_samples_neg = sum( bool_samples_neg );
    all_regions_table(r).NumColonies_Pos = num_samples_pos;
    all_regions_table(r).NumColonies_Neg = num_samples_neg;
    all_regions_table(r).NumColonies_Ambig = num_samples - num_samples_pos - num_samples_neg;
    all_regions_table(r).NamesColonies_Pos = SampleNames( bool_samples_pos );
    all_regions_table(r).NamesColonies_Neg = SampleNames( bool_samples_neg );
    all_regions_table(r).NamesColonies_Ambig = SampleNames( bool_samples_ambig );
    % Length and position
    all_regions_table(r).Region_Length = end_index-start_index+1;
    all_regions_table(r).Genome_IndexStart = start_index;
    all_regions_table(r).Genome_IndexEnd = end_index;
    all_regions_table(r).Contig_Num = contig_index;
    all_regions_table(r).Contig_Length = contig_lengths( contig_index );
    all_regions_table(r).Contig_PosStart = start_contig_pos(2);
    all_regions_table(r).Contig_PosEnd = end_contig_pos(2);
    % Gene content
    contig_seq = genome_contig_seqs{ contig_index };
    region_seq = contig_seq( start_contig_pos(2):end_contig_pos(2) );
    all_regions_table(r).Sequence = region_seq;
    [ genes_annotations, genes_translations, genes_locustags ] = get_gene_annnotations( all_regions_starts(r), all_regions_ends(r), ChrStarts, CDS );
    all_regions_table(r).Gene_Annotations = genes_annotations;
    all_regions_table(r).Gene_AAs = genes_translations;
    all_regions_table(r).Gene_locustags = genes_locustags;

    % Make a fasta file
    fasta_filename = [ dir_save_regions '/' region_name '.fasta' ];
    if exist( fasta_filename, 'file' ) % remove if already exists
        delete( fasta_filename )
    end
    fastawrite( fasta_filename, region_name, region_seq );

    % Make a summary plot
    f1=figure('Position',[0 0 1200 800],'visible','off'); clf(f1); hold on;
    if verLessThan('matlab','9.5')
      st=suptitle([ 'region: ' region_name '  length: ' num2str(end_index-start_index+1)] );
    else
      st=sgtitle([ 'region: ' region_name '  length: ' num2str(end_index-start_index+1)] );
    end
    st.Interpreter = 'none';
    st.FontSize=20;
    subplot_num = 3;
    color_pos = rgb('Blue');
    color_neg = rgb('Red');
    color_ambig = 0.15*[ 1 1 1 ] ;
    % 1. Absolute coverage plot
    subplot_index = 1;
    make_summary_coverage_plot( subplot_num, subplot_index, num_samples, ...
        bool_samples_pos, bool_samples_neg, bool_samples_ambig, ...
        all_coverage_per_bp_copynum, 'copy number', ...
        start_index, end_index, contig_start_indices, GenomeLength, ...
        color_pos, color_neg, color_ambig )
    % 2. Copy number coverage plot
    subplot_index = 2;
    make_summary_coverage_plot( subplot_num, subplot_index, num_samples, ...
        bool_samples_pos, bool_samples_neg, bool_samples_ambig, ...
        all_coverage_per_bp, 'raw coverage', ...
        start_index, end_index, contig_start_indices, GenomeLength, ...
        color_pos, color_neg, color_ambig )
    % 3. Normalized coverage plot
    subplot_index = 3;
    make_summary_coverage_plot( subplot_num, subplot_index, num_samples, ...
        bool_samples_pos, bool_samples_neg, bool_samples_ambig, ...
        cov_mat_doublenorm, 'normalized coverage', ...
        start_index, end_index, contig_start_indices, GenomeLength, ...
        color_pos, color_neg, color_ambig )
    % Save summary plot
    print([dir_save_regions '/' region_name '.png'],'-dpng')

end

%% Check for covariance to make metaregions

fprintf(1,['Forming metaregions from correlated regions...' '\n'] )
f1000=figure('Position',[0 0 1600 800],'visible','off'); clf(f1000); hold on;

% Make matrix of coverage for each region
region_cov=zeros(num_samples,num_regions_final);
num_metaregions_checked=0;
candidate_metaregions=zeros(num_regions_final);

for r=1:numel(all_regions_starts)
    % Basic info on region
    region_name = [ this_clade_name '-reg-' num2str(r) ];
    start_index = all_regions_starts(r);
    end_index = all_regions_ends(r);
    contig_index = contig_num_by_index(start_index);
    region_cov(:,r)=mean(all_coverage_per_bp_copynum(:,start_index:end_index),2);    
end

region_present=region_cov; %all samples, including ambiguous
region_present(region_present<=max_cutoff_for_loss)=0;
region_present(isnan(region_present))=0; %fix in future?
region_present(region_present>=min_cutoff_for_gain)=1;

max_copynum_deviation = 2; %fold
max_copynum_deviation = 1./max_copynum_deviation;

if numel(all_regions_starts)>1 %if multiple regions

    % plot presence heatmap
    subplot(2,2,1)
    h1 = heatmap(region_present);
    xlabel('Region')
    ylabel('Sample')
    h1.YDisplayLabels=SampleNames;
    title(['Region Presence: Copy Number < ' num2str(max_cutoff_for_loss) ' or > ' num2str(min_cutoff_for_gain)] )
    
    % plot coverage heatmap
    subplot(2,2,2)
    h2 = heatmap(round(region_cov,1));
    xlabel('Region')
    ylabel('Sample')
    h2.YDisplayLabels=SampleNames;
    title('Mean Region Coverage')
    
    % calculate correlation based on presence, slope based on coverage
    reg_corr=corrcoef(region_present);
    reg_slope=zeros(num_regions_final,num_regions_final);
    for m=1:num_regions_final
        for n=1:num_regions_final
            reg_slope(m,n)=region_cov(:,m)\region_cov(:,n);
            if reg_slope(m,n) > 1 %convert all to ratio
                reg_slope(m,n) = 1./reg_slope(m,n);
            end
        end
    end
    
    % plot dendrogram from presence
    subplot(2,2,3)
    presence_dist=pdist(region_present','correlation');
    presence_link = linkage(presence_dist);
    dendrogram(presence_link)
    title('Regions Clustered by Presence Correlation')
    % cluster with cutoff of 1 sample mismatch 
    candidate_metaregions = cluster(presence_link,'Criterion','Distance','cutoff',1/num_samples);
    num_metaregions=0;
    for m=1:numel(unique(candidate_metaregions))
        % singletons get assigned 0
        mr_idx = find(candidate_metaregions == m);
        if numel(mr_idx) > 1
            num_metaregions = num_metaregions + 1;
            candidate_metaregions(mr_idx) = num_metaregions;
        else
            candidate_metaregions(mr_idx) = 0;
        end
    end
    
    %relabel correlation heatmap
    reg_labels = {};
    for m=1:num_regions_final
        reg_labels{m}=[num2str(m) '   ' num2str(candidate_metaregions(m)) ''];
    end
    h1.XDisplayLabels=reg_labels;
    h1.NodeChildren(3).XTickLabelRotation = 90;
    

    % plot slope for each metaregion
    subplot(2,2,4)
    hold off
    if num_metaregions %if at least one metaregion  
        all_candidate_metaregions = unique(candidate_metaregions(candidate_metaregions>0));
        m=1;
        while m<=numel(all_candidate_metaregions) %iteratively split metaregions
            mr = all_candidate_metaregions(m);
            mr_idx = find(candidate_metaregions == mr);
            num_metaregions_checked = sum(all_candidate_metaregions>0);
            if numel(mr_idx)>1 && mr>0
                fprintf(1,['Checking metaregion ' num2str(mr) '.\n'] )
                mr_pairs = nchoosek(mr_idx,2); %need to convert to index
                mr_yvals = reg_slope(sub2ind(size(reg_slope),mr_pairs(:,1),mr_pairs(:,2)));
                mr_xvals = m.*ones(size(mr_yvals));
                scatter(mr_xvals,mr_yvals,'filled')
                hold on
                dx = 0.02;
                textlabel = num2str(mr_pairs);
                text(mr_xvals+dx,mr_yvals-dx,textlabel,'Fontsize', 8)

                % split cluster in half if metaregion has variable copy number
                pair_cutoff = (numel(mr_idx)-1); % number of deviating pairs if one region doesn't belong
                if sum(mr_yvals<max_copynum_deviation) >= pair_cutoff
                    metaregion_cov = region_cov(:,mr_idx);
                    split_metaregions = clusterdata(metaregion_cov','Maxclust',2);
                    candidate_metaregions(mr_idx)=mr.*10+split_metaregions;
                    all_candidate_metaregions = [all_candidate_metaregions; mr.*10+1; mr.*10+2];
                end
            elseif numel(mr_idx)==1
                candidate_metaregions(mr_idx)=0; %singletons get assigned 0
            end
            m = m+1;
        end
        ylim([0,1])
        xlim([0,num_metaregions_checked+1])
        xticks(1:num_metaregions_checked)
        xticklabels(num2str(all_candidate_metaregions(all_candidate_metaregions>0)))
 
    end
    plot([0,num_metaregions_checked+1],[max_copynum_deviation,max_copynum_deviation],'r--')
    xlabel('Metaregion')
    ylabel('Slope of Coverage Correlation')
    
    %relabel heatmaps
    reg_labels = {};
    for m=1:num_regions_final
        reg_labels{m}=[num2str(m) '   ' num2str(candidate_metaregions(m)) ''];
    end
    h2.XDisplayLabels=reg_labels;
    h2.NodeChildren(3).XTickLabelRotation = 90;
    print([dir_save_regions '/' this_clade_name '_metaregions.png'],'-dpng')
    metaregions_final = unique(candidate_metaregions(candidate_metaregions>0)); %no singletons
    num_metaregions_final = numel(metaregions_final); %no singletons
    
    for m=1:num_metaregions_final
        mr = metaregions_final(m);
        candidate_metaregions(candidate_metaregions == mr) = m;
    end
    
    fprintf(1,['Number of metaregions identified: ' num2str(num_metaregions_final) '.\n'] )
end


%% Save metaregion info to table and file

for r=1:numel(all_regions_starts)
    mr=candidate_metaregions(r);
    %save to table
    all_regions_table(r).Metaregion = mr;
    all_regions_table(r).per_sample_coverage = {region_cov(:,r)};
    
    if candidate_metaregions(r)~=0
        %save region to metaregion fasta
        metaregion_name = [ this_clade_name '-metareg-' num2str(mr) ];
        fasta_filename = [ dir_save_regions '/' metaregion_name '.fasta' ];
        %if first instance of metaregion
        if r==find(candidate_metaregions==mr,1, 'first')
            if exist( fasta_filename, 'file' ) % remove if already exists
                delete( fasta_filename )
            end
        end
        fastawrite(fasta_filename,all_regions_table(r).Name,all_regions_table(r).Sequence)
    end
    % save region to fasta
    region_name = [ this_clade_name '-reg-' num2str(r) ];
    fasta_filename = [ dir_save_regions '/' region_name '.fasta' ];
    if exist( fasta_filename, 'file' ) % remove if already exists
        delete( fasta_filename )
    end
    fastawrite(fasta_filename,all_regions_table(r).Name,all_regions_table(r).Sequence)
end


%%
% Save table
save( [ dir_save_regions '/' this_clade_name '_regions.mat' ], 'all_regions_table','potential_contaminants')
fprintf(1,['Done!' '\n\n\n'])
end

