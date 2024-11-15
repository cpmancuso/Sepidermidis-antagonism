function [interaction_structure_with_identity_out,fighandle] = plot_AAI_ANI(interaction_structure_with_identity,AAI_filename,ANI_filename, fignum)

ANI = [];
AAI = [];

load(ANI_filename);
load(AAI_filename);

ANI_matrix = [];
num = numel(interaction_structure_with_identity.metadata);
for n=1:num
    for i=1:num
        % If AAI or ANI values exist, keep them
        try 
            AAI_matrix(n,i)=AAI(interaction_structure_with_identity.metadata(n).AAI_index,interaction_structure_with_identity.metadata(i).AAI_index);
        catch
            AAI_matrix(n,i)=0;
        end
        try
            ANI_matrix(n,i)=ANI(interaction_structure_with_identity.metadata(n).ANI_index,interaction_structure_with_identity.metadata(i).ANI_index);
        catch
            ANI_matrix(n,i)=0;
        end
    end
end
interaction_structure_with_identity.AAI_matrix = AAI_matrix;
interaction_structure_with_identity.ANI_matrix = ANI_matrix;

% filter by call
ZOI_call = interaction_structure_with_identity.ZOI_call;
ZOI_AUC = (interaction_structure_with_identity.ZOI_AUC).*ZOI_call;
ANI_matrix = interaction_structure_with_identity.ANI_matrix;

% bin data
binsize = 0.0025;
AAI_lim = 0.45;
ANI_lim = 0.75;
AAI_binedge=[AAI_lim:binsize:1];
ANI_binedge=[ANI_lim:binsize:1];
AAI_calls_per_bin = cell(size(AAI_binedge));
ANI_calls_per_bin = cell(size(ANI_binedge));
AAI_AUC_per_bin = cell(size(AAI_binedge));
ANI_AUC_per_bin = cell(size(ANI_binedge));

% Assign data to each bin
for n=1:num
    for i=1:num
        AAI_idx=floor((AAI_matrix(n,i)-AAI_lim)/binsize);
        ANI_idx=floor((ANI_matrix(n,i)-ANI_lim)/binsize);
        if ZOI_call(n,i)
            if AAI_idx>0
                AAI_calls_per_bin{AAI_idx}=[AAI_calls_per_bin{AAI_idx} 1];
                AAI_AUC_per_bin{AAI_idx}=[AAI_AUC_per_bin{AAI_idx} ZOI_AUC(n,i)];
            end
            if ANI_idx>0
                ANI_calls_per_bin{ANI_idx}=[ANI_calls_per_bin{ANI_idx} 1];
                ANI_AUC_per_bin{ANI_idx}=[ANI_AUC_per_bin{ANI_idx} ZOI_AUC(n,i)];
            end
        else
            if AAI_idx>0
                AAI_calls_per_bin{AAI_idx}=[AAI_calls_per_bin{AAI_idx} 0];
            end
            if ANI_idx>0
                ANI_calls_per_bin{ANI_idx}=[ANI_calls_per_bin{ANI_idx} 0];
            end
        end
    end
end

AAI_bin_interaction_frequency = [];
AAI_bin_weight = [];
AAI_AUC_median = [];
AAI_AUC_weight = [];
ANI_bin_interaction_frequency = [];
ANI_bin_weight = [];
ANI_AUC_median = [];
ANI_AUC_weight = [];

for b=1:numel(AAI_calls_per_bin)
    AAI_bin_interaction_frequency(b) = mean(AAI_calls_per_bin{b});
    AAI_bin_weight(b) = numel(AAI_calls_per_bin{b});
    AAI_AUC_median(b) = median(AAI_AUC_per_bin{b});
    AAI_AUC_weight(b) = numel(AAI_AUC_per_bin{b});
end

for b=1:numel(ANI_calls_per_bin)
    ANI_bin_interaction_frequency(b) = mean(ANI_calls_per_bin{b});
    ANI_bin_weight(b) = numel(ANI_calls_per_bin{b});
    ANI_AUC_median(b) = median(ANI_AUC_per_bin{b});
    ANI_AUC_weight(b) = numel(ANI_AUC_per_bin{b});
end

% remove empty bins
% AAI_binedge = AAI_binedge(AAI_bin_weight>0);
% AAI_bin_interaction_frequency = AAI_bin_interaction_frequency(AAI_bin_weight>0);
% AAI_AUC_median = AAI_AUC_median(AAI_bin_weight>0);
% AAI_bin_weight = AAI_bin_weight(AAI_bin_weight>0);

% ANI_binedge = ANI_binedge(ANI_bin_weight>0);
% ANI_bin_interaction_frequency = ANI_bin_interaction_frequency(ANI_bin_weight>0);
% ANI_AUC_median = ANI_AUC_median(ANI_bin_weight>0);
% ANI_bin_weight = ANI_bin_weight(ANI_bin_weight>0);


%% Plot interaction frequency results
fighandle = figure(fignum);
fighandle.Position = [100 100 1000 700];
clf(fighandle)
subplot(2,3,1)

% AAI vs Interaction Frequency
subplot(2,3,1)
[rsq,pval] = plot_interaction_vs_identity(AAI_binedge,AAI_bin_interaction_frequency,AAI_bin_weight);
xlim([AAI_lim 1])
xlabel('Average Amino Acid Identity')


% ANI vs Interaction Frequency, interspecies
subplot(2,3,2)
[rsq,pval] = plot_interaction_vs_identity(ANI_binedge,ANI_bin_interaction_frequency,ANI_bin_weight);
xlim([ANI_lim 1])
xlabel('Average Nucleotide Identity')

% ANI vs Interaction Frequency, interspecies
subplot(2,3,3)
intra_idxs = (ANI_binedge > 0.95);
[rsq,pval] = plot_interaction_vs_identity(ANI_binedge(intra_idxs),ANI_bin_interaction_frequency(intra_idxs),ANI_bin_weight(intra_idxs));
xlim([0.95 1])
xlabel('Average Nucleotide Identity')

%% Plot interaction AUC results
lw = 1;
% Note that because some bins had all 0 interactions, can't take median
% Interspecies ZOI size
subplot(2,3,4)
scatter_AAI = reshape(AAI_matrix,1,[]);
scatter_ZOI = reshape(ZOI_AUC,1,[]);
idxs = scatter_AAI>AAI_lim & scatter_ZOI>0;

scatter_AAI = scatter_AAI(idxs);
scatter_ZOI = scatter_ZOI(idxs);

[rsq,pval] = plot_AUC_vs_identity(scatter_AAI,scatter_ZOI,AAI_AUC_median(~isnan(AAI_AUC_median)),AAI_AUC_weight(~isnan(AAI_AUC_median)));
xlim([AAI_lim 1])
xlabel('Average Amino Acid Identity')
ylabel('Area Under Curve of ZOI')

% Interspecies ZOI size
subplot(2,3,5)
scatter_ANI = reshape(ANI_matrix,1,[]);
scatter_ZOI = reshape(ZOI_AUC,1,[]);
idxs = scatter_ANI>ANI_lim & scatter_ZOI>0;

scatter_ANI = scatter_ANI(idxs);
scatter_ZOI = scatter_ZOI(idxs);

[rsq,pval] = plot_AUC_vs_identity(scatter_ANI,scatter_ZOI,ANI_AUC_median(~isnan(ANI_AUC_median)),ANI_AUC_weight(~isnan(ANI_AUC_median)));
xlim([ANI_lim 1])
xlabel('Average Nucleotide Identity')
ylabel('Area Under Curve of ZOI')

% Interspecies ZOI size
subplot(2,3,6)
scatter_ANI = reshape(ANI_matrix,1,[]);
scatter_ZOI = reshape(ZOI_AUC,1,[]);
idxs = scatter_ANI>0.95 & scatter_ZOI>0;

scatter_ANI = scatter_ANI(idxs);
scatter_ZOI = scatter_ZOI(idxs);
[rsq,pval] = plot_AUC_vs_identity(scatter_ANI,scatter_ZOI,ANI_AUC_median(intra_idxs&~isnan(ANI_AUC_median)),ANI_AUC_weight(intra_idxs&~isnan(ANI_AUC_median)));
xlim([0.95 1])
xlabel('Average Nucleotide Identity')
ylabel('Area Under Curve of ZOI')

interaction_structure_with_identity_out = interaction_structure_with_identity;
end

%% Subfunctions
function [rsq,pval] = plot_interaction_vs_identity(binedge,bin_interaction_frequency,binweight)
    formatSpec = '%.3f';
    mdl = fitlm(binedge,bin_interaction_frequency,'Weights',binweight);
    tbl = anova(mdl);
    h = plot(mdl);
    delete(h(1))
    hold on
    mksz=10*round(sqrt(binweight));
    mksz(mksz==0)=NaN;
    scatter(binedge,bin_interaction_frequency,mksz,'k','filled','MarkerFaceAlpha',0.8)
    scatter([1 1 1],[0.3,0.4,0.5],10*round(sqrt([1 10 100])),'r','filled','MarkerFaceAlpha',0.8) %for scale
    rsq = mdl.Rsquared.Ordinary;
    pval = tbl.pValue(1);
    title(['R-squared: ' num2str(rsq,formatSpec) '  p-value ' num2str(pval,formatSpec)])
    % legend('Linear Fit','Confidence Bounds','','Data')
    legend('off')

    ylabel('Frequency of Antagonism')
    ylim([0 0.6])
end

function [rsq,pval] = plot_AUC_vs_identity(scatter_identity,scatter_ZOI,AUC_median,AUC_median_weight)
    formatSpec = '%.3f';
    scatter(scatter_identity,scatter_ZOI,'k','filled','MarkerFaceAlpha',0.5)
    hold on
    mdl = fitlm(scatter_identity,scatter_ZOI);
    tbl = anova(mdl);
    h = plot(mdl);
    delete(h(1))
    rsq = mdl.Rsquared.Ordinary;
    pval = tbl.pValue(1);
    title(['R-squared: ' num2str(rsq,formatSpec) '  p-value ' num2str(pval,formatSpec)])
    legend('off')
    % plot(AAI_binedge(~isnan(AAI_AUC_median)),AAI_AUC_median(~isnan(AAI_AUC_median)),'r','LineWidth',lw)
end

