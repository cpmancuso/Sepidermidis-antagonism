function expanded_interaction_structure = expand_interaction_structure(interaction_structure,replacement_data_filename)

% NOTE: this function only replaces the interaction data, not the metadata.
% Metadata should be loaded normally using the updated information for
% replacements.

% Generate list of replacements
replacement_idxs = [];

for n=1:numel(interaction_structure.metadata)
    if strcmp(interaction_structure.metadata(n).Trustworthy,'replacement')
        replacement_idxs = [replacement_idxs, n];
    end
end

% disp(replacement_idxs)
% disp([interaction_structure.metadata(replacement_idxs).Lineage])

% Load replacement data
load(replacement_data_filename) %automated_ZOI_results
[nspot_data,nlawn_data] = size(automated_ZOI_results.ZOI_subarea);
[sz,~] = size(interaction_structure.ZOI_AUC);

if nlawn_data<sz
    % assume that interaction structure contains both datasets
    automated_ZOI_results.all_adj_centers = [automated_ZOI_results.all_adj_centers(1:numel(replacement_idxs),:,:), automated_ZOI_results.all_adj_centers(numel(replacement_idxs)+1:end,:,:)];
    automated_ZOI_results.all_radii = [automated_ZOI_results.all_radii(1:numel(replacement_idxs),:), automated_ZOI_results.all_radii(numel(replacement_idxs)+1:end,:)];
    automated_ZOI_results.norm_int_local = [automated_ZOI_results.norm_int_local(1:numel(replacement_idxs),:,:), automated_ZOI_results.norm_int_local(numel(replacement_idxs)+1:end,:,:)];
    automated_ZOI_results.norm_x = [automated_ZOI_results.norm_x(1:numel(replacement_idxs),:,:), automated_ZOI_results.norm_x(numel(replacement_idxs)+1:end,:,:)];
    automated_ZOI_results.ZOI_depth = [automated_ZOI_results.ZOI_depth(1:numel(replacement_idxs),:), automated_ZOI_results.ZOI_depth(numel(replacement_idxs)+1:end,:)];
    automated_ZOI_results.ZOI_subarea = [automated_ZOI_results.ZOI_subarea(1:numel(replacement_idxs),:), automated_ZOI_results.ZOI_subarea(numel(replacement_idxs)+1:end,:)];
    automated_ZOI_results.ZOI_width = [automated_ZOI_results.ZOI_width(1:numel(replacement_idxs),:), automated_ZOI_results.ZOI_width(numel(replacement_idxs)+1:end,:)];
end
[nspot_data,nlawn_data] = size(automated_ZOI_results.ZOI_subarea);

% disp([nspot_data,nlawn_data])

% Replace spot data with nan
for n=1:numel(replacement_idxs)
    idx = replacement_idxs(n);

    % replace 
    interaction_structure.centers(:,idx,:) = NaN;
    interaction_structure.radii(:,idx) = NaN;
    interaction_structure.norm_int(:,idx,:) = NaN; 
    interaction_structure.norm_x(:,idx,:) = NaN;
    interaction_structure.ZOI_depth(:,idx) = NaN;
    interaction_structure.ZOI_AUC(:,idx) = NaN;
    interaction_structure.ZOI_width(:,idx) = NaN;
end

% Replace lawn data with collected data
for n=1:numel(replacement_idxs)
    idx = replacement_idxs(n);

    % replace with 
    interaction_structure.centers(idx,:,:) = automated_ZOI_results.all_adj_centers(n,:,:);
    interaction_structure.radii(idx,:) = automated_ZOI_results.all_radii(n,:);
    interaction_structure.norm_int(idx,:,:) = automated_ZOI_results.norm_int_local(n,:,:);
    interaction_structure.norm_x(idx,:,:) = automated_ZOI_results.norm_x(n,:,:);
    interaction_structure.ZOI_depth(idx,:) = automated_ZOI_results.ZOI_depth(n,:);
    interaction_structure.ZOI_AUC(idx,:) = automated_ZOI_results.ZOI_subarea(n,:);
    interaction_structure.ZOI_width(idx,:) = automated_ZOI_results.ZOI_width(n,:);
end

expanded_interaction_structure = interaction_structure;