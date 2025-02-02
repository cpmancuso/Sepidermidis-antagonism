function merged_interaction_structure = merge_interaction_structure(interaction_structure1,interaction_structure2,merge_option)

% OPTION: In most cases, when merging interaction structures, one should
% use only isolates that grew in both conditions, but it may be necessary 
% for some plots to include isolates from only one condition 1.
merge_exclusion_data = true;


% merge_option: and, or, only1, only2

merged_interaction_structure = interaction_structure1;
num_spots = size(merged_interaction_structure.metadata,1);

switch merge_option
    case 'and'
        merged_interaction_structure.ZOI_depth = min(interaction_structure1.ZOI_depth(:,:),interaction_structure2.ZOI_depth(:,:));
        merged_interaction_structure.ZOI_AUC = min(interaction_structure1.ZOI_AUC(:,:),interaction_structure2.ZOI_AUC(:,:));
        merged_interaction_structure.ZOI_width = min(interaction_structure1.ZOI_width(:,:),interaction_structure2.ZOI_width(:,:));
        merged_interaction_structure.ZOI_call = min(interaction_structure1.ZOI_call(:,:),interaction_structure2.ZOI_call(:,:));  
    case 'or'
        merged_interaction_structure.ZOI_depth = max(interaction_structure1.ZOI_depth(:,:),interaction_structure2.ZOI_depth(:,:));
        merged_interaction_structure.ZOI_AUC = max(interaction_structure1.ZOI_AUC(:,:),interaction_structure2.ZOI_AUC(:,:));
        merged_interaction_structure.ZOI_width = max(interaction_structure1.ZOI_width(:,:),interaction_structure2.ZOI_width(:,:));
        merged_interaction_structure.ZOI_call = max(interaction_structure1.ZOI_call(:,:),interaction_structure2.ZOI_call(:,:));
    case 'only1'
        merged_interaction_structure.ZOI_depth = interaction_structure1.ZOI_depth(:,:);
        merged_interaction_structure.ZOI_AUC = interaction_structure1.ZOI_AUC(:,:);
        merged_interaction_structure.ZOI_width = interaction_structure1.ZOI_width(:,:);
        merged_interaction_structure.ZOI_call = interaction_structure1.ZOI_call(:,:);
    case 'only2'
        merged_interaction_structure.ZOI_depth = interaction_structure2.ZOI_depth(:,:);
        merged_interaction_structure.ZOI_AUC = interaction_structure2.ZOI_AUC(:,:);
        merged_interaction_structure.ZOI_width = interaction_structure2.ZOI_width(:,:);
        merged_interaction_structure.ZOI_call = interaction_structure2.ZOI_call(:,:);
end

% Propagate exclusion criteria from both experiments
for n=1:num_spots
    if ~strcmp(interaction_structure1.metadata(n).Trustworthy,"TRUE")
        if strcmp(interaction_structure1.metadata(n).Trustworthy,"replacement") %can't merge replacements
            merged_interaction_structure.metadata(n).Trustworthy = "was replaced";
        else
            merged_interaction_structure.metadata(n).Trustworthy = interaction_structure1.metadata(n).Trustworthy;
        end
    else
        merged_interaction_structure.metadata(n).Trustworthy = interaction_structure1.metadata(n).Trustworthy;
    end
    if ~strcmp(interaction_structure2.metadata(n).Trustworthy,"TRUE") & merge_exclusion_data
        if strcmp(interaction_structure2.metadata(n).Trustworthy,"replacement") %can't merge replacements
            merged_interaction_structure.metadata(n).Trustworthy = "was replaced";
        else
            merged_interaction_structure.metadata(n).Trustworthy = interaction_structure2.metadata(n).Trustworthy;
        end
    end
end

% Propagate ZOI data from both experiments
merged_interaction_structure.ZOI_depth1 = interaction_structure1.ZOI_depth(:,:);
merged_interaction_structure.ZOI_AUC1 = interaction_structure1.ZOI_AUC(:,:);
merged_interaction_structure.ZOI_width1 = interaction_structure1.ZOI_width(:,:);
merged_interaction_structure.ZOI_call1 = interaction_structure1.ZOI_call(:,:);
merged_interaction_structure.ZOI_depth2 = interaction_structure2.ZOI_depth(:,:);
merged_interaction_structure.ZOI_AUC2 = interaction_structure2.ZOI_AUC(:,:);
merged_interaction_structure.ZOI_width2 = interaction_structure2.ZOI_width(:,:);
merged_interaction_structure.ZOI_call2 = interaction_structure2.ZOI_call(:,:);

