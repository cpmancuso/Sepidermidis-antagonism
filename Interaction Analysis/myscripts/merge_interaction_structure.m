function merged_interaction_structure = merge_interaction_structure(interaction_structure,additional_interaction_structure,merge_option)

merged_interaction_structure = interaction_structure;


% calculate correlation between interaction structures
ZOI_jaccard = pdist([interaction_structure.ZOI_call(:),additional_interaction_structure.ZOI_call(:)]','jaccard');
disp(['Jaccard similarity is: ' num2str(1-ZOI_jaccard)])

ZOI_SMC = sum(interaction_structure.ZOI_call(:)==additional_interaction_structure.ZOI_call(:))./numel(interaction_structure.ZOI_call(:));
disp(['Simple matching coefficient is: ' num2str(ZOI_SMC)])

switch merge_option
    case 'and'
        merged_interaction_structure.ZOI_depth = min(interaction_structure.ZOI_depth(:,:),additional_interaction_structure.ZOI_depth(:,:));
        merged_interaction_structure.ZOI_AUC = min(interaction_structure.ZOI_AUC(:,:),additional_interaction_structure.ZOI_AUC(:,:));
        merged_interaction_structure.ZOI_width = min(interaction_structure.ZOI_width(:,:),additional_interaction_structure.ZOI_width(:,:));
        merged_interaction_structure.ZOI_call = min(interaction_structure.ZOI_call(:,:),additional_interaction_structure.ZOI_call(:,:));
    case 'or'
        merged_interaction_structure.ZOI_depth = max(interaction_structure.ZOI_depth(:,:),additional_interaction_structure.ZOI_depth(:,:));
        merged_interaction_structure.ZOI_AUC = max(interaction_structure.ZOI_AUC(:,:),additional_interaction_structure.ZOI_AUC(:,:));
        merged_interaction_structure.ZOI_width = max(interaction_structure.ZOI_width(:,:),additional_interaction_structure.ZOI_width(:,:));
        merged_interaction_structure.ZOI_call = max(interaction_structure.ZOI_call(:,:),additional_interaction_structure.ZOI_call(:,:));
end

