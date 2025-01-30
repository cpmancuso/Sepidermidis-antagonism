function subsample_structure = subsample_interaction_structure(interaction_structure,idxs_to_keep)

subsample_structure = struct();
subsample_structure.ZOI_depth = interaction_structure.ZOI_depth(idxs_to_keep,idxs_to_keep);
subsample_structure.ZOI_AUC = interaction_structure.ZOI_AUC(idxs_to_keep,idxs_to_keep);
subsample_structure.ZOI_width = interaction_structure.ZOI_width(idxs_to_keep,idxs_to_keep);
subsample_structure.ZOI_call = interaction_structure.ZOI_call(idxs_to_keep,idxs_to_keep);

% expects metadata in structure form
if any(contains(fieldnames(interaction_structure),'metadata'))
    subsample_structure.metadata = interaction_structure.metadata(idxs_to_keep);
end

% copy image specific information over only if needed
if any(contains(fieldnames(interaction_structure),'centers'))
    subsample_structure.centers = interaction_structure.centers(idxs_to_keep,idxs_to_keep,:); %3D
    subsample_structure.radii = interaction_structure.radii(idxs_to_keep,idxs_to_keep);
    subsample_structure.norm_int = interaction_structure.norm_int(idxs_to_keep,idxs_to_keep,:); %3D
    subsample_structure.norm_x = interaction_structure.norm_x(idxs_to_keep,idxs_to_keep,:); %3D
end

% retain ZOI calls from multiple experiments, if exist, needed for some
% analyses

if any(contains(fieldnames(interaction_structure),'ZOI_call1'))
    subsample_structure.ZOI_depth1 = interaction_structure.ZOI_depth1(idxs_to_keep,idxs_to_keep);
    subsample_structure.ZOI_AUC1 = interaction_structure.ZOI_AUC1(idxs_to_keep,idxs_to_keep);
    subsample_structure.ZOI_width1 = interaction_structure.ZOI_width1(idxs_to_keep,idxs_to_keep);
    subsample_structure.ZOI_call1 = interaction_structure.ZOI_call1(idxs_to_keep,idxs_to_keep);
    subsample_structure.ZOI_depth2 = interaction_structure.ZOI_depth2(idxs_to_keep,idxs_to_keep);
    subsample_structure.ZOI_AUC2 = interaction_structure.ZOI_AUC2(idxs_to_keep,idxs_to_keep);
    subsample_structure.ZOI_width2 = interaction_structure.ZOI_width2(idxs_to_keep,idxs_to_keep);
    subsample_structure.ZOI_call2 = interaction_structure.ZOI_call2(idxs_to_keep,idxs_to_keep);
end