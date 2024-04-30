function interactions_by_group = group_interactions(interaction_structure,groupby_option,minmax_option)

metadata_table = struct2table(interaction_structure.metadata);

%% Group samples and sort
interactions_by_group = struct();

switch groupby_option
    case 'isolate'
        %group by sample
        [group, id] = findgroups(metadata_table.Name);
    case 'lineage'
        % group by lineage
        [group, id] = findgroups(metadata_table.Lineage);
end

% create array of representative indexes for each group
rep_idxs = [];
for g =1:numel(id)
    idxs = find(group==g);
    switch minmax_option
        %choose representative member
        case 'min'
            num_interactions =  sum(interaction_structure.ZOI_call(idxs,:),2);
            [~,i] = min(num_interactions);
            rep_idxs(g) = idxs(i);
        case 'max'
            num_interactions =  sum(interaction_structure.ZOI_call(idxs,:),2);
            [~,i] = max(num_interactions);
            rep_idxs(g) = idxs(i);
    end
    % carry over AAI/ANI/Treedist info from groupmate if necessary
    if interaction_structure.metadata(rep_idxs(g)).ANI_index==0
        interaction_structure.metadata(rep_idxs(g)).ANI_index = max([interaction_structure.metadata(idxs).ANI_index]);
    end
    if interaction_structure.metadata(rep_idxs(g)).AAI_index==0
        interaction_structure.metadata(rep_idxs(g)).AAI_index = max([interaction_structure.metadata(idxs).AAI_index]);
    end
    if interaction_structure.metadata(rep_idxs(g)).Treedist_index==0
        interaction_structure.metadata(rep_idxs(g)).Treedist_index = max([interaction_structure.metadata(idxs).Treedist_index]);
    end

end

% Take representative metadata
interactions_by_group.metadata =  interaction_structure.metadata(rep_idxs);
% interactions_by_group.centers =  interaction_structure.centers(rep_idxs,rep_idxs);
% interactions_by_group.radii =  interaction_structure.radii(rep_idxs,rep_idxs);
% interactions_by_group.norm_int =  interaction_structure.norm_int(rep_idxs,rep_idxs,:); %3D
% interactions_by_group.norm_x =  interaction_structure.norm_x(rep_idxs,rep_idxs,:); %3D

% Group interaction data, complicated splitapply

interactions_by_group.ZOI_depth = zeros(numel(id));
interactions_by_group.ZOI_AUC = zeros(numel(id));
interactions_by_group.ZOI_width = zeros(numel(id));
interactions_by_group.ZOI_call = false(numel(id));

for g =1:numel(id)
    gidxs = find(group==g);
    for h=1:numel(id)        
        hidxs = find(group==h);
        
        if numel(gidxs)==numel(hidxs)&&numel(gidxs)==1
            interactions_by_group.ZOI_depth(g,h) = interaction_structure.ZOI_depth(gidxs,hidxs);
            interactions_by_group.ZOI_AUC(g,h) = interaction_structure.ZOI_AUC(gidxs,hidxs);
            interactions_by_group.ZOI_width(g,h) = interaction_structure.ZOI_width(gidxs,hidxs);
            interactions_by_group.ZOI_call(g,h) = interaction_structure.ZOI_call(gidxs,hidxs);
        else
            group_call = interaction_structure.ZOI_AUC(gidxs,hidxs).*interaction_structure.ZOI_call(gidxs,hidxs);
            % take min or max
            switch minmax_option
                case 'min'
                [~,idx] = min(group_call(:));
                case 'max'
                [~,idx] = max(group_call(:));
            end
            % convert linear index to local subscripts
            [row,col] = ind2sub([numel(gidxs) numel(hidxs)],idx);
            interactions_by_group.ZOI_depth(g,h) = interaction_structure.ZOI_depth(gidxs(row),hidxs(col));
            interactions_by_group.ZOI_AUC(g,h) = interaction_structure.ZOI_AUC(gidxs(row),hidxs(col));
            interactions_by_group.ZOI_width(g,h) = interaction_structure.ZOI_width(gidxs(row),hidxs(col));
            interactions_by_group.ZOI_call(g,h) = interaction_structure.ZOI_call(gidxs(row),hidxs(col));
        end
    end
end
    
