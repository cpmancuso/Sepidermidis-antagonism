function sorted_structure = sort_interaction_structure(interaction_structure,sorting_vars)
lineage_order = [];
load("lineage_order.mat")
table = struct2table(interaction_structure.metadata);
if contains(sorting_vars,'Lineage')
    idxs = [];
    for n=1:numel(lineage_order)
        idxs = [idxs;find(table.Lineage==lineage_order(n))];
    end
    % Also grab other species
    other_species = unique(table.Lineage(table.Lineage>=200));
    if numel(other_species)
        for n=1:numel(other_species)
            idxs = [idxs;find(table.Lineage==other_species(n))];
        end
    end
else
    [~,idxs] = sortrows(table,sorting_vars);
end
sorted_structure = subsample_interaction_structure(interaction_structure,idxs);
