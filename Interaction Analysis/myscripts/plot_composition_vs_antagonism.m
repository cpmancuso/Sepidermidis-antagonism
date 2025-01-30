function plot_composition_vs_antagonism(composition_table,interactions_by_lineage,fignum)

interactions_by_lineage = sort_interaction_structure(interactions_by_lineage,{'Lineage'});

in_collection = false(1);
for n = 1:size(composition_table,2)-3
    col_name = ['Lineage_' num2str(n)];
    % in trustworthy portion of collection
    in_collection(n) = ismember(n,[interactions_by_lineage.metadata.Lineage]);
end

%% Composition plots per family with metagenomes
for f=1:8
    fam_idxs_i = ([composition_table.Family{:}]==f);

    if any(fam_idxs_i)
        fighandle = figure(f);
        clf(f)
        fighandle.Position=[fighandle.Position(1) fighandle.Position(2) 400 375];
        % Make subtable with both
        family_table = [composition_table(fam_idxs_i,:)];
        family_lineages = find(sum(family_table{:,4:end}));
        
        % sort by SID
        family_table = sortrows(family_table,'SID');

        % sort by in collection
        present = ~(in_collection(family_lineages));
        [~,I] = sort(present);
        family_lineages = family_lineages(I);

        family_lineage_labels = family_table.Properties.VariableNames(family_lineages+3);
        family_lineage_labels = strrep(family_lineage_labels, '_', ' ');
        composition_matrix = family_table{:,family_lineages+3}; %header columns
        
        % make color scheme with greys for lineages not in collection
        
        make_grey = ~(in_collection(family_lineages)); 
        color_vals = ones(numel(make_grey),3);
        for n=1:numel(family_lineages)
            if make_grey(n)
                color_vals(n,:) = [0 0 0];
            else
                this_lineage_idx = find([interactions_by_lineage.metadata.Lineage] == family_lineages(n));
                ant_freq = sum(interactions_by_lineage.ZOI_call(:,this_lineage_idx))./size(interactions_by_lineage.ZOI_call,2);
                if ant_freq<=mean(sum(interactions_by_lineage.ZOI_call))./size(interactions_by_lineage.ZOI_call,2);
                    color_vals(n,:) = [1 1 0];
                else
                    color_vals(n,:) = [ant_freq 0 1-ant_freq];
                end
            end
        end
        color_vals(make_grey,:) = ones(size(color_vals(make_grey,:))).*linspace(1,0.7,sum(make_grey))';
        colororder(color_vals)
        
        % plot
        bar(composition_matrix,'Stacked')
        legend(family_lineage_labels,'Interpreter','none','Location','eastoutside')
        xticklabels(family_table.SID)
        xtickangle(90)
        ylim([0 1])
        xlim([0 13])
        ylabel('Lineage Abundance')
    end
end

