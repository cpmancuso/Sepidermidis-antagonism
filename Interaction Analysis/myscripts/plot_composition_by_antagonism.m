function plot_composition_by_antagonism(composition_table,interactions_by_lineage,fignum)

color_option = 'norm'; %rank or freq or norm

% rank lineages by antagonism 
ant_freq = [];
for n=1:size(interactions_by_lineage.ZOI_call,2)
    ant_freq(n) = sum(interactions_by_lineage.ZOI_call(:,n))./size(interactions_by_lineage.ZOI_call,2);
end
ant_rank = tiedrank(ant_freq);
ant_rank = (ant_rank-min(ant_rank))./(max(ant_rank)-min(ant_rank));


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
        fighandle = figure(fignum+f);
        clf(fignum+f)
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
                color_vals(n,:) = [0.5 0.5 0.5];
            else
                this_lineage_idx = find([interactions_by_lineage.metadata.Lineage] == family_lineages(n));
                switch color_option
                    case 'rank'
                        color_vals(n,:) = [ant_rank(this_lineage_idx), 0, 1-ant_rank(this_lineage_idx)];
                        custom_cmap = [linspace(0, 1, 256)', linspace(0, 0, 256)', linspace(1, 0, 256)'];
                    case 'freq'
                        color_vals(n,:) = [ant_freq(this_lineage_idx), 0, 1-ant_freq(this_lineage_idx)];
                        custom_cmap = [linspace(0, 1, 256)', linspace(0, 0, 256)', linspace(1, 0, 256)'];
                    case 'norm'
                        norm_ant_freq = sqrt(ant_freq);
                        color_vals(n,:) = [norm_ant_freq(this_lineage_idx), 0, 1-norm_ant_freq(this_lineage_idx)];
                        custom_cmap = [sqrt(linspace(0, 1, 256)'), linspace(0, 0, 256)', 1-sqrt(linspace(0, 1, 256)')];
                end

            end
        end
        colororder(color_vals)
        
        % plot
        bar(composition_matrix,'Stacked')

        

        % Add text labels
        x=1:size(composition_table,1);
        for i = 1:size(composition_matrix, 1)       % Loop through each bar group
            cumulativeHeight = 0;     % Initialize cumulative height for stacked bars
            for j = 1:size(composition_matrix, 2)   % Loop through each stack
                % Calculate the position for the label
                y = cumulativeHeight + composition_matrix(i, j) / 2; % Middle of the stack
                cumulativeHeight = cumulativeHeight + composition_matrix(i, j); % Update cumulative height
                
                % Add text
                if composition_matrix(i, j) > 0.05     % Only add labels for sufficiently large values
                    text(x(i), y, num2str(family_lineages(j)), 'HorizontalAlignment', 'center','Color','w','FontSize',8);
                end
            end
        end

        legend(family_lineage_labels,'Interpreter','none','Location','eastoutside')
        xticklabels(family_table.SID)
        xtickangle(90)
        ylim([0 1])
        xlim([0 13])
        ylabel('Lineage Abundance')
        

    end

end
% Apply the custom colormap and plot
figure(fignum)
clf(fignum)
colormap(custom_cmap);
colorbar
