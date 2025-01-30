function fighandle = plot_interaction_clustergram(interaction_structure,fignum)

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 1000 800];

ZOI_call = double(interaction_structure.ZOI_call);
labels = [interaction_structure.metadata.Lineage];
labels = cellstr(num2str(labels'));
% convert lineage labels to species for non S. epidermidis
for l = 1:numel(labels)
    if interaction_structure.metadata(l).Lineage>200
        newlabel = char(interaction_structure.metadata(l).Species);
        labels{l} = newlabel;
    end
end

cgo = clustergram(ZOI_call','Colormap',[0 0 0; 0 0 0; 1 1 1],'ColumnLabels',labels,'RowLabels',labels);
plot(cgo,fighandle)
