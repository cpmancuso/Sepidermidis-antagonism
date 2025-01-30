function export_interaction_table(interaction_structure,label_option,filename)

labels = [interaction_structure.metadata.(label_option)]; 
matrix = interaction_structure.ZOI_call;

% Convert numeric labels to strings
label_strings = arrayfun(@num2str, labels, 'UniformOutput', false);

% Create a cell array for CSV export
csv_data = [[' ', label_strings]; [label_strings', num2cell(matrix)]];

% Write to CSV
writecell(csv_data, filename);