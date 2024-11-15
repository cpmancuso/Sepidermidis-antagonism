function read_and_plot_composition(data_filename, header_filename, MG_filename, plot_summary_option,fignum)
% plot_summary_option: True or False, recommend false. Will generate
% summary graphs describing the number of colonies in unused samples.

%% Read composition from datafiles

data = readtable(data_filename,'Delimiter',',');

num = size(data,1);
for n=1:num
    data.Family{n} = char(data.SampleName{n}(1));
    data.Subject{n} = char(data.SampleName{n}(1:3));
    data.SID{n} = char(data.SampleName{n}(1:4));
end

%% Filter out low diversity lineages

diversity_cutoff = 2;

sepi_isolates = data((data.SepiClusters_new~=0)&(strcmp(data.is_singleton,'FALSE')),[3,8,11:13]);
sepi_isolates.Properties.VariableNames{2} = 'Lineage';
for n=1:numel(sepi_isolates.Family)
    sepi_isolates.diverse_lineage(n) = numel(find(sepi_isolates.Lineage==sepi_isolates.Lineage(n)))>diversity_cutoff;
end
sepi_isolates = sepi_isolates(sepi_isolates.diverse_lineage,:);

%% Filter out low coverage samples

isolate_cutoff = 20;
family_isolate_cutoff = 100;

[G,SIDs] = findgroups(sepi_isolates.SID);
num_SID = splitapply(@numel,sepi_isolates.SampleName,G);

[G,subs] = findgroups(sepi_isolates.Subject);
num_Subject = splitapply(@numel,sepi_isolates.SampleName,G);

[G,fams] = findgroups(sepi_isolates.Family);
num_Family = splitapply(@numel,sepi_isolates.SampleName,G);

summary_table = table(SIDs);
num=size(summary_table,1);
for n=1:num
    summary_table.Family{n} = char(summary_table.SIDs{n}(1));
    summary_table.num_Family(n) = num_Family(find(strcmp(fams,summary_table.Family{n})));
    summary_table.Subject{n} = char(summary_table.SIDs{n}(1:3));
    summary_table.num_Subject(n) = num_Subject(find(strcmp(subs,summary_table.Subject{n})));
    summary_table.SID{n} = char(summary_table.SIDs{n}(1:4));
    summary_table.num_SID(n) = num_SID(find(strcmp(SIDs,summary_table.SID{n})));
    
end
summary_table = removevars(summary_table, 'SIDs');

good_coverage_SIDs = {summary_table.SID{(summary_table.num_SID>isolate_cutoff)&(summary_table.num_Family>family_isolate_cutoff)}}';
for n=1:numel(sepi_isolates.Family)
    sepi_isolates.good_coverage(n) = any(strcmp(good_coverage_SIDs,sepi_isolates.SID(n)));
end

sepi_isolates = sepi_isolates(sepi_isolates.good_coverage,:);

%% Compile composition statistics
lineages = 1:92;
composition_table = table();
for s = 1:numel(good_coverage_SIDs)
    sid = good_coverage_SIDs{s};
    idxs = find(strcmp(sepi_isolates.SID,sid));
    composition_table.SID{s} = sid;
    composition_table.Family{s} = char(sid(1));
    composition_table.Subject{s} = char(sid(1:3));
    for n = 1:max(sepi_isolates.Lineage)
        col_name = ['Lineage_' num2str(n)];
        composition_table.(col_name)(s) = sum(sepi_isolates.Lineage(idxs) == n);
    end
end

%%
% Print empty columns
comp_matrix = composition_table{:,4:95};
temp = sum(comp_matrix);
extra_lineages = composition_table.Properties.VariableNames{3+(find(temp==0))};

% List missing lineages from strain compettition experiment
headerdata = readtable(header_filename);
in_collection = false(1);
for n = 1:max(sepi_isolates.Lineage)
    col_name = ['Lineage_' num2str(n)];
    % in trustworthy portion of collection
    in_collection(n) = sum(headerdata.Lineage(strcmp(headerdata.Trustworthy,'TRUE')) == n);
end

%% Recalculate summary stats

[G,SIDs] = findgroups(sepi_isolates.SID);
num_SID = splitapply(@numel,sepi_isolates.SampleName,G);

[G,subs] = findgroups(sepi_isolates.Subject);
num_Subject = splitapply(@numel,sepi_isolates.SampleName,G);

[G,fams] = findgroups(sepi_isolates.Family);
num_Family = splitapply(@numel,sepi_isolates.SampleName,G);

summary_table = table(good_coverage_SIDs);
num=size(summary_table,1);
for n=1:num
    summary_table.Family{n} = char(summary_table.good_coverage_SIDs{n}(1));
    summary_table.num_Family(n) = num_Family(find(strcmp(fams,summary_table.Family{n})));
    summary_table.Subject{n} = char(summary_table.good_coverage_SIDs{n}(1:3));
    summary_table.num_Subject(n) = num_Subject(find(strcmp(subs,summary_table.Subject{n})));
    summary_table.SID{n} = char(summary_table.good_coverage_SIDs{n}(1:4));
    summary_table.num_SID(n) = num_SID(find(strcmp(SIDs,summary_table.SID{n})));
    
end
summary_table = removevars(summary_table, 'good_coverage_SIDs');

%% Calculate collection coverage percent based on lineage
comp_matrix = composition_table{:,4:95};

for s=1:numel(SIDs)
    summary_table.num_lineages(s) = sum(comp_matrix(s,:)>0,2);
    covered_isolates = sum(comp_matrix(s,:).*in_collection);
    total_isolates = sum(comp_matrix(s,:));
    summary_table.collection_coverage(s) = covered_isolates./total_isolates;
end

% Identify samples with low collection coverage
collection_coverage_cutoff = 0.75;
collection_coverage_SIDs = {summary_table.SID{summary_table.collection_coverage>collection_coverage_cutoff}}';

for n=1:numel(sepi_isolates.Family)
    sepi_isolates.in_collection(n) = any(strcmp(collection_coverage_SIDs,sepi_isolates.SID(n)));
end


%% Final filtering of composition and isolate tables

sepi_isolates = sepi_isolates(sepi_isolates.in_collection,:);
composition_table = composition_table(summary_table.collection_coverage>collection_coverage_cutoff,:);
comp_matrix = composition_table{:,4:95};
summary_table = summary_table(summary_table.collection_coverage>collection_coverage_cutoff,:);

[G,SIDs] = findgroups(sepi_isolates.SID);
num_SID = splitapply(@numel,sepi_isolates.SampleName,G);

[G,subs] = findgroups(sepi_isolates.Subject);
num_Subject = splitapply(@numel,sepi_isolates.SampleName,G);

[G,fams] = findgroups(sepi_isolates.Family);
num_Family = splitapply(@numel,sepi_isolates.SampleName,G);

summary_table = table(collection_coverage_SIDs);
num=numel(collection_coverage_SIDs);
for n=1:num
    summary_table.Family{n} = char(summary_table.collection_coverage_SIDs{n}(1));
    summary_table.num_Family(n) = num_Family(find(strcmp(fams,summary_table.Family{n})));
    summary_table.Subject{n} = char(summary_table.collection_coverage_SIDs{n}(1:3));
    summary_table.num_Subject(n) = num_Subject(find(strcmp(subs,summary_table.Subject{n})));
    summary_table.SID{n} = char(summary_table.collection_coverage_SIDs{n}(1:4));
    summary_table.num_SID(n) = num_SID(find(strcmp(SIDs,summary_table.SID{n})));
    summary_table.num_lineages(n) = sum(comp_matrix(n,:)>0,2);
    covered_isolates = sum(comp_matrix(n,:).*in_collection);
    total_isolates = sum(comp_matrix(n,:));
    summary_table.collection_coverage(n) = covered_isolates./total_isolates;
end
summary_table = removevars(summary_table, 'collection_coverage_SIDs');

%% Make summary graphs
if plot_summary_option
    figure(9)
    subplot(3,1,1)
    bar(num_Family)
    xticks(1:numel(unique(sepi_isolates.Family)))
    xticklabels(unique(sepi_isolates.Family))
    
    
    subplot(3,1,2)
    bar(num_Subject)
    xticks(1:numel(unique(sepi_isolates.Subject)))
    xticklabels(unique(sepi_isolates.Subject))
    
    subplot(3,1,3)
    bar(num_SID)
    xticks(1:numel(unique(sepi_isolates.SID)))
    xticklabels(unique(sepi_isolates.SID))
    
    
    figure(10)
    scatter(summary_table.num_SID,summary_table.num_Subject)
    ylabel('Number of Isolates in Subject')
    xlabel('Number of Isolates in SID')
end

%% Create final tables
clearvars -except headerdata summary_table sepi_isolates initial_composition_table in_collection MG_filename
collection_representation_cutoff = 0.75;

% Filter out samples with low coverage
is_represented = summary_table.collection_coverage>collection_representation_cutoff;
good_SIDs = {summary_table.SID{is_represented}}';

for s=1:numel(good_SIDs)
    sid = good_SIDs{s};    
    idxs = find(strcmp(sepi_isolates.SID,sid));
    sepi_isolates.well_represented(idxs) = true;
end

final_summary_table = summary_table(is_represented,:);
final_sepi_isolates = sepi_isolates(sepi_isolates.diverse_lineage & sepi_isolates.good_coverage & sepi_isolates.well_represented,:);
final_composition_table = table();
for s = 1:numel(good_SIDs)
    sid = good_SIDs{s};    
    idxs = find(strcmp(final_sepi_isolates.SID,sid)&final_sepi_isolates.diverse_lineage&final_sepi_isolates.good_coverage);
    final_composition_table.SID{s} = sid;
    final_composition_table.Family{s} = str2num(char(sid(1)));
    final_composition_table.Subject{s} = char(sid(1:3));
    for n = 1:max(sepi_isolates.Lineage)
        col_name = ['Lineage_' num2str(n)];
        final_composition_table.(col_name)(s) = sum(final_sepi_isolates.Lineage(idxs) == n);
    end
end


%% Save composition
composition_table = final_composition_table;
for s=1:numel(good_SIDs)
    comp = composition_table{s,4:end};
    norm_comp = comp./sum(comp);
    composition_table{s,4:end} = norm_comp;
end

% uncomment to regenerate composition_table.mat
%save('composition_table.mat','composition_table')

%% Read metagenomic composition
MG_data = readtable(MG_filename,'Delimiter',',');
MG_data(:,2:end) = fillmissing(MG_data(:,2:end),'constant',0);
numS = size(MG_data,1);

% remove non-included timepoints

MG_idxs = [];
for n=1:numS
    MG_data.Family{n} = str2num(char(MG_data.SID{n}(1)));
    MG_data.Subject{n} = char(MG_data.SID{n}(1:3));
end
numS = size(MG_data,1);
numC = size(MG_data,2);

MG_data = [MG_data(:,1) MG_data(:,(numC-1):(numC)) MG_data(:,2:(numC-2))];
for n=1:numS
    % Include metagenomics if >0.7 of community iscalled
    if ismember(MG_data.SID{n},good_SIDs) && (sum(MG_data{n,4:end})>0.7)
        MG_idxs = [MG_idxs, n];
    end
end
MG_table = MG_data(MG_idxs,:);
MG_table.SID = strcat({MG_table.SID{:}}', '(MG)');

%% Composition plots per family with metagenomes
for f=1:8
    fam_idxs_i = ([composition_table.Family{:}]==f);
    fam_idxs_m = ([MG_table.Family{:}]==f);

    if any(fam_idxs_i)
        fighandle = figure(f);
        clf(f)
        fighandle.Position=[fighandle.Position(1) fighandle.Position(2) 400 375];
        % Make subtable with both
        family_table = [composition_table(fam_idxs_i,:); MG_table(fam_idxs_m,:)];
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
        color_vals(~make_grey,:) = turbo(sum(~make_grey));
        if sum(~make_grey)>10
            color_vals(1:2:end,:)=color_vals(1:2:end,:).*0.7;
        end
        % color_vals(make_grey,:) = [mean(color_vals(make_grey,:),2) mean(color_vals(make_grey,:),2) mean(color_vals(make_grey,:),2)];
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

