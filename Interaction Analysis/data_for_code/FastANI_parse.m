clc
clear all
close all

datatable = readtable('FastANI_output.csv');
varnames = datatable.Properties.VariableNames;

data = [nan(1,159); datatable{:,:}];
idxs=[];
for n=1:159
    data(n,n)=100;
    for m=n+1:159
        data(n,m)=data(m,n);
    end
    if mode(data(n,:)==0)
        idxs = [idxs n];
    end
    varnames{n}=varnames{n}(2:end);
end
data(idxs,:)=[];
data(:,idxs)=[];
ANI=data./100;

missing_pair = sum(data==0,'all');
varnames(:,idxs)=[];
varnames=varnames';

% add distance pairs to header data
full_header_data = readtable('all_header_data.csv');

for n=1:numel(full_header_data.Name)
    idx = find(strcmp(varnames,full_header_data.Name{n}));
    if idx
        full_header_data.Distance_Index(n) = idx;
    else
        full_header_data.Distance_Index(n) = 0;
    end
end