function fighandle = plot_intralineage_variation(interactions_by_isolate,variation_threshold,list_option,fignum)
%list_option list out all potential cases of intralineage varaition? True/False

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 400 400];

metadata_table = struct2table(interactions_by_isolate.metadata);

%% Group samples by lineage
interactions_by_group = struct();
[group, id] = findgroups(metadata_table.Lineage);

% check for intralineage variation
num_variant_ant = 0;
num_invariant_ant = 0;
num_invariant_nonant = 0;
num_variant_sen = 0;
num_invariant_sen = 0;
num_invariant_nonsen = 0;

variant_ant_table = table(0,0,{0},{0},'VariableNames',{'Antagonist','Bait','AUCs','Calls'});
variant_sen_table = table(0,0,{0},{0},'VariableNames',{'Antagonist','Bait','AUCs','Calls'});


for a =1:numel(id) %antagonist
    for b = 1:numel(id) %bait
        aidxs = find(group==a);
        bidxs = find(group==b);
        
        % variation in antagonism
        if numel(aidxs)>1
            for i=1:numel(bidxs)
                ZOI_call_vector = [interactions_by_isolate.ZOI_call(bidxs(i),aidxs)];
                ZOI_call_vector = ZOI_call_vector(:);
                ZOI_AUCs_vector = [interactions_by_isolate.ZOI_AUC(bidxs(i),aidxs)];
                ZOI_AUCs_vector = ZOI_AUCs_vector(:);

                if any(ZOI_call_vector) && ~all(ZOI_call_vector)
                    if max(ZOI_AUCs_vector)-min(ZOI_AUCs_vector) >variation_threshold
                        if list_option
                            disp([num2str(id(a)) ' varies antagonizing ' num2str(id(b))])
                        end
                    end
                    num_variant_ant = num_variant_ant+1;
                    variant_ant_table = [variant_ant_table; table(id(a),id(b),{ZOI_AUCs_vector},{ZOI_call_vector},'VariableNames',{'Antagonist','Bait','AUCs','Calls'})];
                elseif all(ZOI_call_vector)
                    num_invariant_ant = num_invariant_ant+1;
                elseif ~any(ZOI_call_vector)
                    num_invariant_nonant = num_invariant_nonant+1;
                end
            end
        end

        % variation in sensitivity
        if numel(bidxs)>1
            for i=1:numel(aidxs)
                ZOI_call_vector = [interactions_by_isolate.ZOI_call(bidxs,aidxs(i))];
                ZOI_call_vector = ZOI_call_vector(:);
                ZOI_AUCs_vector = [interactions_by_isolate.ZOI_AUC(bidxs,aidxs(i))];
                ZOI_AUCs_vector = ZOI_AUCs_vector(:);

                if any(ZOI_call_vector) && ~all(ZOI_call_vector)
                    if max(ZOI_AUCs_vector)-min(ZOI_AUCs_vector) >variation_threshold
                        if list_option
                            disp([num2str(id(a)) ' causes variable sensitivity in ' num2str(id(b))])
                        end
                    end
                    num_variant_sen = num_variant_sen+1;
                    variant_sen_table = [variant_sen_table; table(id(a),id(b),{ZOI_AUCs_vector},{ZOI_call_vector},'VariableNames',{'Antagonist','Bait','AUCs','Calls'})];
                elseif all(ZOI_call_vector)
                    num_invariant_sen = num_invariant_sen+1;
                elseif ~any(ZOI_call_vector)
                    num_invariant_nonsen = num_invariant_nonsen+1;
                end
            end
        end
    end
end
variant_ant_table(1,:) = [];
variant_sen_table(1,:) = [];

fraction_variant_ant = num_variant_ant./(num_variant_ant+num_invariant_ant+num_invariant_nonant);
fraction_variant_sen = num_variant_sen./(num_variant_sen+num_invariant_sen+num_invariant_nonsen);
disp([num2str(fraction_variant_ant,3) ' of antagonisms vary across lineage'])
disp([num2str(fraction_variant_sen,3) ' of sensitivities vary across lineage'])

%%
figure(fignum)
ant_range = zeros(num_variant_ant,1);
sen_range = zeros(num_variant_sen,1);
for n=1:num_variant_ant
    data = variant_ant_table{n,3}{:};
    ant_range(n) = max(data)-min(data);
end
for n=1:num_variant_sen
    data = variant_sen_table{n,3}{:};
    sen_range(n) = max(data)-min(data);
end

subplot(2,2,1)
plot(1:num_variant_ant,ant_range,'.')
ylim([0 max(max(interactions_by_isolate.ZOI_AUC))])
title('Intralineage Variation in Antagonism')
hold on
plot([0 num_variant_ant],[variation_threshold, variation_threshold],'r')

subplot(2,2,2)
histogram(ant_range,[0:100:3000])
title('Intralineage Variation in Antagonism')
hold on
plot([variation_threshold variation_threshold],[0 100],'r')
ylim([0 100])


subplot(2,2,3)
plot(1:num_variant_sen,sen_range,'.')
ylim([0 max(max(interactions_by_isolate.ZOI_AUC))])
title('Intralineage Variation in Sensitivity')
hold on
plot([0 num_variant_sen],[variation_threshold, variation_threshold],'r')

subplot(2,2,4)
histogram(sen_range,[0:100:3000])
title('Intralineage Variation in Sensitivity')
hold on
plot([variation_threshold variation_threshold],[0 100],'r')
ylim([0 100])

