import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import scipy.stats as stats


df = pd.read_csv('pangenome_output_table.csv')
num_lineages = max(df.loc[:, 'No. isolates'])
core_thresh = 0.95*num_lineages
shell_thresh = 0.15*num_lineages


#%% calculate enrichment of gene group
total_core = sum(df.loc[:, 'No. isolates']>core_thresh)
total_shell = sum(df.loc[:, 'No. isolates']>shell_thresh)-total_core
total_cloud = sum(df.loc[:, 'No. isolates']>0)-total_core-total_shell
group = ['Total']
core_num = [total_core]
shell_num = [total_shell]
cloud_num = [total_cloud]
odds = [0]
p = [0]

columns = ['amrfinder_hit','virulence_hit','defensefinder_hit']
for col in columns:
    non_nan_indices = np.where(df[col].notna())
    core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
    pseudo=False
    if core==0:
        pseudo = True
        core=1 #add pseudo count
    shell = sum(df.loc[non_nan_indices, 'No. isolates']>shell_thresh)-core
    cloud = sum(df.loc[non_nan_indices, 'No. isolates']>0)-core-shell
    contingency_table = [[core,total_core-core],[shell,total_shell-shell],[cloud,total_cloud-cloud]]
    odds_ratio, p_value = stats.fisher_exact([[cloud, core],[total_cloud-cloud,total_core-core]])
    group.append(col)
    if pseudo:
        core=0 
    core_num.append(core)
    shell_num.append(shell)
    cloud_num.append(cloud)
    odds.append(odds_ratio)
    p.append(p_value)

COG_dictionary = pd.read_csv('COG_dictionary.csv')
for code in COG_dictionary['Code']:
    non_nan_indices = np.where(df['COG_category'].fillna('').str.contains(code))
    core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
    shell = sum(df.loc[non_nan_indices, 'No. isolates']>shell_thresh)-core
    cloud = sum(df.loc[non_nan_indices, 'No. isolates']>0)-core-shell
    odds_ratio, p_value = stats.fisher_exact([[cloud, core],[total_cloud-cloud,total_core-core]])
    group.append(code)
    core_num.append(core)
    shell_num.append(shell)
    cloud_num.append(cloud)
    odds.append(odds_ratio)
    p.append(p_value)
    
df_stats = pd.DataFrame({'Group':group,'Core Number':core_num,'Shell Number':shell_num,'Cloud Number':cloud_num,'Odds Ratio':odds,'p_value':p})
df_stats.to_csv('stats_output.csv')


#%% Calculate enrichment in recent gains/losses
total_gainloss = sum(df.loc[:, 'num_gainlosses']>0)
group = ['Total']
core_num = [total_core]
gainloss_num = [total_gainloss]
odds = [0]
p = [0]

for col in columns:
    non_nan_indices = np.where(df[col].notna())
    core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
    pseudo=False
    if core==0:
        pseudo = True
        core=1 #add pseudo count
    gainloss = sum((df.loc[non_nan_indices, 'num_gainlosses']>0))
    odds_ratio, p_value = stats.fisher_exact([[gainloss, core],[total_gainloss-gainloss,total_core-core]])
    group.append(col)
    if pseudo:
        core=0        
    core_num.append(core)
    gainloss_num.append(gainloss)
    odds.append(odds_ratio)
    p.append(p_value)

for code in COG_dictionary['Code']:
    non_nan_indices = np.where(df['COG_category'].fillna('').str.contains(code))
    core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
    gainloss = sum((df.loc[non_nan_indices, 'num_gainlosses']>0))
    odds_ratio, p_value = stats.fisher_exact([[gainloss, core],[total_gainloss-gainloss,total_core-core]])
    group.append(code)
    core_num.append(core)
    gainloss_num.append(gainloss)
    odds.append(odds_ratio)
    p.append(p_value)
    
    
df_stats = pd.DataFrame({'Group':group,'Core Number':core_num,'Gainloss Number':gainloss_num,'Odds Ratio':odds,'p_value':p})
df_stats.to_csv('gainloss_stats_output.csv')

#%% Calculate enrichment in recent gains/losses with duplication
total_gainloss = sum(df['num_gainlosses'])
group = ['Total']
core_num = [total_core]
gainloss_num = [total_gainloss]
odds = [0]
p = [0]

for col in columns:
    non_nan_indices = np.where(df[col].notna())
    core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
    pseudo=False
    if core==0:
        pseudo = True
        core=1 #add pseudo count
    gainloss = sum(df.loc[non_nan_indices, 'num_gainlosses'])
    odds_ratio, p_value = stats.fisher_exact([[gainloss, core],[total_gainloss-gainloss,total_core-core]])
    group.append(col)
    if pseudo:
        core=0  
    core_num.append(core)
    gainloss_num.append(gainloss)
    odds.append(odds_ratio)
    p.append(p_value)

for code in COG_dictionary['Code']:
    non_nan_indices = np.where(df['COG_category'].fillna('').str.contains(code))
    core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
    gainloss = sum(df.loc[non_nan_indices, 'num_gainlosses'])
    odds_ratio, p_value = stats.fisher_exact([[gainloss, core],[total_gainloss-gainloss,total_core-core]])
    group.append(code)
    core_num.append(core)
    gainloss_num.append(gainloss)
    odds.append(odds_ratio)
    p.append(p_value)
    
    
df_stats = pd.DataFrame({'Group':group,'Core Number':core_num,'Gainloss Number':gainloss_num,'Odds Ratio':odds,'p_value':p})
df_stats.to_csv('gainloss_dupstats_output.csv')


#%% Calculate enrichment in non-core
total_noncore = total_shell+total_cloud
group = ['Total']
core_num = [total_core]
noncore_num = [total_noncore]
odds = [0]
p = [0]

for col in columns:
    non_nan_indices = np.where(df[col].notna())
    core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
    pseudo=False
    if core==0:
        pseudo = True
        core=1 #add pseudo count
    noncore = sum(((df.loc[non_nan_indices, 'No. isolates']>0)&(df.loc[non_nan_indices, 'No. isolates']<=core_thresh))|(df.loc[non_nan_indices, 'num_gainlosses']>0))
    odds_ratio, p_value = stats.fisher_exact([[noncore, core],[total_noncore-noncore,total_core-core]])
    group.append(col)
    if pseudo:
        core=0  
    core_num.append(core)
    noncore_num.append(noncore)
    odds.append(odds_ratio)
    p.append(p_value)

for code in COG_dictionary['Code']:
    non_nan_indices = np.where(df['COG_category'].fillna('').str.contains(code))
    core = sum(df.loc[non_nan_indices, 'No. isolates']>core_thresh)
    noncore = sum(((df.loc[non_nan_indices, 'No. isolates']>0)&(df.loc[non_nan_indices, 'No. isolates']<=core_thresh))|(df.loc[non_nan_indices, 'num_gainlosses']>0))
    odds_ratio, p_value = stats.fisher_exact([[noncore, core],[total_noncore-noncore,total_core-core]])
    group.append(code)
    core_num.append(core)
    noncore_num.append(noncore)
    odds.append(odds_ratio)
    p.append(p_value)
    
    
df_stats = pd.DataFrame({'Group':group,'Core Number':core_num,'Noncore Number':noncore_num,'Odds Ratio':odds,'p_value':p})
df_stats.to_csv('combined_stats_output.csv')

#%% Calculate enrichment in BGCs
total_core_BGC = 5
total_noncore_BGC = 7 #pseudocount
noncore_BGC = 6
core_BGC = 1 #pseudocount

odds_ratio, p_value = stats.fisher_exact([[noncore_BGC, core_BGC],[total_noncore_BGC-noncore_BGC,total_core_BGC-core_BGC]])

