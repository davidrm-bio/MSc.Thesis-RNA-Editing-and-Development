'''
@author: David
@date: May 2023
@Description: Processing of REDItools2 output

'''
# Modules
import pandas as pd
from MyFunctions import *

# Part 0 - Set-up 
## Load Metadata and separate the sample ID for the 2 conditions
metadata = pd.read_csv('Mice_KO/MetaData_PRJNA593199.txt', sep=',')
column_group = 'genotype/variation'
metadata_condition =  metadata.groupby(column_group)
conditions = [metadata_condition.get_group(x) for x in metadata_condition.groups]

condD_id = list(conditions[0]['Run']) # mice KO
condD_name = list(conditions[0][column_group])[0]
condC_id = list(conditions[1]['Run']) # mice WT
condC_name = list(conditions[1][column_group])[0]


# Part 1 - Processing of Samples
print ('##### - Processing of Samples in KO - #####')
condD_stats = process_output (condD_id, path='Mice_KO/', depth=5, freq=0.1)
print ('\n##### - Processing of Samples in WT - #####')
condC_stats = process_output (condC_id, path='Mice_KO/', depth=5, freq=0.1)


# Part 2 - Statistical analysis
'''
Binomial test. Given a RES, 
1) Present in at least 1 Sample (WT) AND in at least 1 Sample (KO) >> Binomial test
2) Present in 2 Sample (WT) AND 0 sample in KO >> Significant
3) Present in 2 Sample (KO) AND 0 sample in WT >> Significant
4) The Rest (1/0 WT and 1/0 KO) >> Not significant
We only keep the RNA editing sites that are significant (Pval < 0.05 and cases 2 and 3 from above)
'''
RES_mice = statistical_analysis(condD_stats, condD_id, condC_stats, condC_id,1,1,0, 2)
RES_mice_clean, SubMice_clean = clean_results(RES_mice)
RES_mice_sig = significant_plot_all(RES_mice_clean, condA_name = 'mice KO condition', condB_name= 'mice WT condition', filename='miceall.png')

Mice_sig, Mice_WTSig, Mice_KOSig = get_sig_EP(RES_mice_clean)

###############################################

# Part 1 - PCA
df_mice = principal_component_analysis ('Mice_KO/Bedtools_count/',  D_id = condD_id, C_id = condC_id, filename='PCAMouse.png')

# Part 2 - Cleaning & Processing
df_mice_clean = df_mice[df_mice['color'] != 'grey']
df_mice_clean['Pval'] = ttest_ind(df_mice_clean[condD_id], df_mice_clean[condC_id], axis=1)[1] # Compute Pval

tomato = df_mice_clean[(df_mice_clean['Pval'] < 0.05) & (df_mice_clean['control'] > df_mice_clean['disease'])]
gold = df_mice_clean[(df_mice_clean['Pval'] < 0.05) & (df_mice_clean['control'] < df_mice_clean['disease'])]

# Part 3 - Save
gene_tomato = [gene for gene in tomato.GeneName if 'Gm' not in gene]
gene_gold = [gene for gene in gold.GeneName if 'Gm' not in gene]

pd.DataFrame(gene_tomato).to_csv('RawCounts-WT>KO.txt', index=False, header=False)
pd.DataFrame(gene_gold).to_csv('RawCounts-WT<KO.txt', index=False, header=False)

