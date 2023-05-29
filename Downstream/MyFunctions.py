#### Set-up: Modules and Functions

### Modules
import pandas as pd
import numpy as np
from scipy import stats as sp
import matplotlib.pyplot as plt
import ujson as json
import statsmodels.api as sm
import os
import time
import subprocess
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

############################################
# Modules to use R functions within python
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

# Needed for statistical_analysis to use the
# glm R function in python
stats = importr('stats')
glm = stats.glm
############################################

import warnings
warnings.filterwarnings('ignore')



### Functions included in the analysis

def process_output(cond, path='./', suffix = '_noSNP.out', depth=5, freq=0.1):
    '''
    Downstream processing for REDItools2 output. 
    Input: 
        - cond: list of SRR of the samples. Only the basename. Example: ['SRR10581813', 'SRR10581814']
        - path: path where the files are located. 
        - suffix: suffix that follows the SRR of the samples specify in cond. 
        - depth: Minimum coverage to keep a site.
        - freq: Minimum frequency / RNA editing proportion to keep a site.
    
    Only RNA editing sites  (A > G & T > C) are kept. These sites have a coverage >= 5 (default)
    and frequency >= 1 (default).
    
    Output example:
        - Dictionary with the following format:
            {'chr1|4565256|AG': # RES ID: chr | position | RES type
                {
                    'Samples': ['SRR10581813', 'SRR10581814'], # Sample ID
                    'Strand': ['2', '-'],  # Strand Information
                    'Coverage': [7, 0],  # Site total coverage
                    'EP': [0.14, 0],  # Site editing proportion
                    'A_count': [6, 0],  # Number of A's
                    'C_count': [0, 0],  # Number of C's
                    'G_count': [1, 0],  # Number of G's
                    'T_count': [0, 0]   # Number of T's
                },
            }
    
    '''
    
    # Local Variables
    entries = ['Samples', # Sample ID
               'Strand',  # Strand: 2 (unknown), 0 (- strand), 1 (+ strand)
               'Coverage', 
               'EP', # Editing proportion == frequency
               'A_count', 
               'C_count', 
               'G_count', 
               'T_count']
    
    res = ['AG', 'TC']
    all_res = {}
    
    s0 = time.time() # Record time
    
    sampl_cont = 0 # count number of samples 
    for basename in cond:
        
        with open(os.path.join(path, basename + suffix)) as Read:
            for lines in Read:
                line = lines.rstrip().split('\t')
                
                # Positions to keep:
                # Have AG or TC 
                # Minimum Coverage of 5 (default)
                # Minimum Frequency of 0.1 (default)
                
                if line[7] in res and int (line[4]) >= depth and float(line[8]) >= freq:
                    RES_id = line[0] + '|' + line[1] + '|' + line[7] # ID example: Chr1|1763826|AG
                    
                    if RES_id not in all_res:
                        all_res.update({RES_id:{entry:[] for entry in entries}})
                
                    # Save data
                    all_res[RES_id]['Samples'].append(basename) # Sample Name
                    all_res[RES_id]['Strand'].append(line[3]) # Strand 
                    all_res[RES_id]['Coverage'].append(int (line[4])) # Coverage 
                    all_res[RES_id]['EP'].append(float (line[8])) # Frequency 
                    
                    # Needed for significant test
                    base_counts = line[6].replace('[', '').replace(']', '').split(',') # [A,C,G,T] -> list(A, C, G, T)
                    all_res[RES_id]['A_count'].append(int (base_counts[0])) # A count
                    all_res[RES_id]['C_count'].append(int (base_counts[1])) # C count
                    all_res[RES_id]['G_count'].append(int (base_counts[2])) # G count
                    all_res[RES_id]['T_count'].append(int (base_counts[3])) # T count
                    
        sampl_cont += 1
        
    ## Add samples that do not have the RES and set values to 0
    
    for RES_id in all_res:
        if len (all_res[RES_id]['Samples']) < sampl_cont:
            s_not = [s for s in cond if s not in all_res[RES_id]['Samples']] # Add sample name 
            
            all_res[RES_id]['Samples'] = all_res[RES_id]['Samples'] + s_not
            all_res[RES_id]['Strand'] = all_res[RES_id]['Strand'] + ['-']*len(s_not)
            all_res[RES_id]['Coverage'] = all_res[RES_id]['Coverage'] + [0]*len(s_not)
            all_res[RES_id]['EP'] = all_res[RES_id]['EP'] + [0]*len(s_not)
            
            all_res[RES_id]['A_count'] = all_res[RES_id]['A_count'] + [0]*len(s_not)
            all_res[RES_id]['C_count'] = all_res[RES_id]['C_count'] + [0]*len(s_not)
            all_res[RES_id]['G_count'] = all_res[RES_id]['G_count'] + [0]*len(s_not)
            all_res[RES_id]['T_count'] = all_res[RES_id]['T_count'] + [0]*len(s_not)
    
    check = len (all_res[RES_id]['A_count'])
    assert check == sampl_cont, f'The number of records per site is different to the number of samples {check} and {sampl_cont}'
    
    s1 = time.time()
    
    ## Check-point
    print ('Time needed to process {} samples, {} s'.format(len(cond), round (s1-s0, 5)))
    
    return all_res


def statistical_analysis (condA, condA_id, condB, condB_id, cutoffA, cutoffB, altA, altB):
    '''
    Statistical analysis for each RNA editing position. 
    Input:
        - condA / condB: dictionary with the format generated by the function process_output
        - cutoffA / cutoffb: minimum number of samples that need to have RES present to apply 
        the binomial test
    
    For each RNA editing site a binomial test is perform that compares the number of As/Ts and Gs/Cs 
    between condA and condB. For instance, consider the position chr1|58885457|AG where we have the
    following:
        df = {edit: [1,0,5,1], ref:[4,0,6,0], condition: [KO,KO,WT,WT]}
    
    In condA (KO) we have 1 sample with G and A counts. In the condB (WT) we have 2 samples with G and A
    counts. Using the r function (We use a module that allow the use of R functions from python) GLM 
    we apply the binomial test and return the Pvalue. 
    
    In case of having only counts in 1 condition we set the pvalue to np.nan. The binomial test is apply
    only if both condA and condB has counts from at least 1 sample.
    
    NOTE CONDB SHOULD ALWAYS BE WT SAMPLES WHILE CONDA THE DISEASE STATE OR KO CONDITION
    
    Output: 
        - dict_data: dictionary where key == RES_ids and values == Pvalue

    '''
    s0 = time.time() # Check time of the analysis
    
    # Get Total Number of Unique RES in Both conditions
    total_res = set(condA.keys())
    total_res.update (list (condB.keys()))
    
    # Local Variables
    cont = 0
            
    entries = ['mean', 'std']
    conditions = ['condA', 'condB']
    big_data = {res: {cond:
                      {entry:np.nan for entry in entries} 
                      for cond in conditions} 
                for res in total_res}
    
    condA_s = len(condA_id) # Get number of samples in condA
    condB_s = len (condB_id) # Get number of samples in condB
    total_samples = condA_s + condB_s 
    condition = ['KO']*condA_s + ['WT']*condB_s
    
    print ('\nThere are {} RES in total. Proceeding with statistical analysis...'.format(len(total_res)))
    
    for res in total_res:
        condA_present, condB_present = 0, 0 # To check if RES is present in both conditions
        A_pass, B_pass= 0,0
        save_A, save_B = 0,0
        save_Av2, save_Bv2 = 0,0
        
        ### Part 1 - ref and edt are a lists of counts for the ref (A/T) and edited (G/C) base counts
        
        # Define an empty list of counts 
        ref = [0] * total_samples
        edt = [0] * total_samples
        
        # Update with condA counts 
        if res in condA:
            
            # Save mean and standard deviation of res
            big_data[res]['condA']['mean']= np.mean(np.array(condA[res]['EP'])[np.array(condA[res]['EP']) !=0])
            big_data[res]['condA']['std']= np.std(np.array(condA[res]['EP']) [np.array(condA[res]['EP']) !=0])
            
            condA_present = 1
            if 'AG' in res:
                ref[:condA_s] = condA[res]['A_count']
                edt[:condA_s] = condA[res]['G_count']
            else:
                ref[:condA_s] = condA[res]['T_count']
                edt[:condA_s] = condA[res]['C_count']
                
            # Extra check to decide on doing binomial test
            checkA = condA[res]['Strand']
            if len(checkA) - checkA.count('-') >= cutoffA:
                A_pass = 1
            
            # Save RES - 
            # For 2 samples / conditon
            ## Present in 2 in WT and None in KO
            # For 4 samples / condition
            ## Present in >=3 in WT and <=1 in KO
            if len(checkA) - checkA.count('-') <= altA:
                save_A = 1
            
            if len(checkA) - checkA.count('-') >= altB:
                save_Av2 = 1
        else:
            save_A = 1
            
        # Update with condB counts
        if res in condB:
            
            # Save mean and standard deviation of res
            big_data[res]['condB']['mean']= np.mean(np.array(condB[res]['EP'])[np.array(condB[res]['EP']) != 0])
            big_data[res]['condB']['std']= np.std(np.array(condB[res]['EP'])[np.array(condB[res]['EP']) != 0])
            
            condB_present = 1
            if 'AG' in res:
                ref[condA_s:] = condB[res]['A_count']
                edt[condA_s:] = condB[res]['G_count']
            else:
                ref[condA_s:] = condB[res]['T_count']
                edt[condA_s:] = condB[res]['C_count']
                
            
            # Extra check to decide on doing binomial test
            checkB = condB[res]['Strand']
            if len(checkB) - checkB.count('-') >= cutoffB:
                B_pass = 1
            
            if len(checkB) - checkB.count('-') >= altB:
                save_B = 1
            
            if len(checkB) - checkB.count('-') <= altA:
                save_Bv2 = 1
        else:
            save_Bv2 = 1
        
        df = pd.DataFrame ({'edt':edt, 'ref':ref, 'condition':condition})
    
        ### Part 2 - Statistical analysis
        if condB_present + condA_present == 2 and A_pass + B_pass == 2:
            # Using R functions
            with (ro.default_converter + pandas2ri.converter).context():
                df2 = ro.conversion.get_conversion().py2rpy(df)
            my_glm = glm('cbind(edt, ref) ~ condition', data=df2, family='binomial')
            pval = float (ro.r.summary(my_glm).rx2("coefficients")[7])
        
        elif save_A + save_B == 2:
            pval = -1
        
        elif save_Av2 + save_Bv2 == 2:
            pval = -2

        else:
            pval = np.nan # RES is only present in 1 condition or not sufficient samples for the test

        big_data[res]['Pval']=  pval

        cont +=1 
        
        if cont % 50000 ==0:
            print ('{} So far...'.format(cont))

    s1 = time.time ()
    
    print ('Time needed to process {} samples, {} s'.format(total_samples, round (s1-s0, 5)))
    
    return big_data


def clean_results (data):
    '''
    Remove RNA editing sites with Pval == NaN.
    - Input:
        data: dictionary generated by statistical_analysis
    
    This function remove RNA editing sites whose Pval is NaN. The Pval depends on the number of sample
    a RNA editing site is present
    
    + For 2 sample per condition:
        - 1 sample of WT and 0 sample of KO >> NaN
        - 2 sample of WT and 0 sample of KO >> -1
        - 1 or 2 sample of WT and 1 or 2 sample of KO >> Pval (Binomial test)
        
    + For 4 sample per condition: (cutoff 2, 2; Saving 0/1, 3/4)
        - 0 or 1 of WT and 0 or 1 of KO >> NaN
        - 3 or 4 of WT and 0 or 1 of KO >> -1
        - 2, 3 or 4 of WT and 2, 3 or 4 of KO >> Pval (Binomial test)
    
    In this fraction only cases with a Pval != NaN are kept
    '''
    # Local variables
    data_clean = {}
    subdata_clean = {}

    total = 0
    subtotal = 0
    
    print ('There are {} RES in total\n'.format(len(data)))
    
    for res in data:
        if str(data[res]['Pval']) != 'nan':
            data_clean[res] = data[res]
            total +=1

            condA = data[res]['condA']['mean'] 
            condB = data[res]['condB']['mean']

            if data[res]['Pval'] not in [-1, -2]:
                subdata_clean[res] = data[res]

            if data[res]['Pval'] == -1:
                subtotal +=1
                if str(condA) == 'nan':
                    data_clean[res]['condA']['mean'] = 0
            if data[res]['Pval'] == -2:
                subtotal +=1
                if str(condB) == 'nan':
                    data_clean[res]['condB']['mean'] = 0
    
    print ('Cleaning finished. There are {} RES with significant difference ({} cases were saved)'.format(total, subtotal))
    return data_clean, subdata_clean
    

def significant_plot_all (data_clean, condA_name, condB_name, filename):
    '''
    Scatterplot of Editing Propotion (EP) of 2 conditions (condA Vs conB). Significant values are colored red.
    
    Input:
        - data_clean: dictionary generated by the function statistical analysis. 
    The input data can be pre-process to remove positions with a pvalue of nan. These positions correspond to RNA
    editing sites only present in 1 condition.
    
    NO FDR CORRECTION IS APPLIED FOR THIS FUNCTION. 
    '''
    
    # Local Variables
    A_sig, B_sig = [], [] 
    A_notsig, B_notsig = [], []
    
    RES_sig = []

    # Pre-plotting - Separate Significant and non-significant RES for each condition
    for res in data_clean:
        if data_clean[res]['Pval'] <= 0.05:
            A_sig.append(data_clean[res]['condA']['mean'])
            B_sig.append(data_clean[res]['condB']['mean'])
            RES_sig.append(res)
        else:
            A_notsig.append(data_clean[res]['condA']['mean'])
            B_notsig.append(data_clean[res]['condB']['mean'])

    # Plot
    plt.figure(figsize=(20, 13), dpi=300)
    plt.title('Scatterplot of the average Editing Proportion for\nreported RNA editing sites', fontsize=30)
    plt.ylabel(condA_name, fontsize = 23)
    plt.xlabel(condB_name, fontsize = 23)
    
    plt.scatter(B_sig, A_sig, color='tomato', label='<0.05')
    plt.scatter(B_notsig, A_notsig, color='grey', alpha=0.5, label='>0.05')
    plt.axline((0.1,0.1), slope=1, color='black', lw=2)
    
    plt.legend(prop={'size': 20}, loc=4, title= 'Significantly edited', title_fontsize='large')
    plt.savefig(filename)
    return  RES_sig


def principal_component_analysis (path, D_id, C_id, filename, suffix = '_clean_RawCounts.bed', greaterX=15, greaterY=10):
    '''
    PCA analysis.
    - Input:
        path: path to GTF files with the counts of RES per gene
        D/C_id: the sample ID for disease and control condition

    For coloring the PCA, the mean value of raw nimber of RES per condition was calculated.
    '''
    # Part 0 - Set up | Local variables
    files = [f for f in os.listdir(path) if 'RawCounts' in f]
    exist = 0

    # Part 1 - Analysis

    ## Part 1.1 - Create df for PCA
    for f in files:
        basename = f.replace('_clean_RawCounts.bed', '') # Sample Name
        # Structure data
        if exist == 0:
            df = pd.read_csv(os.path.join(path, f), sep='\t', header=None) # Load GTF with Raw counts
            df['GeneName'] = df[8].str.extract(r'gene_name ([.?(?)?\w-]+);') # Extract Gene_name
            df.rename(columns={9: basename}, inplace=True)
            df = df[[0,3,4,6,'GeneName',basename]]
            exist=1
        else:
            newdf = pd.read_csv(os.path.join(path, f), sep='\t', header=None) 
            newdf.rename(columns={9: basename}, inplace=True)
            df[basename] = newdf[basename]
    df.rename(columns={0: 'chr', 3: 'start', 4:'end', 6:'strand',}, inplace=True)  

    ## Part 1.2 - PCA
    filt = D_id + C_id
    df_clean = df[filt]
    # Normalize data
    scaling=StandardScaler()
    df_norm = scaling.fit_transform(df_clean)
    df_norm = pd.DataFrame(df_norm)
    # PCA
    pca = PCA(n_components=df_norm.shape[1])
    pca.fit(df_norm)
    x = pca.transform(df_norm)

    # Part 2 - Plotting
    plt.figure(figsize=(10,7), dpi=300)

    # Categories for number of raw counts
    # disease > control 
    # control > disease
    # disease == control
    df['control'] = df.loc[:, C_id].mean(axis=1)
    df['disease'] = df.loc[:, D_id].mean(axis=1)
    df['color'] = 'white'

    df.loc[df['control']> df['disease'], 'color'] = 'tomato'
    df.loc[df['control']== df['disease'], 'color'] = 'grey'
    df.loc[df['control']< df['disease'], 'color'] = 'gold'

    # Add the labels
    genes = df['GeneName'].iloc[(x[:,0] > greaterX) & (abs (x[:,1]) >greaterY) ]
    names = df['GeneName'].to_list()
    for name in genes:
        # Get the index of the name
        i = names.index(name)
        # Add the text label
        labelpad = 0.01   # Adjust this based on your dataset
        plt.text(x[i,0]+labelpad, x[i,1]+labelpad, name, fontsize=9)

    plt.scatter(x[np.where (df['color']=='tomato'), 0],x[np.where (df['color']=='tomato'), 1], color='tomato', label='WT mice > KO mice', alpha=0.5)
    plt.scatter(x[np.where (df['color']=='gold'), 0],x[np.where (df['color']=='gold'), 1], color='gold', label='WT mice < KO mice', alpha=0.5)
    plt.scatter(x[np.where (df['color']=='grey'), 0],x[np.where (df['color']=='grey'), 1], color='grey', label='WT mice == KO mice', alpha=0.5)

    plt.xlabel('PC 1 (%.2f%%)' % (pca.explained_variance_ratio_[0]*100), fontsize = 18)
    plt.ylabel('PC 2 (%.2f%%)' % (pca.explained_variance_ratio_[1]*100), fontsize=18) 
    plt.title('PCA on counts of RES per gene', fontsize = 25)
    plt.legend(prop={'size': 10}, title= 'Higher Number of RES', title_fontsize='large')
    plt.savefig(filename)
    return df


def get_sig_EP (data, tomato_name='MiceEP-WT>KO.bed', gold_name = 'MiceEP-WT<KO.bed'):
    '''
    Generate BED file for RNA editing sites with differential
    edited proportion
    '''
    clean = {}
    for res in data:
        if data[res]['Pval'] < 0.05:
            clean[res] = data[res]
    
    clean_tomato = [res for res in clean if clean[res]['condA']['mean'] < clean[res]['condB']['mean']]
    clean_gold = [res for res in clean if clean[res]['condA']['mean'] > clean[res]['condB']['mean']]
    
    def save2Bed(l, name):
        with open (name, 'w') as out:
            for res in l:
                case = res.split('|')
                strand = '+'
                if 'TC' in res:
                    strand = '-'
                line = [case[0], case[1], case[1], res, '0', strand]
                out.write('\t'.join(line)+'\n')
        return
    
    save2Bed(clean_tomato, tomato_name)
    save2Bed(clean_gold, gold_name)
    
    return clean, clean_tomato, clean_gold



# Other functions (Might be useful)

def statistical_analysis_old (condA, condA_id, condB, condB_id, cutoffA, cutoffB):
    '''
    Statistical analysis for each RNA editing position. 
    Input:
        - condA / condB: dictionary with the format generated by the function process_output
        - cutoffA / cutoffb: minimum number of samples that need to have RES present to apply 
        the binomial test
    
    For each RNA editing site a binomial test is perform that compares the number of As/Ts and Gs/Cs 
    between condA and condB. For instance, consider the position chr1|58885457|AG where we have the
    following:
        df = {edit: [1,0,5,1], ref:[4,0,6,0], condition: [KO,KO,WT,WT]}
    
    In condA (KO) we have 1 sample with G and A counts. In the condB (WT) we have 2 samples with G and A
    counts. Using the r function (We use a module that allow the use of R functions from python) GLM 
    we apply the binomial test and return the Pvalue. 
    
    In case of having only counts in 1 condition we set the pvalue to np.nan. The binomial test is apply
    only if both condA and condB has counts from at least 1 sample.
    
    NOTE CONDB SHOULD ALWAYS BE WT SAMPLES WHILE CONDA THE DISEASE STATE OR KO CONDITION
    
    Output: 
        - dict_data: dictionary where key == RES_ids and values == Pvalue

    '''
    s0 = time.time() # Check time of the analysis
    
    # Get Total Number of Unique RES in Both conditions
    total_res = set(condA.keys())
    total_res.update (list (condB.keys()))
    
    # Local Variables
    cont = 0
    #pvalues = []
    #RES_ids = []
            
    entries = ['mean', 'std']
    conditions = ['condA', 'condB']
    big_data = {res: {cond:
                      {entry:np.nan for entry in entries} 
                      for cond in conditions} 
                for res in total_res}
    
    condA_s = len(condA_id) # Get number of samples in condA
    condB_s = len (condB_id) # Get number of samples in condB
    total_samples = condA_s + condB_s 
    condition = ['KO']*condA_s + ['WT']*condB_s
    
    print ('\nThere are {} RES in total. Proceeding with statistical analysis...'.format(len(total_res)))
    
    for res in total_res:
        condA_present, condB_present = 0, 0 # To check if RES is present in both conditions
        A_pass, B_pass= 0,0
        
        ### Part 1 - ref and edt are a lists of counts for the ref (A/T) and edited (G/C) base counts
        
        # Define an empty list of counts 
        ref = [0] * total_samples
        edt = [0] * total_samples
        
        # Update with condA counts 
        if res in condA:
            
            # Save mean and standard deviation of res
            big_data[res]['condA']['mean']= np.mean(np.array(condA[res]['EP'])[np.array(condA[res]['EP']) !=0])
            big_data[res]['condA']['std']= np.std(np.array(condA[res]['EP']) [np.array(condA[res]['EP']) !=0])
            
            condA_present = 1
            if 'AG' in res:
                ref[:condA_s] = condA[res]['A_count']
                edt[:condA_s] = condA[res]['G_count']
            else:
                ref[:condA_s] = condA[res]['T_count']
                edt[:condA_s] = condA[res]['C_count']
                
            # Extra check to decide on doing binomial test
            checkA = condA[res]['Strand']
            if len(checkA) - checkA.count('-') >= cutoffA:
                A_pass = 1
            

        # Update with condB counts
        if res in condB:
            
            # Save mean and standard deviation of res
            big_data[res]['condB']['mean']= np.mean(np.array(condB[res]['EP'])[np.array(condB[res]['EP']) != 0])
            big_data[res]['condB']['std']= np.std(np.array(condB[res]['EP'])[np.array(condB[res]['EP']) != 0])
            
            condB_present = 1
            if 'AG' in res:
                ref[condA_s:] = condB[res]['A_count']
                edt[condA_s:] = condB[res]['G_count']
            else:
                ref[condA_s:] = condB[res]['T_count']
                edt[condA_s:] = condB[res]['C_count']
                
            
            # Extra check to decide on doing binomial test
            checkB = condB[res]['Strand']
            if len(checkB) - checkB.count('-') >= cutoffB:
                B_pass = 1
        
        df = pd.DataFrame ({'edt':edt, 'ref':ref, 'condition':condition})
    
        ### Part 2 - Statistical analysis
        if condB_present + condA_present == 2 and A_pass + B_pass == 2:
            # Using R functions
            with (ro.default_converter + pandas2ri.converter).context():
                df2 = ro.conversion.get_conversion().py2rpy(df)
            my_glm = glm('cbind(edt, ref) ~ condition', data=df2, family='binomial')
            pval = float (ro.r.summary(my_glm).rx2("coefficients")[7])

        else:
            pval = np.nan # RES is only present in 1 condition or not sufficient samples for the test

        #pvalues.append(pval)
        #RES_ids.append(res)

        big_data[res]['Pval']=  pval

        cont +=1 
        #if cont % 50000 == 0:
        #    print ('50,000 checked')
        
        if cont % 50000 ==0:
            print ('{} So far...'.format(cont))

    s1 = time.time ()
    
    print ('Time needed to process {} samples, {} s'.format(total_samples, round (s1-s0, 5)))
    
    #dict_data = dict(zip(RES_ids, pvalues))
    return big_data


def significant_plot (data_clean, subdata_clean, condA_name, condB_name, filename):
    '''
    Scatterplot of Editing Propotion (EP) of 2 conditions (condA Vs conB). Significant values are colored red.
    
    Input:
        - data_clean: dictionary generated by the function statistical analysis. Format example:
            {RES_id:
                condA:
                    mean: 0.1 # mean Editing proportion
                    std: 0.005 # standard deviation of the mean
                condB:
                    mean: 0.2
                    std: 0.1
                Pval: 0.6
            }
    
    The input data can be pre-process to remove positions with a pvalue of nan. These positions correspond to RNA
    editing sites only present in 1 condition. 
    '''

    pvalues = [subdata_clean[res]['Pval'] for res in subdata_clean]
    newpvals = fdrcorrection(pvalues, alpha=0.05, method='indep', is_sorted=False)[0]

    print ('After FDR correction, there are {} significant RES'.format(len (np.where(newpvals == True)[0])))
    for idx, res in enumerate (subdata_clean):
        subdata_clean[res]['NewPval'] = newpvals[idx]
    
    # Local Variables
    A_sig, B_sig = [], [] 
    A_notsig, B_notsig = [], []
    
    RES_sig = []


    # Pre-plotting - Separate Significant and non-significant RES for each condition
    for res in subdata_clean:
        if subdata_clean[res]['NewPval'] == True:
            A_sig.append(subdata_clean[res]['condA']['mean'])
            B_sig.append(subdata_clean[res]['condB']['mean'])
            RES_sig.append(res)
        else:
            A_notsig.append(subdata_clean[res]['condA']['mean'])
            B_notsig.append(subdata_clean[res]['condB']['mean'])
    
    for res in data_clean:
        if data_clean[res]['Pval'] in [-1, -2]:
            RES_sig.append(res)

    # Plot
    plt.figure(figsize=(21, 18), dpi=300)
    plt.title('Scatterplot of the Average Editing Proportion', fontsize=25)
    plt.ylabel(condA_name, fontsize = 20)
    plt.xlabel(condB_name, fontsize = 20)
    
    plt.scatter(B_sig, A_sig, color='tomato', label='Significant')
    plt.scatter(B_notsig, A_notsig, color='grey', alpha=0.5, label='Not Significant')
    plt.axline((0.1,0.1), slope=1, color='black', lw=2)
    
    plt.legend(prop={'size': 25})
    plt.savefig(filename)
    return  RES_sig


def pca_significant (df, D_id, C_id, adj=False):
    '''
    PCA analysis

    - Input:
        df: df generated by principal_component_analysis
        D/C_id: the sample ID for disease and control condition
        adj = if the Pvalue adjusted with FDR should be used

    In this PCA analysis, only genes that have counts in both conditions and these counts are not equal
    are kept. Then, a t-test is perform between the counts in disease and control. If the adj is set to True, the FDR
    correction is applied.
    The coloring of the PCA is based on the Pvalue, using a significance level of 0.05
    '''
    
    # Part 1 - Filter and statistical test
    # Remove Genes without RES in both control and disease
    # Remove Genes where the number of RES in control == disease
    df_clean =  df[(df['control'] != 0) & (df['disease'] !=0) & (df['control'] != df['disease'])]
    df_clean['Pval'] = ttest_ind (df_clean[C_id], df_clean[D_id], axis=1,  alternative='two-sided')[1]
    df_clean['Pval.adj'] = fdrcorrection(ttest_ind (df_clean[D_id], df_clean[C_id], axis=1, alternative='two-sided')[1])[1]
    
    
    # Part 2 - PCA Analysis
    filt = D_id + C_id
    df_tmp = df_clean[filt]
    ## Normalize
    scaling=StandardScaler()
    df_norm = scaling.fit_transform(df_tmp)
    df_norm = pd.DataFrame(df_norm)
    ## PCA
    pca = PCA(n_components=df_norm.shape[1])
    pca.fit(df_norm)
    x = pca.transform(df_norm)
    
    # Part 3 - Plotting
    plt.figure(figsize=(10,10), dpi=300)

    # Categories for number of raw counts
    # Significant  (Pval < 0.05)
    # Not Significant (Pval > 0.05)
    df_clean['color'] = 'white'

    if adj:
        df_clean.loc[df_clean['Pval.adj']> 0.05, 'color'] = 'grey'
        df_clean.loc[df_clean['Pval.adj']< 0.05, 'color'] = 'tomato'
    else:
        df_clean.loc[df_clean['Pval']> 0.05, 'color'] = 'grey'
        df_clean.loc[df_clean['Pval']< 0.05, 'color'] = 'tomato'

    # Add the labels
    genes = df_clean['GeneName'].iloc[(x[:,0] > 5) & (abs (x[:,1]) >2) ]
    names = df_clean['GeneName'].to_list()
    for name in genes:
        # Get the index of the name
        i = names.index(name)
        # Add the text label
        labelpad = 0.01   # Adjust this based on your dataset
        plt.text(x[i,0]+labelpad, x[i,1]+labelpad, name, fontsize=9)

    plt.scatter(x[np.where (df_clean['color']=='grey'), 0],x[np.where (df_clean['color']=='grey'), 1], color='grey', label='Not Significant', alpha=0.5)
    plt.scatter(x[np.where (df_clean['color']=='tomato'), 0],x[np.where (df_clean['color']=='tomato'), 1], color='tomato', label='Significant', alpha=0.5)

    plt.xlabel('PC 1 (%.2f%%)' % (pca.explained_variance_ratio_[0]*100), fontsize = 17)
    plt.ylabel('PC 2 (%.2f%%)' % (pca.explained_variance_ratio_[1]*100), fontsize=17) 
    plt.title('PCA on Raw counts of RES', fontsize = 20)
    plt.legend()
    return df_clean

   
def saving_forR (data, data_sig, name, path='./', greater = 'control'):
    '''
    Save into Bed format significant RES
    
    - Input:
        data: dictionary generated by the function statistical_analysis
        data_sig: dictionary generated by the function significant_plot
        name: name of the file 
        path: path to the file
        greater:
            If greater == control >> The output file will have RES whose mean Editing proportion is 
            greater than the disease state
    
    - Output: 
        File in Bed format. Each line contain information of a RES.
    '''
    # Local variables
    forR = []
    
    # Analysis
    
    for res in data_sig:
        ## CondB is always control. 
        ## Keep RES where Editing prortion in control is greater than in disiease
        
        if greater == 'control':
            if data[res]['condB']['mean'] > data[res]['condA']['mean']:
                forR.append(res) 
        else:
            if data[res]['condB']['mean'] < data[res]['condA']['mean']:
                forR.append(res) 
    
    print ('There are {} significant RES'.format(len(forR)))
    
    # Save into file
    
    with open (os.path.join(path, name), 'w') as out:
        for res in forR:
            line = []
            case = res.split('|') 
            line.append (case[0]) # chrom 
            line.append (case[1]) # pos
            line.append (case[1]) # pos
            line.append (res) # name
            line.append('0') # score
            if 'TC' in res:
                line.append ('-') # strand
            else:
                line.append('+')
            newline = '\t'.join (line)
            out.write (newline + '\n')
    return
