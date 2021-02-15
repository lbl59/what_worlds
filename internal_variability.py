from __future__ import division
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns
from scipy import stats
import os

'''
Plots boxplots characterizing the internal variability as a function of the number of realization in
a synthetic inflow of site s.
1) Obtain the 25% lowest/highest flows (depending on drought/flood scenario) 
2) Obtain the mean and std deviation of the 25% flows across n-years of the synthetic record
3) Obtain the mean and std deviation of the means and std deviations across m-realizations
4) Repeat for increasing numbers of realizations (100,200,...1000)

Site to be printed can be changed on line 171.
'''

all_sites = ['trainingLittleRiverRaleighInflow','trainingOWASAInflow','trainingClaytonGageInflow',
             'trainingCrabtreeCreekInflow','trainingFallsLakeInflow','trainingJordanLakeInflow',
             'trainingLakeWBInflow','trainingLillingtonInflow','trainingLittleRiverInflow','trainingMichieInflow']

all_sitenames = ['Little River Raleigh', 'OWASA', 'Clayton', 'Crabtree Creek', 'Falls Lake', 'Jordan Lake',
                 'Lake Wheeler/Benson', 'Lillington', 'Little River', 'Lake Michie']

# https://justgagan.wordpress.com/2010/09/22/python-create-path-or-directories-if-not-exist/
# this is necessary so the workflow does not have to ensure that the figures are plotted in a particular order
# something nice to have
def assure_path_exists(path):
    '''Creates directory if it doesn't exist'''
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)
                
assure_path_exists(os.getcwd() + '/figures/')

def init_plotting():
    '''Sets plotting characteristics'''
    sns.set_style('whitegrid')
    plt.rcParams['figure.figsize'] = (12, 12)
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Source Sans Pro'
    plt.rcParams['axes.labelsize'] = 1.1*plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.1*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']

def set_box_color(bp, color):
    '''Sets colors of boxplot elements'''
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color, linestyle='solid')
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color='k')

'''
Returns a figure with 2 sub-plots: the boxplot of means and the boxplot of std deviations across
increasing realizations.

@param s The index of the site of interest
@param space A string variable denoting real or log space
@param quantile The %lowest/highest flows to examine between [0,1]. Default value is 1.0 (examine all values).
@param drought True if testing for lowest 25% of flows. False if testing for highest 25% of flows. Default True

'''

def internal_variability(s, space, quantile=1.0, drought=True):
    assure_path_exists(os.getcwd() + '/figures/')
    site = all_sites[s]
    sitename = all_sitenames[s]

    S = np.loadtxt('synthetic-data-stat/' + site + '_SYN60.csv', delimiter=',') # modify depending on stationary or dynamic dataset
    S = S.reshape((np.shape(S)[0], int(np.shape(S)[1]/52),52)) # n_realizations x n_syn_years x 52
    
    # sort all years in each realization by increasing/decreasing inflows
    if space == 'log':
        S = np.log(S)
    for i in range(S.shape[0]):     # for i in range 1000
        S_i = S[i,:,:]
        for j in range(S_i.shape[0]):    # for j in range n_syn_years
            S_ij = S_i[j,:]
            # if drought, sort in order of increasing flows
            # # else if flood, sort in order of decreasing flows
            if (drought == True):
                S_ij = np.sort(S_ij)
            else:
                S_ij = np.sort(S_ij)[::-1]
            S_i[j,:] = S_ij
        S[i,:,:] = S_i

    S25 = S[:, :,:int(S.shape[2]*quantile)]
    idx = np.arange(50, 1050, 50)   # realization indices
    S25plot = np.zeros((len(idx), S25.shape[1], S25.shape[2]), dtype=float)
    for m in range(len(idx)):
        S25plot[m,:,:] = S25[idx[m]-1,:,:]
    S25plot = np.reshape(S25plot, (S25plot.shape[1]*S25plot.shape[2], S25plot.shape[0]))
    S25 = np.reshape(S25, (S25.shape[1]*S25.shape[2], S25.shape[0]))
    means_per_real = np.mean(S25, axis=0)
    std_per_real = np.std(S25, axis=0)
    S25_means = np.zeros(len(idx))
    S25_std = np.zeros(len(idx))

    for k in range(len(idx)):
        over = idx[k]
        S25_means[k] = np.mean(means_per_real[:over-1])
        S25_std[k] = np.std(std_per_real[:over-1])
    
    # make the plots
    fig = plt.figure()
    ax = fig.add_subplot(3,1,1)    # box plots
    bpl_syn = plt.boxplot(S25plot, sym='', widths=0.3, patch_artist=True, meanline=True)
    set_box_color(bpl_syn, 'lightskyblue')
    plt.plot([], c='lightskyblue')      
    
    if space == 'log':
        plt.ylabel('Log of inflow', fontsize=8)
        plt.ylim([0,10])
    else:
        plt.ylabel('Inflow ($10^{6}$ gal/week)', fontsize=8)
        plt.yticks(np.arange(0,30000, 5000))
    
    plt.xticks([])
    ax.set_xticklabels([])
    plt.grid(axis='y')
    ax_title = 'default'

    if drought == True:
        ax_title = 'Range of the lower 25% of annual ' + space + ' flows'
    else:
        ax_title = 'Range of the higher 25% of annual ' + space + 'flows' 
        
    plt.title(ax_title, fontsize=11)

    ax = fig.add_subplot(3,1,2)    # mean convergence
    plt.plot(idx, S25_means, 'b-')
    if space == 'log':
        plt.ylim([0,10])
    else:
        plt.yticks(np.arange(6700, 7100, 100))
    plt.xticks([])
    ax.set_xticklabels([])
    plt.grid(axis='y')
    plt.ylabel('$\mu$', fontsize=8)
    plt.title('$\mu$ through n-realizations', fontsize=11)
        
    ax = fig.add_subplot(3,1,3)    # std dev convergence
    plt.plot(idx, S25_std, 'b-')    
    if space == 'log':
        plt.ylim([0,0.2])
    else:
        plt.yticks(np.arange(1200, 2500, 500))
    #plt.xticks(np.arange(0, 20))
    #ax.set_xticklabels(idx)
    plt.ylabel('$\sigma$', fontsize=8)
    plt.grid(axis='y')
    plt.xlabel('Realizations', fontsize=8)
    plt.title('$\sigma$ through n-realizations', fontsize=11)
    
    sns.despine(left=True)
    
    main_title = sitename + ' ' + space + '-space inflows'
    plt.suptitle(main_title, fontsize=14)
    #fig.tight_layout()
    fig.subplots_adjust(top=0.9,wspace=0.1, hspace=0.5)
    
    fig_name = 'figures/internal-variability-' + space + '-stat.pdf'     # modify depending on stationary or dynamic dataset
    fig.savefig(fig_name)
    fig.clf()

internal_variability(5, 'log', quantile=1.0, drought=False)