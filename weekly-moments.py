'''Makes boxplots of bootstrapped historical weekly flows (pink) and synthetic
weekly flows (blue) as well as their means and standard deviations at
Marietta. Also plots p-values from rank-sum test for differences in the median
between historical and synthetic flows and from Levene's test for differences
in the variance between historical and synthetic flows. The site being plotted
can be changed on line 84.'''

from __future__ import division
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns
from scipy import stats
import os

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
    plt.rcParams['figure.figsize'] = (14, 14)
    plt.rcParams['font.size'] = 18
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.labelsize'] = 1.1*plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.1*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']

init_plotting()

def set_box_color(bp, color):
    '''Sets colors of boxplot elements'''
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color, linestyle='solid')
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color='k')

# thanks to http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
def boxplots(syn, hist, xticks=True, legend=True, loc='upper right'):
  '''Makes boxplots'''
  # bpl = boxplots of synthetic data, bpr = boxplots of bootstrapped historical data
  bpl = plt.boxplot(syn, positions=np.arange(1,53)-0.15, sym='', widths=0.25, patch_artist=True)
  bpr = plt.boxplot(hist, positions=np.arange(1,53)+0.15, sym='', widths=0.25, patch_artist=True)
  set_box_color(bpl, 'lightskyblue')
  set_box_color(bpr, 'lightcoral')

  plt.plot([], c='lightskyblue', label='Synthetic')
  plt.plot([], c='lightcoral', label='Historical') # remember means and stdevs here are bootstrapped
  plt.gca().xaxis.grid(False)

  if xticks:
    points = range(1,13,13)
    plt.gca().set_xticks(points)
    plt.gca().set_xticklabels(points)
  else:
    plt.gca().set_xticks([])
  plt.gca().set_xlim([0,13])

  if legend:
    plt.legend(ncol=2, loc=loc)

  plt.locator_params(axis='y', nbins=5)


# Make statistical validation plots of weekly moments
s = 5       # index of the site to be plotted
space = ['real', 'log']
legend_loc = ['upper right', 'lower left']
site = all_sites[s]
sitename = all_sitenames[s]

H = np.loadtxt('historical-data/' + site + '.csv', delimiter=',')
S = np.loadtxt('synthetic-data-stat/' + site + '_SYN60.csv', delimiter=',')     # modify based on stationary or dynamic dataset
S = S.reshape((np.shape(S)[0], int(np.shape(S)[1]/52),52)) # n_realizations x n_syn_years x 52

for j in range(len(space)):
    if j == 1:
        H = np.log(H)
        S = np.log(S)
    N = H.shape[0]
    num_resamples = np.shape(S)[0]
    r = np.random.randint(N, size=(N, num_resamples))

    fig = plt.figure()

    # Plot boxplot of weekly totals from n_realizations*n_syn_years
    # and all historical years
    ax = fig.add_subplot(5,1,1)
    boxplots(S.reshape((np.shape(S)[0]*np.shape(S)[1], 52)), H, xticks=False, legend=True,
                loc=legend_loc[j])
    if j == 0:
        ax.set_ylabel('Q ($10^{6}$ week)')
        ax.set_yticks(np.arange(0, 55000, 10000))
    else:
        ax.set_ylabel('Log(Q)')
        ax.set_yticks(np.arange(5, 20, 5))
        
    ax = fig.add_subplot(5,1,2)
    # get weekly means throughout all n_syn_years and n_realizations
    boxplots(S.mean(axis=1), H[r].mean(axis=0), xticks=False, legend=False)
    ax.set_ylabel('$\hat{\mu}_Q$')
    
    if j == 1:
        ax.set_yticks(np.arange(8, 11, 1))
    else:
        ax.set_yticks(np.arange(0, 40000, 5000))

    ax = fig.add_subplot(5,1,3)
    # get weekly std deviations throughout all n_syn_years and n_realizations
    boxplots(S.std(axis=1), H[r].std(axis=0), xticks=False, legend=False)
    ax.set_ylabel('$\hat{\sigma}_Q$')
    if j == 1:
        ax.set_yticks(np.arange(0, 3, 1))
    else:
        ax.set_yticks(np.arange(0, 35000, 5000))

    # Wilcoxon's rank-sum test for weekly medians and Levene's test for weekly variances
    # Ideally Wilconxon's p ~= 1.0
    # Ideally Levene's p ~= 1.0
    wilcoxon_pvals = np.zeros(52)
    levene_pvals = np.zeros(52)

    for i in range(52):
        wilcoxon_pvals[i] = stats.ranksums(H[:,i], S.reshape((np.shape(S)[0]*np.shape(S)[1], 52))[:,i])[1]
        levene_pvals[i] = stats.levene(H[:,i], S.reshape((np.shape(S)[0]*np.shape(S)[1], 52))[:,i])[1]

    ax = fig.add_subplot(5,1,4)
    ax.bar(np.arange(1,53)-0.4, wilcoxon_pvals, facecolor='0.7', edgecolor='None')
    ax.set_xlim([0,53])
    ax.plot([0, 54], [0.05, 0.05], color='k')
    ax.set_xticks([])
    ax.set_ylabel('Wilcoxon $p$')
    ax.set_yticks(np.arange(0, 1.5, 0.5))

    ax = fig.add_subplot(5,1,5)
    ax.bar(np.arange(1,53)-0.4, levene_pvals, facecolor='0.7', edgecolor='None')
    ax.set_xlim([0,53])
    ax.plot([0, 54], [0.05, 0.05], color='k')
    ax.set_xticks([])
    ax.set_ylabel('Levene $p$')
    ax.set_yticks(np.arange(0, 1.5, 0.5))

    if j == 0:
        fig.suptitle('Real space ' + '(' + sitename + ')')
    else:
        fig.suptitle('Log space ' + '(' + sitename + ')')

    fig.tight_layout()
    fig.savefig('figures/moments_pvalues_' + space[j] + '_' + sitename + '-stat.pdf')   # modify based on stationary or dynamic dataset

    fig.clf()
