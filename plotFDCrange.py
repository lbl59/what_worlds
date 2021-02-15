'''Plots the ranges spanned by historical (pink) and synthetic (blue) flow
duration curves for each year of the historical and synthetic records at all
sites. 

The list of sites can be changed on line 11 and 15.

The assumption of stationarity can be changed on lines 101, 103 and 106. '''

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os

all_sites = ['trainingLittleRiverRaleighInflow','trainingOWASAInflow','trainingClaytonGageInflow',
             'trainingCrabtreeCreekInflow','trainingFallsLakeInflow','trainingJordanLakeInflow',
             'trainingLakeWBInflow','trainingLillingtonInflow','trainingLittleRiverInflow','trainingMichieInflow']

all_sitenames = ['Little River Raleigh', 'OWASA', 'Clayton', 'Crabtree Creek', 'Falls Lake', 'Jordan Lake',
                 'Lake Wheeler/Benson', 'Lillington', 'Little River', 'Lake Michie']

# https://justgagan.wordpress.com/2010/09/22/python-create-path-or-directories-if-not-exist/
# a just-in-case step
def assure_path_exists(path):
    '''Creates directory if it doesn't exist'''
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)

# getcwd() returns the current working directory of a process
assure_path_exists(os.getcwd() + '/figures/')

def init_plotting():
    '''Sets plotting characteristics'''
    sns.set_style('whitegrid')
    plt.rcParams['figure.figsize'] = (20,12)  # changes made here
    plt.rcParams['font.size'] = 18
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.labelsize'] = 1.1*plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.1*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']

init_plotting()

# FDC: flow duration curve
def plotFDCrange(Qdaily_syn, Qdaily_hist, sites):

    n = 52
    M = np.array(range(1, n+1))
    P = (M-0.5)/n

    fig = plt.figure()
    for j in range(len(sites)):
        # load data for site j
        histData_j = Qdaily_hist[:, j]
        synData_j = Qdaily_syn[:, j]

        histData_newshape = (len(histData_j)//n, n)
        synData_newshape = (len(synData_j)//n, n)
        q_hist = np.reshape(histData_j, histData_newshape)
        q_syn = np.reshape(synData_j, synData_newshape)

        fdc_hist = np.empty(np.shape(q_hist))
        fdc_hist[:] = np.NaN
        for k in range(np.shape(fdc_hist)[0]):
            fdc_hist[k,:] = np.sort(q_hist[k,:],0)[::-1]
        
        fdc_syn = np.empty(np.shape(q_syn))
        fdc_syn[:] = np.NaN
        for k in range(np.shape(fdc_syn)[0]):
            fdc_syn[k,:] = np.sort(q_syn[k,:], 0)[::-1]

        ax = fig.add_subplot(2,5,j+1)
        # plotting the borders of the historical and synthetic FDC curves
        ax.semilogy(P, np.min(fdc_syn, 0), c='lightskyblue', label='Synthetic')
        ax.semilogy(P, np.max(fdc_syn, 0), c='lightskyblue', label='Synthetic')
        ax.semilogy(P, np.min(fdc_hist, 0), c='lightcoral', label='Historical')
        ax.semilogy(P, np.max(fdc_hist, 0), c='lightcoral', label='Historical')
        if j == 0 or j == 5:
            ax.set_ylabel('Q ($10^{6}$ gal/week)')

        ax.fill_between(P, np.min(fdc_syn, 0), np.max(fdc_syn, 0), color='lightskyblue')
        ax.fill_between(P, np.min(fdc_hist, 0), np.max(fdc_hist, 0), color='lightcoral')
        
        ax.set_title(sites[j], fontsize=18)
        ax.tick_params(axis='both', labelsize=14)

        if j <= 1:
            ax.tick_params(axis='x', labelbottom='off')

        handles, labels = plt.gca().get_legend_handles_labels()
        labels, ids = np.unique(labels, return_index=True)
        handles = [handles[i] for i in ids]
        plt.grid(True, which='both', ls='-')

    fig.subplots_adjust(bottom=0.2, wspace=0.6)
    fig.legend(handles, labels, fontsize=18, loc='lower center', ncol=2, 
               frameon=True, bbox_to_anchor=(0.5, 0.07))
    fig.text(0.5, 0.14, 'Probability of exceedance', ha='center', size=20)
    fig.suptitle('Flow duration curves assuming extreme flooding', fontsize=25)     # modify based on stationary or dynamic dataset
    plt.subplots_adjust(top=0.9)
    fig.savefig('figures/FDCs-dyn.pdf')     # modify based on stationary or dynamic dataset
    fig.clf()
    
Qdaily_syn = np.loadtxt('synthetic-data-dyn/Qdaily-syn-dyn.csv', delimiter=',')    # modify based on stationary or dynamic dataset
Qdaily_hist = np.loadtxt('historical-data/Qdaily-hist.csv', delimiter=',')

plotFDCrange(Qdaily_syn, Qdaily_hist, all_sitenames)