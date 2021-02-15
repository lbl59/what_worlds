import numpy as np
import os

'''
Reshape each historical and one-year synthetic inflow record into a 1D vector
Concatenate all vectors into a n x 10 matrix (one column for each site)
Output the Qdaily matrix for historical and synthetic flows for FDC validation.

'''

hist_datadir = "historical-data"
syn_datadir = "synthetic-data-stat"     # modify depending on stationary or dynamic dataset

# modify based on test case
nweeks = 52
nsites = 10
nhist_years = 81
nsyn_real = 1000

Qdaily_hist = np.zeros((nhist_years*nweeks, nsites), dtype=float)
Qdaily_syn = np.zeros((nsyn_real*nweeks, nsites), dtype=float)

# load the historical and synthetic data
# Thanks to https://stackoverflow.com/questions/33503993/read-in-all-csv-files-from-a-directory-using-python
hist_count = 0
for filename in os.listdir(hist_datadir):
    if filename.endswith("Inflow.csv"):
        filepath = hist_datadir + '/' + filename
        q_file = np.loadtxt(filepath, delimiter=',')
        q_file2 = np.reshape(q_file, (q_file.shape[0]*q_file.shape[1],))
        Qdaily_hist[:, hist_count] = q_file2
        hist_count = hist_count + 1
        print("hist_count = ", hist_count)

syn_count = 0
for filename in os.listdir(syn_datadir): 
    if filename.endswith("_SYN01.csv"):
        filepath = syn_datadir + '/' + filename
        q_file = np.loadtxt(filepath, delimiter=',')
        q_file2 = np.reshape(q_file, (q_file.shape[0]*q_file.shape[1],))
        Qdaily_syn[:, syn_count] = q_file2
        syn_count = syn_count + 1
        print("syn_count = ", syn_count)

hist_datapath = hist_datadir + "/Qdaily-hist.csv"  # modify depending on stationary or dynamic dataset
syn_datapath = syn_datadir + "/Qdaily-syn-stat.csv"     # modify depending on stationary or dynamic dataset
np.savetxt(hist_datapath, Qdaily_hist, delimiter=',')
np.savetxt(syn_datapath, Qdaily_syn, delimiter=',')
