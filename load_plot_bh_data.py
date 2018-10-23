# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 11:12:06 2018

@author: RyanDataPC
"""

#%%
import matplotlib as mpl 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.io as sio

n = 22 #we have modelling results for 22 subjects.

#%%
# import the behavioural data excel sheet as a pandas dataframe.
# We ask it to import the first sheet in the file, 'Trial by trial'

# The 'r' before the path indicates a raw string. We want the first sheet only: 'Trial by trial'
bh_data = pd.read_excel(r'C:\Users\RyanDataPC\Dropbox\DARPAK\material\Behavioural_Data.xlsx', 
                        sheet_name='Trial by trial')

#%%
# Pull out just the averages for plotting. They are in row 23
# There are missing data for some subjects after trial 15, so we'll only plot those
rewrd_trials = bh_data.iloc[23,1:16]   # The 15 reward trials first
instr_trials = bh_data.iloc[23,26:41]  # The 15 instructive trials

# Also pull out the standard errors; they are row 24
rewrd_trials_se = bh_data.iloc[24,1:16]
instr_trials_se = bh_data.iloc[24,26:41] 

#%%
# Now make the line plots in a new figure
fig = plt.figure()
# Convert performance data into a percentage contrast when plotting
plt.errorbar(np.arange(1,16), rewrd_trials*100, yerr = rewrd_trials_se*100, 
             color = 'r')
plt.errorbar(np.arange(1,16), instr_trials*100, yerr = instr_trials_se*100,
             color = 'k')

# Adjust the axis limits
plt.xlim(0,16)
plt.ylim(0,30)

# Adjust the tick labels on the x axis
plt.xticks(np.arange(0,16, 5))


#%%  ***** Now load and plot the initial modelling results. *****
#
# These have been processed from matlab data in the try_to_unpack_bh_data_DARPAK.m script, in the form of a .mat file.

# Load the mat file using scipy.io functionality:
model_results = sio.loadmat (r'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\behavioural\DARKPAK_proc_model_results.mat')
# The mat file is loaded as a numpy ndarray.

# Compute mean MI across subjects for instructive condition
MI_mean_instr = np.nanmean(model_results['MI']['mean_instructive'][0][0],axis=0)
MI_mean_monet = np.nanmean(model_results['MI']['mean_monetary'][0][0],axis=0)

# And now do the same for entropy
entropy_mean_instr = np.nanmean(model_results['entropy']['mean_instructive'][0][0],axis=0)
entropy_mean_monet = np.nanmean(model_results['entropy']['mean_monetary'][0][0],axis=0)

# Now compute the standard errors of these values, for plotting.
MI_SE_instr = np.nanstd(model_results['MI']['mean_instructive'][0][0],axis=0) / np.sqrt(n)
MI_SE_monet = np.nanstd(model_results['MI']['mean_monetary'][0][0],axis=0) / np.sqrt(n)

entropy_SE_instr = np.nanstd(model_results['entropy']['mean_instructive'][0][0],axis=0) / np.sqrt(n)
entropy_SE_monet = np.nanstd(model_results['entropy']['mean_monetary'][0][0],axis=0) / np.sqrt(n)


#%% Plot the modelling results.
# Only plot the first 15 trials, because they are missing for some subjects for trials >15
fig = plt.figure(figsize=(8*1.6,8)) #open a new figure and set its size (values in inches)
ax = fig.subplots(2,1) # , sharex=True)

# Plot entropy.
# To insert a legend that does not include error bars, we need to name each parameter 
# generated when constructing the errorbar plots:
# see here for how: http://scientificpythonsnippets.com/index.php/distributions/4-scientific-plotting-with-python-plot-with-error-bars-using-pyplot
mon_line, bars, caps = ax[0].errorbar(np.arange(1,16), entropy_mean_monet[0:15], yerr = entropy_SE_monet[0:15], 
             color = 'r')
ins_line, bars, caps = ax[0].errorbar(np.arange(1,16), entropy_mean_instr[0:15], yerr = entropy_SE_instr[0:15], 
             color = 'k')

# then label the relevant parameters, namely the lines, accordingly
plt.setp(mon_line,label='Monetary')
plt.setp(ins_line,label='Instructive')

# Now we can insert this legend; only need 1 for the 2 subplots:
ax[0].legend()

# Plot mutual information
ax[1].errorbar(np.arange(1,16), MI_mean_monet[0:15], yerr = MI_SE_monet[0:15], 
             color = 'r')
ax[1].errorbar(np.arange(1,16), MI_mean_instr[0:15], yerr = MI_SE_instr[0:15], 
             color = 'k')

ax[1].set_ylim([0,0.45])
ax[1].set_xlim([0,16])

ax[0].set_xlim([0,16])

# Put in labels for the axes
ax[1].set_xlabel('Trial Number')
ax[0].set_ylabel('Belief Uncertainty (a.u.)')
ax[1].set_ylabel('Belief Update Size (a.u.)')

