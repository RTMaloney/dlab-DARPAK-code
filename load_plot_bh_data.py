# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 11:12:06 2018

@author: RyanDataPC
"""

import matplotlib as mpl 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# import the behavioural data excel sheet as a pandas dataframe.
# We ask it to import the first sheet in the file, 'Trial by trial'

# The 'r' before the path indicates a raw string. We want the first sheet only: 'Trial by trial'
bh_data = pd.read_excel(r'C:\Users\RyanDataPC\Dropbox\DARPAK\material\Behavioural_Data.xlsx', 
                        sheet_name='Trial by trial')

# Pull out just the averages for plotting. They are in row 23
# There are missing data for some subjects after trial 15, so we'll only plot those
rewrd_trials = bh_data.iloc[23,1:16]   # The 15 reward trials first
instr_trials = bh_data.iloc[23,26:41]  # The 15 instructive trials

# Also pull out the standard errors; they are row 24
rewrd_trials_se = bh_data.iloc[24,1:16]
instr_trials_se = bh_data.iloc[24,26:41] 

# Now make the line plots in a new figure
fig = plt.figure()
# Convert performance data into a percentage contrast when plotting
plt.plot(np.arange(1,16), rewrd_trials*100, color='r') 
plt.plot(np.arange(1,16), instr_trials*100, color='k')

# Adjust the axis limits
plt.xlim(0, 16)
plt.ylim(0,30)

plt.errorbar(np.arange(1,16), rewrd_trials*100, yerr = rewrd_trials_se, color='r')


