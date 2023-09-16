from pylab import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

# Read the transition output files:

transition_btw_states = pd.read_csv("Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_networkStates.csv", header=None)
print(transition_btw_states.head())
print(transition_btw_states.shape)

xthicklabels = ['ASCL1', 'FOXA1', 'FOXA2', 'ELF3', 'RBPJ', 'FLI1',
                'SMAD4', 'NR0B2', 'NR0B1', 'BCL3', 'STAT6', 'ISL1',
                'SOX11', 'CEBPD', 'EBF1', 'TCF4', 'RCOR2', 'TCF3',
                'NEUROD2', 'OLIG2', 'MITF', 'SIX5', 'TEAD4', 'ZNF217',
                'KLF2', 'GATA4', 'REST']

fig, ax = plt.subplots(figsize=(6.5, 9))
plt.subplot(2,1,1)
sns.heatmap(transition_btw_states.iloc[0:50], cmap="Reds", vmin=0, vmax=1, cbar=False, linewidths=0.1, linecolor='gray')
plt.xticks([])
plt.subplot(2,1,2)
sns.heatmap(transition_btw_states.iloc[80567:], cmap="Reds", vmin=0, vmax=1, cbar=False, linewidths=0.1, linecolor='gray', xticklabels=xthicklabels)

fig.text(0.01, 0.5, 'Asynchronous iterations', va='center', rotation='vertical')
fig.savefig('NON-NE to NE transition_summary of the steps', dpi=600)

fig, ax = plt.subplots(figsize=(2, 9))
sns.heatmap(pd.DataFrame(transition_btw_states.iloc[0].transpose()), cmap="Reds", vmin=0, vmax=1, cbar=False, linewidths=0.1, linecolor='gray', yticklabels=xthicklabels, square=True)
plt.xticks([])
plt.title('Subtype: \n NON-NE')
fig.savefig('NON-NE to NE transition_initialState', dpi=600)

fig, ax = plt.subplots(figsize=(2, 9))
sns.heatmap(pd.DataFrame(transition_btw_states.iloc[-1].transpose()), cmap="Reds", vmin=0, vmax=1, cbar=False, linewidths=0.1, linecolor='gray', yticklabels=xthicklabels, square=True)
plt.xticks([])
plt.title('Subtype: \n NE')
fig.savefig('NON-NE to NE transition_lastState', dpi=600)

plt.show()