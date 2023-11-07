# Script to generate the following figures

# washout, AC jumps + nuk dynamics (in percentage)
# cumulative distribution: slippage vs exit
# + percentage (slippage vs exit)

# APC-A washout
# ecdf, slippage vs exit, apca vs wt

# GAL Mad2 washout
# ECDF + percentages slippage vs exit

import matplotlib
matplotlib.use('Agg') 
import os # for folder creation in code
import pickle # for saving out large res files 

import stochpy
smod = stochpy.SSA()
smod.DoStochSim()
smod.PlotSpeciesTimeSeries

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.pyplot import cm
import pandas as pd

plt.rcParams['figure.dpi'] = 100

smod.Model('model_washout.psc', dir='./model_files/')

df_ind = ['Condition', 'Slippage (%)', 'Exit (%)', 'Nr of adapters before washout']
df_washout = pd.DataFrame(columns=df_ind)

## Functions
# function to define the number of unattached kineotchores and the time point the SAC on - SAC off switch happens. We follow AC transition from the lower 
# (SAC on) to upper steady state (SAC off)
# input:
## traj - stochastic simulation result trajectories
## species_col - the column index of the species from species selection, from the steady states of which we calculate the threshold
# output[0] - timepoint of adaptation
# output[1] - number of unattached kinetochores at the time of adaptation
def switch_nuk(traj, species_col, threshold):#, transient_time):
   # traj = traj[transient_time:]
    t_switch_up = np.argmax(traj[:,species_col]>threshold)
    t_switch_down = np.argmax(traj[:,species_col]<threshold)
    if t_switch_up==0:
        t_switch_up=np.inf
        return[np.inf, traj[t_switch_down,0], traj[t_switch_down,2]]
    return [traj[t_switch_up,0], traj[t_switch_up,2]]
   
# function to compute empirical cumulative distribution function
def ecdf(data):
    """ Compute ECDF """
    x = np.sort(data)
    n = x.size
    y = np.arange(1, float(n)+1) / float(n)
    return(x,y)

### Simulation settings

traj_nr = 100 # simulate 100 trajectories corresponding to 100 single cells
time = 1000 # simulate for 1000 time units
rescale = 3.33 #1 time unit will be considered ~ 20 seconds
wo_time=180 # time of washout

# Do stochastic simulation for wild-type
smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", species_selection=['ac','nuk', 'cycb', 'mstar','c','acmc'])
# Save the results
res_wo_wt = list()
for n in np.arange(1, traj_nr+1):
    smod.GetTrajectoryData(n)
    res_wo_wt.append(smod.data_stochsim.getSpecies())

### save res trajectories
f=open('./Results_trajectories/res_wt_wo.pickle', 'w')
pickle.dump(res_wo_wt, f)
f.close()

# Lists to store plot handles and labels for the legend
handles = []
labels = []
fig, ax= plt.subplots()
# plot result trajectories, AC vs time
# color-coding to differentiate between cells that adapt or exit regularly
for i in res_wo_wt:
    if switches_wt_nuk[count][1] < 10: #remove early adapters
        if switches_wt_nuk[count][1] > 1: #the ones which switch steady states before the saddle node (=1) is reached are the adapters
            line, = ax.plot(i[:, 0]/rescale, i[:, 1], color='deeppink', alpha=0.5)
            if not slippage_legend_added:
                handles.append(line)
                labels.append('Slippage')
                slippage_legend_added = True
        else: # the rest which switch steady states with nuk=1 or nuk=0 are the ones which exit regularly
            line, = ax.plot(i[:, 0]/rescale, i[:, 1], color='g', alpha=0.35)
            if not mitotic_exit_legend_added:
                handles.append(line)
                labels.append('Mitotic exit')
                mitotic_exit_legend_added = True

        ax.plot(i[:, 0], i[:, 2] * 3.125, 'k', alpha=0.5)
    count += 1

# Add the "washout" line to the legend
handles.append(plt.axvline(x=wo_time/rescale, color='grey', linestyle='--'))
labels.append('Washout')
ax.set_title('Nocodazole washout at t = 180 min')
ax.set_xlabel('Time (min)')
ax.set_ylabel('APC/CCdc20 (mol./cell) and % of UK(black)')
# set timescale ticks 
start_value = wo_time / rescale
x_values = [start_value] + np.arange(0, 241, 60)
tick_labels = [0,60,120,180,240]
ax.set_xticks(x_values)
ax.set_xticklabels(tick_labels)  
# Create the legend using custom handles and labels
ax.legend(handles, labels)
fig.savefig('./figures/washout/wt_wo_apc.png')  
plt.close(fig)


## create cumulative distribution functions for slippage/regular exit cells
slippage_wt=list()
exit_wt=list()


for i in switches_wt_nuk:
    if i[1]<10:
        if i[1]>1.0:
            slippage_wt.append(i[0])
        else:
            exit_wt.append(i[0])

    fig, ax = plt.subplots()
        
x_slp_wt,y_slp_wt = ecdf(slippage_wt)
ax.step(x=x_slp_wt/rescale, y=y_slp_wt,color='m', label='Slippage WT');

x_exit_wt,y_exit_wt = ecdf(exit_wt)
ax.step(x=x_exit_wt/rescale, y=y_exit_wt, color='g', label='Mitotic exit, WT')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Fraction of adapted cells')
ax.legend()

fig.savefig('/storage/home/bstier/figures/washout/wt_noco_wo_ECDF.png') 
plt.close(fig)

### Count slippage and exit
without_early_adapt = len(slippage_wt)+len(exit_wt)
slippage_percentage = float(len(slippage_wt))/float(without_early_adapt)*100
exit_percentage = float(len(exit_wt))/float(without_early_adapt)*100
earlies_wt = traj_nr-without_early_adapt

## save to a table

row = {'Condition': 'Wild-type', 'Slippage (%)': slippage_percentage, 'Exit (%)': exit_percentage, 'Nr of adapters before washout': earlies_wt}
df_washout = df_washout.append(row, ignore_index=True)

######################################################################################################################

### APC-A WASHOUT
new_kassac = 0.002
smod.ChangeParameter("kassac", new_kassac)

### INITIAL CONDITIONS APCA MUTANT
c = 136
mstar = 34
mc = 113
ac = 8
acmc = 27
m = 1
a = 114
cycb = 1188
initial_conditions = {
    "a":a,
    "ac":ac,
    "acmc":acmc,
    "mc":mc,
    "mstar":mstar,
    "m":m,
    "c":c,
    "cycb":c
}
for var, val in initial_conditions.items():
    smod.ChangeInitialSpeciesCopyNumber(var, val)

smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", species_selection=['ac','nuk', 'cycb', 'mstar','c','acmc'])

res_wo_apca = list()
for n in np.arange(1, traj_nr+1):
    smod.GetTrajectoryData(n)
    res_wo_apca.append(smod.data_stochsim.getSpecies())

### save res trajectories
f=open('./Results_trajectories/res_apca_wo.pickle', 'w')
pickle.dump(res_wo_apca, f)
f.close()

switches_apca_nuk = list()

for i in res_wo_apca:
    switches_apca_nuk.append(switch_nuk(i,1,threshold))
    
slippage_apca=list()
exit_apca=list()


for i in switches_apca_nuk:
    if i[1]<10:
        if i[1]>=1.0: #since the saddle node is lower here we consider slippage the cells which switch with nuk=1
            slippage_apca.append(i[0])
        else:
            exit_apca.append(i[0])

        

x_slp_apca,y_slp_apca = ecdf(slippage_apca)
x_exit_apca,y_exit_apca = ecdf(exit_apca)


### CREATE ECDF SLIPPAGE VS EXIT, WT VS APC-A

fig, ax = plt.subplots()

ax.step(x=(x_slp_wt-wo_time)/rescale, y=y_slp_wt,color='deeppink', label='Slippage WT');

ax.step(x=(x_exit_wt-wo_time)/rescale, y=y_exit_wt, color='g', label='Mitotic exit, WT')

ax.step(x=(x_slp_apca-wo_time)/rescale, y=y_slp_apca,color='deeppink',linestyle='dashed', label='Slippage APC-A');

ax.step(x=(x_exit_apca-wo_time)/rescale, y=y_exit_apca, color='g',linestyle='dashed', label='Mitotic exit, APC-A')

ax.axvline(x=0, color='k') # mark washout time

ax.set_xlabel('Time after washout (min)')
ax.set_ylabel('Fraction of adapted cells')
ax.set_xlim([-10,240])
ax.set_xticks(range(0,240,60))
ax.set_xticks(range(0,240,30), minor=True)
ax.set_yticks(np.arange(0,1.1,0.125), minor=True)
ax.set_yticks(np.arange(0,1.1,0.25), minor=False)
ax.grid(which='both', alpha=0.25)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2),fancybox=False, shadow=False, ncol=4)

fig.savefig('./figures/washout/wt_apca_wo_ECDF.png')
plt.close(fig)

### Count slippage and exit

without_early_adapt_apca = len(slippage_apca)+len(exit_apca)
slippage_percentage_apca = float(len(slippage_apca))/float(without_early_adapt_apca)*100
exit_percentage_apca = float(len(exit_apca))/float(without_early_adapt_apca)*100

# set back the wild-type AC association rate 
smod.ChangeParameter("kassac",kassac)
# add to table
earlies_apca = traj_nr-without_early_adapt_apca
row = {'Condition': 'APC-A mutant', 'Slippage (%)': slippage_percentage_apca, 'Exit (%)': exit_percentage_apca, 'Nr of adapters before washout': earlies_apca}
df_washout = df_washout.append(row, ignore_index=True)

######################################################################################################################

### GAL MAD2 WASHOUT

## Set initial conditions, which simulate Mad2 overexpression
c = 19
mstar = 162
mc = 57
ac = 18
acmc = 30
m = 1
a = 101
cycb = 623

initial_conditions = {
    "a":a,
    "ac":ac,
    "acmc":acmc,
    "mc":mc,
    "mstar":mstar,
    "m":m,
    "c":c,
    "cycb":c
}

for var, val in initial_conditions.items():
    smod.ChangeInitialSpeciesCopyNumber(var, val)

# do stochastic simulation
smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", species_selection=['ac','nuk', 'cycb', 'mstar','c','acmc'])
# save the results
res_wo_gm2 = list()
for n in np.arange(1, traj_nr+1):
    smod.GetTrajectoryData(n)
    res_wo_gm2.append(smod.data_stochsim.getSpecies())

### save res trajectories
f=open('./Results_trajectories/res_gm2_wo.pickle', 'w')
pickle.dump(res_wo_gm2, f)
f.close()

switches_gm2=list()
for i in res_wo_gm2:
    switches_gm2.append(switch_nuk(i,1,threshold))

fig, ax=plt.subplots()
x,y = ecdf(switches_gm2)
ax.step(x=(x-wo_time)/rescale, y=y,color='darkorange', label='GAL-MAD2 mutant')
ax.axvline(x=0, color='k', label='washout')
ax.set_xlim([-10,240])
ax.set_xticks(range(0,240,60))
ax.set_xticks(range(0,240,30), minor=True)
ax.set_yticks(np.arange(0,1.1,0.125), minor=True)
ax.set_yticks(np.arange(0,1.1,0.25), minor=False)
ax.set_ylabel("Fraction of adapted cells")
ax.set_xlabel("Time after washout (min)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend()
fig.savefig("./figures/washout/gm2_ecdf.svg")
