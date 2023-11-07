### Script to generate the following plots:

# APC/Cdc20 trajectories vs time, nuk constant (nocodazole arrest)
# cumulative distribution: fraction of adapted cells
   # Will-type cells
   # APC-A mutant cells
   # Gal-Mad2 mutant cells
import matplotlib
matplotlib.use('Agg')
import os
import pickle
import stochpy
smod = stochpy.SSA()
smod.DoStochSim()
smod.PlotSpeciesTimeSeries

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.pyplot import cm

plt.rcParams['figure.dpi'] = 150

smod.Model('model_arrest.psc', dir='./model_files/')

## Functions

# function to define the number of unattached kineotchores and the time point the SAC on - SAC off switch happens. We follow AC transition from the lower 
# (SAC on) to upper steady state (SAC off)
# input:
## traj - stochastic simulation result trajectories
## species_col - the column index of the species from species selection, from the steady states of which we calculate the threshold
# output[0] - timepoint of adaptation
# output[1] - number of unattached kinetochores at the time of adaptation
def switch_nuk(traj, species_col, threshold):
    t_switch_up = np.argmax(traj[:,species_col]>threshold)
    t_switch_down = np.argmax(traj[:,species_col]<threshold)
    if t_switch_up==0:
        t_switch_up=np.inf
        return[np.inf, traj[t_switch_down,0], traj[t_switch_down,2]]
    return [traj[t_switch_up,0], traj[t_switch_up,2]]

# function to compute empirical cumulative distribution function
def ecdf(data):
    x = np.sort(data)
    n = x.size
    y = np.arange(1, float(n)+1) / float(n)
    return(x,y)

### Simulation settings

traj_nr = 100 # simulate 100 trajectories corresponding to 100 single cells
time = 2000 # simulate for 2000 time units
rescale = 3.33 #1 time unit will be considered ~ 20 seconds

# Do stochastic simulation for wild-type
smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", species_selection=['ac','nuk', 'cycb', 'mstar','c','acmc'])
# Save the results
res_wt = list()
for n in np.arange(1, traj_nr+1):
    smod.GetTrajectoryData(n)
    res_wt.append(smod.data_stochsim.getSpecies())

### save res trajectories
f=open('./Results_trajectories/res_wt_arrest.pickle', 'w')
pickle.dump(res_wt, f)
f.close()
    
# APC/CCdc20 trajectories vs time, wildtype
fig, ax = plt.subplots()
for i in res_wt:
    ax.plot(i[:,0]/rescale, i[:,1],alpha=0.75)
    ax.set_xlabel('Time (min)',fontsize=15)
    ax.set_ylabel('APC/CCdc20 (molecules/cell)', fontsize=15)
    ax.set_title('Constant nocodazole arrest, wild-type')
   
fig.savefig('../figures/arrest/wt_constant_noco_apc.png')
plt.close(fig)

### GENERATE ECDF: FRACTION OF ADAPTED CELLS VS TIME

switches_wildtype = list()
threshold = 55 # calculated as the middle value between SAC_off steady state - SAC_on steady state for APC/CCdc20. Steady states go

# use the switch_nuk function outputs first argument to list the timepoints where transitions happen
for r in res_wt:
    switches_wildtype.append(switch_nuk(r,1,threshold)[0])
    
fig, ax = plt.subplots()
x,y = ecdf(switches_wildtype)
ax.step(x=x/rescale, y=y,color='r')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Fraction of adapted cells')
ax.set_title('Constant nocodazole arrest, wild-type')
fig.savefig('./figures/arrest/wt_constant_noco_ECDF.png')
plt.close(fig)
######################################################################################################################

### APC-A MUTANT 

new_kassac = 0.015 # 40% decrease of APC/C - Cdc20 binding affinity
smod.ChangeParameter("kassac",new_kassac)

### INITIAL CONDITIONS APCA MUTANT (steady states from integration with XPPAUT)
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
# change the initial conditions
for var, val in initial_conditions.items():
    smod.ChangeInitialSpeciesCopyNumber(var, val)

# do stochastic simulation for APC-A mutant cells
# time and traj_nr remain the same
smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", species_selection=['ac','nuk', 'cycb', 'mstar','c','acmc'])

# get the trajectories 
res_apca = list()
for n in np.arange(1, traj_nr+1):
    smod.GetTrajectoryData(n)
    res_apca.append(smod.data_stochsim.getSpecies())

### save res trajectories to file
f=open('../Results_trajectories/res_apca_arrest.pickle', 'w')
pickle.dump(res_apca, f)
f.close()
    
# APC/CCdc20 trajectories vs time, APC-A mutant
fig, ax= plt.subplots()
for i in res_apca:
    ax.plot(i[:,0]/rescale, i[:,1],alpha=0.75)
    ax.set_xlabel('Time (min)',fontsize=15)
    ax.set_ylabel('APC/CCdc20 (molecules/cell)', fontsize=15)
    ax.set_title('Constant nocodazole arrest, APC-A mutant')

fig.savefig('../figures/arrest/apca_constant_noco_apc.png')
plt.close(fig)

### GENERATE ECDF: FRACTION OF ADAPTED CELLS VS TIME

switches_apca = list()
for r in res_apca:
    switches_apca.append(switch_nuk(r,1,threshold)[0])
    

fig, ax = plt.subplots()
x,y = ecdf(switches_apca)
ax.step(x=x/rescale, y=y,color='r')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Fraction of adapted cells')
ax.set_title('Constant nocodazole arrest, APC-A mutant')
fig.savefig('./figures/arrest/apca_constant_noco_ECDF.png')
plt.close(fig)

# set back the original value
smod.ChangeParameter("kassac",kassac)

######################################################################################################################


### GAL MAD2 mutant

## Set initial conditions, which give MAD2 x 2, given by steady states from integration with XPPAUT   
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
# set initial conditions
for var, val in initial_conditions.items():
    smod.ChangeInitialSpeciesCopyNumber(var, val)

# do stochastic simulation    
smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", species_selection=['ac','nuk', 'cycb', 'mstar','c','acmc'])

# get trajectories
res_gm2 = list()
for n in np.arange(1, traj_nr+1):
    smod.GetTrajectoryData(n)
    res_gm2.append(smod.data_stochsim.getSpecies())

### save res trajectories to file
f=open('./Results_trajectories/res_gm2_arrest.pickle', 'w')
pickle.dump(res_gm2, f)
f.close()
fig, ax = plt.subplots()    
# APC/CCdc20 trajectories vs time, GAL-Mad2 mutant
for i in res_gm2:
    ax.plot(i[:,0]/rescale, i[:,1],alpha=0.75)
    ax.set_xlabel('Time (min)',fontsize=15)
    ax.set_ylabel('APC/CCdc20 (molecules/cell)', fontsize=15)
    ax.set_title('Constant nocodazole arrest, GAL-Mad2 mutant')

fig.savefig('./figures/arrest/gm2_constant_noco_apc.png')
plt.close(fig)

### GENERATE ECDF: FRACTION OF ADAPTED CELLS VS TIME

switches_gm2 = list()
    
for r in res_gm2:
    switches_gm2.append(switch_nuk(r,1,threshold)[0])
   

fig, ax = plt.subplots()
x,y = ecdf(switches_gm2)
ax.step(x=x/rescale, y=y,color='r')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Fraction of adapted cells')
ax.set_title('Constant nocodazole arrest, GAL-Mad2 mutant')
fig.savefig('./figures/arrest/gm2_constant_noco_ECDF.png')

plt.close(fig)
