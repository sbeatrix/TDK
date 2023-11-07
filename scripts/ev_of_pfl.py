## script to simulate the experiment which provided evidence of a positive feedback loop in the system

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

smod.Model('model_washout.psc', dir='./model_files/')

traj_nr = 50
time = 1000
wo_time = 240
rescale = 3.33
# simulate wild-type cells
smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", species_selection=['mstar','nuk', 'cycb'])
res_wt = list()
for n in np.arange(1, traj_nr+1):
    smod.GetTrajectoryData(n)
    res_wt.append(smod.data_stochsim.getSpecies())

# simulate CDC20D YPD cells
smod.ChangeParameter("ksyn_c20", 0.001) # decrease Cdc20 synthesis to a small value

smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", species_selection=['mstar','nuk', 'cycb'])
res_c20ypd = list()
for n in np.arange(1, traj_nr+1):
    smod.GetTrajectoryData(n)
    res_c20ypd.append(smod.data_stochsim.getSpecies())

wt_label = False
ypd_label = False

handles = []
labels = []

fig, ax = plt.subplots()
for i in res_wt:
    line, = ax.plot((i[:,0]-wo_time)/rescale, i[:,1],'darkblue', alpha=0.75)
    ax.set_xlabel('Time after washout (min)')
    ax.set_ylabel('Active Mad2 (mol./cell)')
    if not wt_label:
        handles.append(line)
        labels.append('Wild-type YPD')
        wt_label = True
for i in res_c20ypd:
    line, = ax.plot((i[:,0]-wo_time)/rescale, i[:,1],'coral', alpha=0.75)
    if not ypd_label:
        handles.append(line)
        labels.append('Cdc20D YPD')
        ypd_label = True
# mark washout time
handles.append(plt.axvline(x=0, color='grey', linestyle='--'))
labels.append('Washout')
ax.set_xticks(range(0,240,60))
ax.set_xticks(range(0,240,30), minor=True)
ax.legend(handles, labels)
fig.savefig("./figures/ev_of_pfl.svg")
plt.close(fig)
