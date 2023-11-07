### Script to generate trajectories for Cdc20 overexpression

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

traj_nr = 100 # simulate 100 trajectories corresponding to 100 single cells
time = 2000 # simulate for 2000 time units
rescale = 3.33 #1 time unit will be considered ~ 20 seconds

ksyn_c20 = [37, 38, 39, 40]

for k in ksyn_c20:
    smod.ChangeParameter("ksyn_c20", int(k))
 
    res_wt = list()
    for n in np.arange(1, traj_nr+1):
        smod.DoStochSim(method="direct",end=time, mode="time", species_selection=['ac','nuk'])

        # smod.GetTrajectoryData(n)
        res_wt.append(smod.data_stochsim.getSpecies())

    f=open('/storage/home/bstier/Results_trajectories/res_ksyn1030'+str(int(k))+'.pickle', 'w')
    pickle.dump(res_wt, f)
    f.close()
