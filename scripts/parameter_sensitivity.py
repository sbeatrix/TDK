# Script to test the model robustness
# starting from an initial parameter set each parameter value is disturbed in a given interval. For each changed parameter a stochastic simulation
# is carried out. The script tests whether the system kept bistability, if so how many of the cells adapt to the checkpoint under constant drug arrest.
#Additionally I also kept track of how the steady state levels (species concentrations) change and whether the fluxes of background degradations 
#was smaller than the fluxes of main degradations. The results are saved to a csv file


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import stochpy
smod = stochpy.SSA()
smod.DoStochSim()
smod.PlotSpeciesTimeSeries
import os
import numpy as np

from SloppyCell.ReactionNetworks import *
import SCfunctions as sc
import libsbml as lib

import pandas as pd

plt.rcParams['figure.dpi'] = 100

# import model for stochpy
smod.Model('model_arrest.psc', dir='./model_files/')
# import model for sloppycell
net = IO.from_SBML_file("./model_files/model_det_ext.xml","net")

## possible sloppycell installation error fix:
# try:
#     net.compile()
# except:
#     net.compile()


## starting parameter set
ksyn_c20 = 36.97731502237753
ksyn_cy = 316.94570558508684
kdegbg = 0.022521686855703566
kinact = 0.45009544921120104
J = 0.06168169109783722
kact = 0.0014200271790012738
k_nuk = 78
J1 = 0.5
kassmc = 0.021973634967776146
kdissmc = 0.6256358653434726
kassac = 0.023088755284985777
kassacmc = 0.04477093007618703
kdissac = 0.664021730464077
kdissacmc = 0.43441487879976926
kdeg = 1.117610105994406
kdeg_cy = 0.023952471280495782
kdegbg_cy = 0.0677970961252346

nuk = 32

parameters = {
    "ksyn_c20": ksyn_c20,
    "ksyn_cy": ksyn_cy,
    "kdegbg": kdegbg,
    "kinact": kinact,
    "J": J,
    "kact": kact,
    "k_nuk": k_nuk,
    "J1": J1,
    "kassmc": kassmc,
    "kdissmc": kdissmc,
    "kassac": kassac,
    "kassacmc": kassacmc,
    "kdissac": kdissac,
    "kdissacmc": kdissacmc,
    "kdeg": kdeg,
    "kdeg_cy": kdeg_cy,
    "kdegbg_cy": kdegbg_cy
}

c = 25
mstar = 106
mc = 38
ac = 30
acmc = 30
m = 1
a = 80
cycb = 430

#JUST FOR SC
Ctot=150
Mtot=175
Atot=150

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

### set starting parameter values and initial conditions stochpy
for var, val in parameters.items():
    smod.ChangeParameter(var, val)
    
for var, val in initial_conditions.items():
    smod.ChangeInitialSpeciesCopyNumber(var, val)

### set ICs for sloppycell
net.set_var_ics({"ksyn_c20": ksyn_c20, "ksyn_cy": ksyn_cy, "kassmc": kassmc, "kdissmc":kdissmc,
                          "kassac": kassac,"kdissac": kdissac, "kassacmc":kassacmc,
                          "kdissacmc": kdissacmc,"kdegbg": kdegbg, "kdeg":kdeg,
                          "kdegbg_cy": kdegbg_cy, "kdeg_cy": kdeg_cy, "kact":kact,
                          "kinact": kinact, "Atot":Atot, "Mtot":Mtot, "J1":J1, "k_nuk":k_nuk, "nuk":nuk, 
                          "AC":ac, "Mstar":mstar, "MC":mc, "AC":ac, "ACMC":acmc, "M":m, "A":a, 
                          "CycB":cycb, "J":J, "Ctot":Ctot})

## Define functions for transitions
def switch(traj, species_col, threshold):
    t_switch_up = np.argmax(traj[:,species_col]>threshold)
    if t_switch_up==0:
        t_switch_up=np.inf
        return np.inf
    return traj[t_switch_up,0]
# Define function for generating empirical cumulative distributions
def ecdf(data):
    x = np.sort(data)
    n = x.size
    y = np.arange(1, float(n)+1) / float(n)
    return(x,y)

# value by which the iterated interval for each parameter is determined [parval/scale -> parval*scale]
scale = 10 ## 

# nr of values evaluated in the interval
nr = 10 ## 

# nr of trajectories in one simulation
traj_nr = 50

# time of one simulation
time = 1000


## create dataframe for results
df_ind = ['par_changed', 'new_val', 'input_change%',  'bistability', 'nr_adapted_cells','ECDF_filename','flux_cy', 'flux_c20', 'A', 'MC', 'CycB', 'ACMC']

results_rates = pd.DataFrame(columns=df_ind)


### loop through parameters
for param, val in parameters.items():
    ## create values on which the actual parameter will iterate through
    parvals = np.logspace(np.log10(val/scale), np.log10(val*scale), nr) 
    folder_name = param
    filename=0
    for i in parvals:
        i_formatted = "{:.4f}".format(i)
        filename+=1
        smod.ChangeParameter(param, i)
        net.set_var_ic(param, i)
        ## get SAC on steady states with the new parameter (nuk high) for ICs
        net.set_var_ic("nuk", 32)
        net.integrate([0,2000])

        # keep track of the ratio of flux of degradation/ bg degradation for cdc20 and cycb
        flux_cy = net.get_var_val("flux_ratio_cycb")
        flux_c20 = net.get_var_val("flux_ratio_c20")
    
        new_ic = {
            "a":int(net.get_var_val("A")),
            "ac":int(net.get_var_val("AC")),
            "acmc":int(net.get_var_val("ACMC")),
            "mc":int(net.get_var_val("MC")),
            "mstar":int(net.get_var_val("Mstar")),
            "m":int(net.get_var_val("M")),
            "c":int(net.get_var_val("C")),
            "cycb":int(net.get_var_val("CycB"))
        }
        
        for ic_var, ic_val in new_ic.items():
            smod.ChangeInitialSpeciesCopyNumber(ic_var, ic_val)
        # get SAC off steady states (nuk low) for AC (for transition def)    
        net.set_var_ic("nuk", 0.1)
        net.integrate([0,2000])
        AC_SAC_off = int(net.get_var_val("AC"))
        AC_SAC_on = new_ic.get("ac")
        smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", 
                        species_selection=['ac','nuk', 'cycb', 'mstar'])
        res = list()
        for n in np.arange(1, traj_nr+1):
            smod.GetTrajectoryData(n)
            res.append(smod.data_stochsim.getSpecies())
        threshold = AC_SAC_on + (AC_SAC_off - AC_SAC_on)/2
                            
        switches = list()
        for r in res:
            switches.append(switch(r,1,threshold))
        switches.sort()    
        nr_adapted_cells = len([s for s in switches if s<=time])
        
        x,y = ecdf(switches)
        # plt.figure()
        fig, ax = plt.subplots()
        ax.step(x=x, y=y,color='r')
        if (AC_SAC_off-AC_SAC_on) < 20:
            ax.set_title("{} changed to {}, bistability is lost".format(param, i_formatted))
            bist = 0
            nr_adapted_cells = 'NaN'
        else:
            ax.set_title("{} changed to {}, {} cell(s) adapted".format(param, i_formatted, nr_adapted_cells))        
            bist = 1
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Fraction of adapted cells')
        
        if not os.path.exists('./figures/parameters/'+folder_name):
            os.makedirs('./figures/parameters/'+folder_name)
        
        ecdf_file = folder_name+'/res'+str(filename)+'.png'
        fig.savefig('/storage/home/bstier/figures/parameters/'+ecdf_file)
        plt.close(fig)
        input_change = val/i
        new_res = {'par_changed': param, 'new_val': i, 'input_change%': input_change, 'bistability':bist, 'nr_adapted_cells':nr_adapted_cells, 'ECDF_filename':ecdf_file, 'flux_cy':flux_cy, 'flux_c20':flux_c20, 'A':new_ic.get("a"), 'MC':new_ic.get("mc"), 'CycB':new_ic.get("cycb"), 'ACMC':new_ic.get("acmc")}
        results_rates = results_rates.append(new_res, ignore_index=True)
    smod.ChangeParameter(param, val)
    net.set_var_ic(param, val)
    

#############################################################################################################   
#############################################################################################################   

### Varying Mtotal and A total (2-fold)

multiplier = np.linspace(0.5, 2, nr)
results_totals = pd.DataFrame(columns=df_ind)
## Atotal
folder_name = 'Atot'
filename = 0
for m in multiplier:
    filename+=1
    net.set_var_ics({"Atot":Atot*m, "A":Atot*m, "AC":0, "ACMC":0, "MC":0, "M":0, "Mstar":0, "C":0, "CycB":0})
                    
    # get SAC on steady states
    
    net.set_var_ic("nuk",32)
    net.integrate([0,2000])



    new_ic = {
            "a":int(net.get_var_val("A")),
            "ac":int(net.get_var_val("AC")),
            "acmc":int(net.get_var_val("ACMC")),
            "mc":int(net.get_var_val("MC")),
            "mstar":int(net.get_var_val("Mstar")),
            "m":int(net.get_var_val("M")),
            "c":int(net.get_var_val("C")),
            "cycb":int(net.get_var_val("CycB"))
        }
    for ic_var, ic_val in new_ic.items():
            smod.ChangeInitialSpeciesCopyNumber(ic_var, ic_val)

      # get SAC off steady states (nuk low) for AC (for transition def)  
    net.set_var_ic("nuk", 0.1)
    net.integrate([0,2000])
    AC_SAC_off = int(net.get_var_val("AC"))
    AC_SAC_on = new_ic.get("ac")
    
    ## keep track of degradation flux ratios
    flux_cy = net.get_var_val("flux_ratio_cycb")
    flux_c20 = net.get_var_val("flux_ratio_c20")

    smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", 
                    species_selection=['ac','nuk', 'cycb', 'mstar'])
    res = list()

    for n in np.arange(1, traj_nr+1):
        smod.GetTrajectoryData(n)
        res.append(smod.data_stochsim.getSpecies())
    threshold = AC_SAC_on + (AC_SAC_off - AC_SAC_on)/2

    switches = list()
    for r in res:
        switches.append(switch(r,1,threshold))
    switches.sort()    
    nr_adapted_cells = len([s for s in switches if s<=time])

    x,y = ecdf(switches)

    fig, ax = plt.subplots()
    ax.step(x=x, y=y,color='r')
    if (AC_SAC_off-AC_SAC_on) < 20:
        ax.set_title("Atot changed to {}, bistability is lost".format(Atot*m))
        bist = 0
        nr_adapted_cells = 'NaN'
    else:
        ax.set_title("Atot changed to {}, {} cell(s) adapted".format(Atot*m, nr_adapted_cells))        
        bist = 1
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Fraction of adapted cells')

    if not os.path.exists('./figures/parameters/'+folder_name):
            os.makedirs('./figures/parameters/'+folder_name)
    ecdf_file = folder_name+'/res'+str(filename)+'.png'
    fig.savefig('./figures/parameters/'+ecdf_file)
    plt.close(fig)
    new_res = {'par_changed': 'Atot', 'new_val': Atot*m, 'input_change%': m, 'bistability':bist, 'nr_adapted_cells':nr_adapted_cells, 'ECDF_filename':ecdf_file, 'flux_cy':flux_cy, 'flux_c20':flux_c20, 'A':new_ic.get("a"), 'MC':new_ic.get("mc"), 'CycB':new_ic.get("cycb"), 'ACMC':new_ic.get("acmc")}
    results_totals = results_totals.append(new_res, ignore_index=True)
    
## Varying Mtotal 

multiplier = np.linspace(0.5, 2, nr)

folder_name = 'Mtot'
filename = 0
for m in multiplier:
    filename+=1
    net.set_var_ics({"Atot":Atot,"Mtot":Mtot*m, "A":0, "AC":0, "ACMC":0, "MC":0, "M":0, "Mstar":Mtot*m, "C":0, "CycB":0})
                    
    # get SAC on steady states
    
    net.set_var_ic("nuk",32)
    net.integrate([0,2000])



    new_ic = {
            "a":int(net.get_var_val("A")),
            "ac":int(net.get_var_val("AC")),
            "acmc":int(net.get_var_val("ACMC")),
            "mc":int(net.get_var_val("MC")),
            "mstar":int(net.get_var_val("Mstar")),
            "m":int(net.get_var_val("M")),
            "c":int(net.get_var_val("C")),
            "cycb":int(net.get_var_val("CycB"))
        }
    for ic_var, ic_val in new_ic.items():
            smod.ChangeInitialSpeciesCopyNumber(ic_var, ic_val)

      # get SAC off steady states (nuk low) for AC (for transition def)  

     #########
    net.set_var_ic("nuk", 0.1)
    net.integrate([0,2000])
    AC_SAC_off = int(net.get_var_val("AC"))
    AC_SAC_on = new_ic.get("ac")
     ##########   
    flux_cy = net.get_var_val("flux_ratio_cycb")
    flux_c20 = net.get_var_val("flux_ratio_c20")

    smod.DoStochSim(method="direct",end=time,trajectories=traj_nr,mode="time", 
                    species_selection=['ac','nuk', 'cycb', 'mstar'])
    res = list()

    for n in np.arange(1, traj_nr+1):
        smod.GetTrajectoryData(n)
        res.append(smod.data_stochsim.getSpecies())
    threshold = AC_SAC_on + (AC_SAC_off - AC_SAC_on)/2

    switches = list()
    for r in res:
        switches.append(switch(r,1,threshold))
    switches.sort()    
    nr_adapted_cells = len([s for s in switches if s<=time])

    x,y = ecdf(switches)
    # plt.figure()
    fig, ax = plt.subplots()
    ax.step(x=x, y=y,color='r')
    if (AC_SAC_off-AC_SAC_on) < 20:
        ax.set_title("Mtot changed to {}, bistability is lost".format(Mtot*m))
        bist = 0
        nr_adapted_cells = 'NaN'
    else:
        ax.set_title("Mtot changed to {}, {} cell(s) adapted".format(Mtot*m, nr_adapted_cells))        
        bist = 1
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Fraction of adapted cells')

    if not os.path.exists('./figures/parameters/'+folder_name):
        os.makedirs('./figures/parameters/'+folder_name)
        
    ecdf_file = folder_name+'/res'+str(filename)+'.png'
    fig.savefig('./figures/parameters/'+ecdf_file)
    plt.close(fig)
    new_res = {'par_changed': 'Mtot', 'new_val': Mtot*m, 'input_change%': m, 'bistability':bist, 'nr_adapted_cells':nr_adapted_cells, 'ECDF_filename':ecdf_file, 'flux_cy':flux_cy, 'flux_c20':flux_c20, 'A':new_ic.get("a"), 'MC':new_ic.get("mc"), 'CycB':new_ic.get("cycb"), 'ACMC':new_ic.get("acmc")}
    results_totals = results_totals.append(new_res, ignore_index=True)


### Save dataframes to csv files

results_rates.to_csv('./Results_csv/results_rates.csv', index=False)
results_totals.to_csv('./Results_csv/results_totals.csv', index=False)
