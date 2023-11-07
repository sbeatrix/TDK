### Import deterministic xml model
# Define constraints for the model
# Optimization L-M algorithm
# Parameter ensemble

import matplotlib
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *
import SCfunctions as sc
import libsbml as lib
import numpy as np

plt.rcParams['figure.dpi'] = 100

net_wt = IO.from_SBML_file("./model_files/model_det_ext.xml","net_wt")

try:
    net.compile()
except:
    net.compile()
    
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

## Initial conditions
c = 25
mstar = 106
mc = 38
ac = 30
acmc = 30
m = 1
a = 80
cycb = 430

## Total amounts
Ctot=150
Mtot=175
Atot=150

net.set_var_ics({"ksyn_c20": ksyn_c20, "ksyn_cy": ksyn_cy, "kassmc": kassmc, "kdissmc":kdissmc,
                          "kassac": kassac,"kdissac": kdissac, "kassacmc":kassacmc,
                          "kdissacmc": kdissacmc,"kdegbg": kdegbg, "kdeg":kdeg,
                          "kdegbg_cy": kdegbg_cy, "kdeg_cy": kdeg_cy, "kact":kact,
                          "kinact": kinact, "Atot":Atot, "Mtot":Mtot, "J1":J1, "k_nuk":k_nuk, "nuk":nuk, 
                          "AC":ac, "Mstar":mstar, "MC":mc, "AC":ac, "ACMC":acmc, "M":m, "A":a, 
                          "CycB":cycb, "J":J, "Ctot":Ctot})

net_wt.set_var_optimizable("Mtot",False)
net_wt.set_var_optimizable("Atot",False)
net_wt.set_var_optimizable("nuk",False)

data = {
    "net_wt": {
               "MC": {10: (36,1),
                    1000:(36,1)},
                "ACMC":{10: (30,0.1),
                        1000: (30,0.1)},
                "A":{10: (60,10),
                     1000: (60,10)},
                "CycB":{10:(400,10),
                        1000:(400,10)},
                "Ctot":{10:(150,5),
                     1000:(150,5)},
                "flux_ratio_cycb":{10:(10,0.1),
                                   1000:(10,0.1)},
                "flux_ratio_c20":{10:(10,0.1),
                                   1000:(10,0.1)},
    }
        
    }


exp_concentrations = Experiment(name = "concentrations", data = data)
exp_concentrations.set_fixed_sf({"MC":1, "ACMC":1, "A":1,"CycB":1,"Ctot":1, "flux_ratio_cycb":1,"flux_ratio_c20":1}) 

## Defining initial conditions corresponding to SAC off and SAC on

# #lower SS high nuk
var_ics_1 = KeyedList([
    ("C", 0),
    ("MC", 0),
    ("AC", 0),
    ("ACMC", 0),
    ("CycB", 0),
    ("Mstar", 175),
    ("nuk", 32),
])

# upper SS high nuk
var_ics_2 = KeyedList([
    ("C", 0),
    ("MC", 0),
    ("AC", 150),
    ("ACMC", 0),
    ("CycB", 0),
    ("Mstar", 0),
    ("nuk", 32),
])

# lower SS low nuk
var_ics_3 = KeyedList([
    ("C", 0),
    ("MC", 0),
    ("AC", 0),
    ("ACMC", 0),
    ("CycB", 0),
    ("Mstar", 175),
    ("nuk", 1),
])

# upper SS low nuk
var_ics_4 = KeyedList([
    ("C", 0),
    ("MC", 0),
    ("AC", 150),
    ("ACMC", 0),
    ("CycB", 0),
    ("Mstar", 0),
    ("nuk", 1),
])


net_1 = net_wt.copy(new_id = "net_1")
net_1.set_var_ics(var_ics_1)
traj_1 = net_1.integrate([0, 2000])

net_2 = net_wt.copy(new_id = "net_2")
net_2.set_var_ics(var_ics_2)
traj_2 = net_2.integrate([0, 2000])

net_3 = net_wt.copy(new_id = "net_3")
net_3.set_var_ics(var_ics_3)
traj_3 = net_3.integrate([0, 2000])

net_4 = net_wt.copy(new_id = "net_4")
net_4.set_var_ics(var_ics_4)
traj_4 = net_4.integrate([0, 2000])


## Defining the data that impose to stable steady states
data = {
    "net_1": {
        "AC": {
            500: (10, 10),
            1500: (10, 10)
        },
    },
    "net_2": {
        "AC": {
            500: (100, 20),
            1500: (100, 20)
        }
    },
     "net_3": {
        "AC": {
            500: (10, 10),
            1500: (10, 10)
        }
    },
    "net_4": {
        "AC": {
            500: (100, 20),
            1500: (100, 20)
        } 
      
      
exp_bistability = Experiment(name='bistability',data = data)
exp_bistability.set_fixed_sf({"AC": 1.0})
      


m = Model([exp_concentrations, exp_bistability], [net_wt, net_1,net_2,net_3,net_4])
params = m.get_params()
for var, val in params.items():
    res = Residuals.PriorInLog(var+'_prior', var, np.log(val), np.log(np.sqrt(10)))
    m.AddResidual(res)
      
      
params = KeyedList([(key, params.getByKey(key)) for key in m.params.keys()])
params_opt = Optimization.fmin_lm_log_params(m, params, maxiter=1000, disp=True)
      
## run ensemble analysis
Network.full_speed()
j = m.jacobian_log_params_sens(np.log(params_opt))
jtj = np.dot(np.transpose(j), j)
ens, gs, r = Ensembles.ensemble_log_params(m, params_opt, jtj, steps=100000, temperature=1)
      
ens_pruned=ens[::1000]
      
## example of plotting parameter distributions
n_rows = 2
n_cols = 7

fig, axes = plt.subplots(n_rows, n_cols, figsize=(11,3), sharey = True)
logbins = np.logspace(np.log10(1e-2), np.log10(1e2), 30)
for i, ax in enumerate(axes.flat):
    if i==len(params_opt):
        break
    var = params_opt.keys()[i]
    val = np.array([e[i] for e in ens_pruned])/params.getByKey(var)
    ax.set_xscale('log')
    ax.hist(val, color = "white", edgecolor = "black", bins = logbins);
    ax.set_title(var)
plt.tight_layout()
fig.savefig("./figures/ensemble.svg")
