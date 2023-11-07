## model_files/
This folder contains all the files required for running the provided scripts. 
### model_det_ext.xml
this is an xml file producing the deterministic model. It is used by parameter_ensemble.py and parameter_sensitivity.py scripts. It is converted to SBML file.

### model_arrest.psc
this is a .psc file required by the stochastic package.  It is used for simulating cells under constant nocodazole arrest.

### model_washout.psc
this is a .psc file required by the stochastic package.  It is used for simulating cells in condition where the drug is washed out

## scripts/
This folder contains the scripts I created to fit the model to experimental data, for studying the behavior of the model and to perform stochastic simulations to reproduce experiments
### parameter_ensemble.py

This script uses the deterministic model, introduces the constraints we defined based on measurements and experimental data. First it looks for the local minimum, then starting from the optimized parameter set it build a parameter ensemble which consists of ~100 individual parameter sets which fulfill the given constraints.

### figures_arrest.py

This script simulates the constant presence of the drug, and creates the empirical cumulative distribution function of adaptation times. Wild-type, APC-A mutant and GAL-Mad2 mutant cells are simulated and the result plots are generated for all three conditions.

### figures_washout.py

In this script the effect of washing out the drug is simulated. It creates ECDFs of adaptaiton times, and plots APC/CCdc20 trajectories vs time. Wild-type, APC-A mutant and GAL-Mad2 mutant cells are simulated and the result plots are generated for all three conditions.

### parameter_sensitivity.py

This is a script to test the model robustness: starting from an initial parameter set each parameter value is disturbed in a given interval. For each changed parameter a stochastic simulation is carried out. The script tests whether the system kept bistability, if so how many of the cells adapt to the checkpoint under constant drug arrest.
Additionally I also kept track of how the steady state levels (species concentrations) change and whether the fluxes of background degradations was smaller than the fluxes of main degradations. The results are saved to a csv file

