## Simulation experiments

We planned several simulation experiments using our existing simulation models to see whether phylogenetic tree shape varies in response to certain parameters (e.g. environmental filtering, niche conservatism, dispersal) in the same way across simulation models. The full set of potential experiments is described here, and includes testing the effect of environmental filtering, niche conservatism, dispersal, speciation/mutation rate and time. The idea is that each model would be run at whatever were deemed "low", "medium" and "high" settings for each of those parameters, with other parameters fixed at "average" values. See [this spreadsheet](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1549290990) for the list of planned experiments and the set of models that can be compared for each one.

If you are the author of one of the models we have been comparing, please provide simulation output (trees only for now) and a parameters file that provides a key to the complete set of parameter settings for each model run according to the procedures listed below.  

### 1. Upload files ###  

Upload simulated tree files to the "experiments" folder of this github repo.  
Please use the following filenaming convention:  

mm_eee_L_rr.tre  

Where  
**mm**  is whatever model abbreviation we have been using for your model(e.g., "etienne", "hs", "ga", "pontarp", "oh", etc. -- see this [link](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=0))  
**eee**  is a 3-letter abbreviation for the experiment (e.g. "env", "nic", "dis", "mut", "tim"; see this [link](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1549290990))  
**L** is the 1-letter treatment level abbreviation ("L", "M", or "H" for low, medium, and high)  
**rr** is the simID/replicate number, a unique integer within the set of simulation output for model mm  

Low, Medium, and High levels are relative designations that span what is qualitatively deemed a relevant range of parameter space for a given model. This will clearly be a judgment call, especially if there are strong interactions between the parameter of interest and other model parameters.  
For the "time" experiment, what Low, Medium, and High correspond to depends on whether the simulation model is equilibrial (choose trees at time slices where t = 0.5\*t_e, t_e, and 1.5\*t_e, where t_e is the number of timesteps until richness reaches 0.95\*S_equilibrium) or non-equilibrial (break time series into 3 equal pieces, and extract time slices at t = 0.33\*t_final, 0.67\*t_final, and t_final).  

Upload trees from 10 replicates per experimental treatment. E.g., for the "strength of niche conservatism" experiment, I would upload 

hs_nic_L_1.tre  
hs_nic_L_2.tre  
hs_nic_L_3.tre  
hs_nic_L_4.tre  
hs_nic_L_5.tre  
hs_nic_L_6.tre  
hs_nic_L_7.tre  
hs_nic_L_8.tre  
hs_nic_L_9.tre  
hs_nic_L_10.tre  
hs_nic_M_1.tre  
hs_nic_M_2.tre  
hs_nic_M_3.tre  
hs_nic_M_4.tre  
hs_nic_M_5.tre  
hs_nic_M_6.tre  
hs_nic_M_7.tre  
hs_nic_M_8.tre  
hs_nic_M_9.tre  
hs_nic_M_10.tre  
hs_nic_H_1.tre  
hs_nic_H_2.tre  
hs_nic_H_3.tre  
hs_nic_H_4.tre  
hs_nic_H_5.tre  
hs_nic_H_6.tre  
hs_nic_H_7.tre  
hs_nic_H_8.tre  
hs_nic_H_9.tre  
hs_nic_H_10.tre  

Note: The replicate number need not be 1:10, repeating for each level, but could simply be unique, unrepeated simulation IDs.  

### 2. Upload parameter file ###  

Upload a parameter file in *.csv format titled "mm_experiment_parameters.csv" where **mm** is your model's abbreviation for this project.  
This file will be similar to the parameter file you may have already uploaded to the [parameters folder](https://github.com/sELDIG/SimulationStudy/tree/master/parameters).  

The first column should be called 'filename', in which the name of the simulation experiment tree is provided (e.g. 'hs_nic_L_1'; no need for the ".tre" suffix). The rest of the columns will be for the parameters specific to YOUR model, and their values.  

