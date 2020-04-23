## Simulation experiments

We planned several simulation experiments to see whether phylogenetic tree shape varies in response to certain parameters (e.g. environmental filtering, niche conservatism, dispersal) in the same way across simulation models. The full set of potential experiments is described here, and includes testing the effect of environmental filtering, niche conservatism, dispersal, speciation/mutation rate and time. The idea is that each model would be run at whatever were deemed "low", "medium" and "high" settings for each of those parameters, with other parameters fixed at "average" values.  

#EXPERIMENTAL FILES#  
Upload simulated tree files to the "experiments" folder of the github repo.  
Please use the following filenaming convention:  

xx_eee_L_mm.tre  

Where
xx  is whatever model abbreviation we have been using for your model(e.g., etienne, hs, ga, pontarp, oh, etc. -- see this [link](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=0))
eee  is a 3-letter abbreviation for the experiment (e.g. "env", "nic", "dis", "mut", "tim"; see this [link](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1549290990))
L is the 1-letter treatment level abbreviation ("L", "M", or "H" for low, medium, and high)
mm is the simID/replicate number, a unique integer within the set of simulation output for model xx

Low, Medium, and High levels are relative designations that span what is qualitatively deemed a relative range of parameter space for a given model. For the "time" experiment, what Low, Medium, and High correspond to depends on whether the simulation model is equilibrial (choose trees at time slices where t = 0.5*t_e, t_e, and 1.5*t_e, where t_e is the number of timesteps until richness reaches 0.95\*S_equilibrium) or non-equilibrial (break time series into 3 equal pieces, and extract time slices at t = 0.33\*t_final, 0.67\*t_final, and t_final).

How many replicates per treatment? 10? 20?
