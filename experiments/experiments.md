## Simulation experiments

We planned several simulation experiments using our existing simulation models to see whether phylogenetic tree shape varies in response to certain parameters (e.g. environmental filtering, niche conservatism, dispersal, etc.) in the same way across simulation models. The full set of potential experiments is described [here](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1549290990). The idea is that each model would be run at whatever were deemed "low", "medium" and "high" settings for each of those parameters, with other parameters fixed at "average" values. 

If you are the author of one of the models we have been comparing, please provide simulation output (trees only for now) and a parameters file that provides a key to the complete set of parameter settings for each model run according to the procedures listed below.  

### Choose the relevant simulations ###  

We will analyze the following "experiments":

1. **Strength of environmental filtering** - "env"  
2. **Strength of niche conservatism** - "nic"  
3. **Strength of competition** - "com"  
4. **Strength of dispersal** - "dis"  
5. **Mutation/speciation rate** - "mut"  
6. **Time** - "tim"  

If your simulation model has a parameter related to one of the above experiments, choose Low, Medium, and High values for that parameter while holding other parameters constant at "average" values. Low, Medium and High levels are relative designations that span what is qualitatively deemed a relevant range of parameter space for a given model. This will clearly be a judgment call, especially if there are strong interactions between the parameter of interest and other model parameters.  

For the "time" experiment, what Low, Medium, and High correspond to depends on whether the simulation model is equilibrial (choose trees at time slices where t = 0.5\*t_e, t_e, and 1.5\*t_e, where t_e is the number of timesteps until richness reaches 0.95\*S_equilibrium) or non-equilibrial (break time series into 3 equal pieces, and extract time slices at t = 0.33\*t_final, 0.67\*t_final, and t_final).  



