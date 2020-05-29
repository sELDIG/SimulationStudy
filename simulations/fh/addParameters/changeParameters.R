pars = read.csv("fh_parameters-old.csv")

pars$env = ifelse(pars$environment == 0, "Low", "High")
pars$nic = "Medium"
pars$com = ifelse(pars$density == 0, "Low", "High")
pars$dis = ifelse(pars$dispersal == 0, "Low", "High")
pars$mut = ifelse(pars$speciationRate == 2, "Low", "High")
pars$tim = "Medium"

pars$tim = ifelse(pars$environment == 0, "Low", "High")

env,tim,com,nic,dis,mut

'env', 'nic', 'dis', 'com', 'mut', and 'tim',

# We will analyze the following "experiments":
#   
# Strength of environmental filtering - "env"
# Strength of niche conservatism - "nic"
# Strength of competition - "com"
# Strength of dispersal - "dis"
# Mutation/speciation rate - "mut"
# Time - "tim"

# If your simulation model has a parameter related to one of the above experiments, choose Low, Medium, and High values for that parameter while holding other parameters constant at "average" values. Low, Medium and High levels are relative designations that span what is qualitatively deemed a relevant range of parameter space for a given model. This will clearly be a judgment call, especially if there are strong interactions between the parameter of interest and other model parameters.
# 
# For the "time" experiment, what Low, Medium, and High correspond to depends on whether the simulation model is equilibrial (choose trees at time slices where t = 0.5*t_e, t_e, and 1.5*t_e, where t_e is the number of timesteps until richness reaches 0.95*S_equilibrium) or non-equilibrial (break time series into 3 equal pieces, and extract time slices at t = 0.33*t_final, 0.67*t_final, and t_final).



# Questions - absent treated as low?
# speciation rate and speciation type (fission, point mutation)
# "average" 
# niche conservatism != mutation?
