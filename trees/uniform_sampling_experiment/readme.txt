Upload trees that have been created by uniformly sampling parameter space for the sELDiG Simulation Experiment.

As before, we are only interested in the phylogenetic tree output of the simulation, stored in Newick format (with a .tre extension). 
Tree files should be saved with a consistent abbreviation for the model that we’ve been using (e.g. hs, pontarp, fh, etc.) followed by an underscore and then an integer representing the simID. 

Also to this folder, please upload a parameters file 
labeled "[model]_USE_parameters.csv" where [model] is the model abbreviation
Column 1 is named "model" and is filled with the model abbreviation
Column 2 is called "simID"
Subsequent columns correspond to the parameters of your model

If you have multiple simulation "scenarios" (e.g. conceptually distinct sets of rules based on one or a few categorical settings) which you would like to examine effectively as distinct submodels, then add a "scenarios" column to the end of your parameters file, and use an abbreviation to distinguish which scenario each simulation ID represents.

Finally, please specify which parameters of your model are relevant to the 6 processes/experiments we are examining by filling in this spreadsheet:
https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1171496897
