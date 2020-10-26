# Working title: 
Inferring modes of diversification from simulation models: identifying robust patterns and reliable methods

## Interested group members: (add your name here!)
Allen*, Mikael, Juliano, Florian, Susanne, Shan, Thiago, Loïc, Oskar, Catherine, Rampal, Antonin, Pedro Neves, Liang Xu, David

## Project Links
[GoogleDoc](https://docs.google.com/document/d/1F9rXWp_DAleZarrXXYzAgbGqgjtE2dTUe7VqpxMdBW4/edit)  
[Simulation Spreadsheet](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit?usp=sharing)

## Summary
The complete set of processes governing the diversification of clades over time and space is unknown, yet many simulation models exist that try to codify subsets of possible processes and examine the consequences for various spatial and phylogenetic patterns. In some cases it is unknown how much the output from a given simulation model is sensitive to the specific implementation decisions used to model the processes of interest. Comparing simulation models that use very different implementations of similar processes can shed light on which patterns are most robust and which processes are most reliably inferred.

## Experiment  
Upload simulation output to the relevant folders within this repo (trees > uniform_sampling_experiment). For the currently prescribed methods for our simulation experiment, refer to the project [Google Doc](https://docs.google.com/document/d/1F9rXWp_DAleZarrXXYzAgbGqgjtE2dTUe7VqpxMdBW4/edit).   

## Shiny app

To make your simulation output visualizable in our [shiny app](https://hurlbertlab.shinyapps.io/simulationstudy/), follow the instructions below. I will eventually make sure all (or at least some subset of) simulated trees provided in the **Experiment** get added here.  

1) Place phylogenies in the 'trees' folder in Newick format, labeled **mm_xxxxx.tre** where **mm** is a two-letter code identifying the model that generated it and 
**xxxxx** = simID (as many digits as necessary, should be an integer) that links to the set of parameters of that simulation.

2) (IGNORE FOR NOW) Place trait information in the 'traits' folder (if applicable). Files should be labeled mm_xxxxx_traits.csv and include a species column (species id's 
should match the tip labels in the tree file), and one column for each trait (environmental traits first, resource/competition traits next).

3) (IGNORE FOR NOW) Place site x species info for each simulation in the 'sitexspecies' folder (if applicable). Files should be labeled mm_xxxxx_sites.csv, and include a 
site (or grid cell ID) column and a species column (i.e. long format), followed by an abundance column if applicable.

4) Upload a parameters file labeled **mm_parameters.csv** to the **parameters** folder of this repo, where **mm** is the abbreviated code identifying the model that generated it (e.g. 'hs', 'etienne', 'fh', etc.).   
* This file should have one row for every simulation run for which output is provided.  
* In this parameters file, the first column should be called **model** and should be filled with the model abbreviation.  
* The second column should be called **simID** and will hold the unique simulation ID (integer only please) you are using to distinguish the different model runs.  
* The next columns will be for the parameter settings of your model, thus the number of columns is expected to vary from model to model. (This will allow us to eventually explore within a given model the effect of these settings on metrics of tree shape within the [shiny app](https://hurlbertlab.shinyapps.io/simulationstudy/)).  
* The last six columns should be called **'env', 'nic', 'dis', 'com', 'mut',** and **'tim'**, and they refer to how a given tree fits into the six **simulation experiments** described [here](https://github.com/sELDIG/SimulationStudy/blob/master/experiments/experiments.md). If a given simulation output has been designated as one of the three treatment levels (Low, Medium, High) for a given experiment, then it will have a 'L', 'M', or 'H' in the column for the relevant experiment. **NOTE1:** You may choose to upload many simulations that explore parameter space that are not part of the specific experiments. In this case, leave these fields blank for those simulations. **NOTE2:** A single simulation run may be a part of multiple experiments, although this will be most likely for 'Medium' treatment levels since we are attempting to hold all other parameters at "average" levels while varying the one of interest. **NOTE3:** See example [here](https://github.com/sELDIG/SimulationStudy/blob/master/parameters/hs_parameters.csv).    

5) Classify your simulation model according to the modeling approach and processes included on this [Google Spreadsheet](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=2047946073). Importantly, there is a separate row for every unique simulation ID in case some of these values differ based on parameter settings, but in some cases all rows for a given model could potentially be identical. This information is for comparing between models on the [shiny app](https://hurlbertlab.shinyapps.io/simulationstudy/).

6) Contact me (Allen: hurlbert AT bio.unc.edu) to let me know if you have uploaded trees, parameter file, and model classification. I will update the shiny app so that new models will then be visualizable.
