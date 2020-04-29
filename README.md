# Working title: 
Inferring modes of diversification from simulation models: identifying robust patterns and reliable methods

## Interested group members: (add your name here!)
Allen*, Mikael, Juliano, Florian, Roland, Thorsten, Susanne, Shan, Thiago, Odile, Loic, Oskar

## Project Links
[GoogleDoc](https://docs.google.com/document/d/1F9rXWp_DAleZarrXXYzAgbGqgjtE2dTUe7VqpxMdBW4/edit)  
[Simulation Spreadsheet](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit?usp=sharing)

## Summary
The complete set of processes governing the diversification of clades over time and space is unknown, yet many simulation models exist that try to codify subsets of possible processes and examine the consequences for various spatial and phylogenetic patterns. In some cases it is unknown how much the output from a given simulation model is sensitive to the specific implementation decisions used to model the processes of interest. Comparing simulation models that use very different implementations of similar processes can shed light on which patterns are most robust and which processes are most reliably inferred.

## Questions
1) Which empirical phylogenetic metrics are most variable and which are robust to the different assumptions made across distinct simulation modeling frameworks? 
2) How well do common analytical frameworks (RPANDA, BAMM, DDD, etc) capture actual diversification rates of data not simulated under the model generating processes?
3) What are the primary axes of variation of existing diversification simulation models, and how do those axes map onto distinct phylogenetic and spatial patterns?

## Steps
Upload simulation output to the relevant folders within this repo. Presumably you will upload output from a range of parameter space over which you ran your model. See separate instructions on the [Experiments Page](https://github.com/sELDIG/SimulationStudy/blob/master/experiments/experiments.md) for some of the variation in parameter space we are especially interested in. If simulation output is highly stochastic, you might consider providing output from replicate simulations with identical parameter settings, but in most cases 10-20 replicates should be plenty.  

To make your simulation output visualizable in our [shiny app](https://hurlbertlab.shinyapps.io/simulationstudy/), follow the instructions below. This may include areas of parameter space not highlighted in the Experiments. Simulation output that will ALSO be part of the Experiments needs to be uploaded separately according to [those instructions](https://github.com/sELDIG/SimulationStudy/blob/master/experiments/experiments.md). (Sorry for the redundancy in effort; this is a locally adaptive maximum for me at the moment!)  


1) Place phylogenies in the 'trees' folder in Newick format, labeled mm_xxxxx.tre where mm is a two-letter code identifying the model that generated it and 
xxxxx = simulation ID (as many digits as necessary) that links to the set of parameters of that simulation.

2) (IGNORE FOR NOW) Place trait information in the 'traits' folder (if applicable). Files should be labeled mm_xxxxx_traits.csv and include a species column (species id's 
should match the tip labels in the tree file), and one column for each trait (environmental traits first, resource/competition traits next).

3) (IGNORE FOR NOW) Place site x species info for each simulation in the 'sitexspecies' folder (if applicable). Files should be labeled mm_xxxxx_sites.csv, and include a 
site (or grid cell ID) column and a species column (i.e. long format), followed by an abundance column if applicable.

4) Upload a parameters file, labeled mm_parameters.csv where mm is the two-letter code identifying the model that generated it. In this file, the first column will be called simID and will hold the unique simulation ID you are using to distinguish the different model runs. The other columns will be for the parameter settings of your model, so that a single row will describe the parameters used to generate a given set of output. (This will allow us to eventually explore within a given model the effect of these settings on metrics of tree shape within the [shiny app](https://hurlbertlab.shinyapps.io/simulationstudy/)).

5) Classify your simulation model according to the modeling approach and processes included on this [Google Spreadsheet](https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=2047946073). Importantly, there is a separate row for every unique simulation ID in case some of these values differ based on parameter settings, but in some cases all rows for a given model could potentially be identical. This information is for comparing between models on the [shiny app](https://hurlbertlab.shinyapps.io/simulationstudy/).

6) Contact me (Allen: hurlbert AT bio.unc.edu) to let me know if you have uploaded trees, parameter file, and model classification. I will update the shiny app so that new models will then be visualizable.
