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
Upload simulation output to the relevant folders within this repo:

1) Place phylogenies in the 'trees' folder in Newick format, labeled mm_xxxxx.tre where mm is a two-letter code identifying the model that generated it and 
xxxxx = simulation ID (as many digits as necessary) that links to the set of parameters of that simulation.

2) Place trait information in the 'traits' folder (if applicable). Files should be labeled mm_xxxxx_traits.csv and include a species column (species id's 
should match the tip labels in the tree file), and one column for each trait (environmental traits first, resource/competition traits next).

3) Place site x species info for each simulation in the 'sitexspecies' folder (if applicable). Files should be labeled mm_xxxxx_sites.csv, and include a 
site (or grid cell ID) column and a species column (i.e. long format), followed by an abundance column if applicable.

## Exploring trees in PCA space
To explore trees in PCA space, either across all models or for just a single model, first clone this repository to your local machine.

Open the repo and click on the simulationstudy.Rproj file which will open the project in RStudio. The file `treeAnalysis.r` provides a guideline for how you can explore the data yourself.

You can view the knitted Rmarkdown file by opening `simulation_comparison.Rmd` and click the `Knit` button.
