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
Compile phylogenies from all simulation models  
Trees in Newick format  
Filename convention?  
Suggestion: mm_xxxxx.tre  
Where mm is a two-letter code identifying the model that generated it:  
hs = Hurlbert-Stegen  
po = Pontarp  
ha = Hartig  
ra = Rangel  
ca = Cabral  
Etc  
xxxxx = simulation ID (as many digits as necessary)  
Each modeler creates table with simulation ID and relevant parameters  
E.g., something like [this]()  
Calculate empirical phylogenetic metrics (gamma, beta, Rampal analysis, clade age-richness correlations, etc)  
Need to identify complete list of metrics  
Conduct PCA (or equivalent) and identify how these simulated trees cluster (or not) with respect to:  
Simulation method (Florian, Allen, Mikael, Juliano, Thiago, etc)  
Process categories (e.g. with vs without assemblage constraint, etc)  
Key parameters  
Run diversification models (RPANDA, BAMM, DDD, others)  
Evaluating success/performance of each method across sims in relation to “true” speciation and extinction rates observed over time  
