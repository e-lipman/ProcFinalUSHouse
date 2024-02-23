# ProcFinalUSHouse

This repository produces the results and figures from [Explaining Differences in Voting Patterns Across Voting Domains Using Hierarchical Bayesian Models (Lipman, Moser, & Rodriguez, 2023)](<https://arxiv.org/abs/2312.15049>) 

## The `IdealPointsCompare` Submodule
The source code for the analysis is contained in a submodule called [IdealPointsCompare](<https://github.com/e-lipman/IdealPointsCompare/tree/11e6be83530b39dcfd18dc1bac1d410dc372154d>). 
To clone the reposity with the submodule, run the following:
```
git clone --recurse-submodules git@github.com:e-lipman/ProcFinalUSHouse.git
```

Documentation for the source code and data is contained in the README for the `IdealPointsCompare` repository

## Scripts
`run_all`

- Top level bash script to run all analyses and produce all figures. Global hyperparameters are set to small values to run very short test runs for a few sessions of the House. Parameters needed to replicate the results in the manuscript are indicated via comments in this script.

- Running the script for the full-length analyses requires parallelization at the level of congress and chain by modifying the for loops. Runtimes for the full model runs are on the order of days. 

`make_data_figures.R`: Makes Figures 1-3 (descriptive analysis)

`make_results_figures.R`: Makes Figures 4-8 (results)
