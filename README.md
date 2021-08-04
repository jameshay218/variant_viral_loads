# Viral loads observed under competing strain dynamics
This Git repository accompanies the preprint available [here](https://www.medrxiv.org/content/10.1101/2021.07.27.21261224v1). In summary, these analyses investigate how viral loads in an infected population change over the course of a two-strain epidemic. 

**NOTE TO REVIEWERS:** please do report an issue or send me an email if you run into issues (assuming this acceptable in with blinded peer review).

## Test status
This code has been tested on two separate machines, both running macOS Big Sur version 11.1. The code was developed using R version 4.0.2, but should run on early versions. No non-standard hardware is required. Install time should take around 10-15 minutes, accounting for reading instructions and download times. Installation may take a bit longer if you need to install a compiler and configure your Rtools.

The `main.R` script will take ~30 minutes to run, but the run time can be shortened by reducing the number of repeats/trials being simulated. The `main_infer_dynamics.R` script will take a couple of hours, but can be shortened by reducing the number of MCMC iterations and adaptive period in the `mcmc_pars_pt` variable.

## Setup
The two main analyses scripts, `main.R` and `main_infer_dynamics.R` can be run simply by sourcing the scripts. However, there are two R packages on Github that must first be installed:

1. `lazymcmc` is used for the MCMC procedure. This is can be installed using `devtools::install_github("jameshay218/lazymcmc@parallel_tempering")`, and is available [here](https://github.com/jameshay218/lazymcmc). NOTE it is important to use the `parallel_tempering` branch here.
2. `virosolver` is used for the viral kinetics model and Ct-based inference. This can be installed using `devtools::install_github("jameshay218/virosolver")`, with further installation and usage instructions available [here](https://github.com/jameshay218/virosolver/).

There are a number of package from CRAN also required:

``` r
c("tidyverse","ggthemes","ggpubr","patchwork","moments",
"fitdistrplus","deSolve","doParallel","coda")
```

Note that this code uses compiled C++ code (using Rcpp), so you will need to have a C++ compiler (Rtools on Windows, Xcode on Mac) installed.

## Repository structure

- All R functions for plotting, simulation and analyses are in the `code` directory. 
- Functions not found in there will be part of `lazymcmc` or `virosolver`. 
- Figures will be saved to the `figures` folder. 
- Tables of parameters used for simulation and `virosolver` fitting are in the `pars` folder.
- Scripts needed to run the full analyses in Figure 4 and Figures S8-11 are in the `scripts` folder.
- In general, MCMC chains will be saved to the `chains` folder when the script is run. But you can change this to any file path you like.

## Part 1: simulating and comparing Ct values
To regenerate the analyses in Figures 1-3 and Figures S1-S7, just source the script `main.R`. You will need to change the `HOME_WD` variable at the top of the script, but after than this should run straight out of the box. Please do log an issue if you encounter any problems. Note that there are many places in the script to change parameter assumptions: for example, you can change all of the viral kinetics parameters at the top of the script, or the assumed sample times on eg. line 116 and line 164.

## Part 2: inferring differences in epidemic dynamics and viral kinetics
An example run of simulating a dataset and using `virosolver` to re-estimate variant-specific viral kinetics and epidemic dynamics is given in `main_infer_dynamics.R`. As above, just change the `HOME_WD` variable at the top of the script and step through the code. This script will simulate a dataset with cross-sectional samples for two variants with different dynamics, and then use `virosolver` to estimate the posterior distribution of the viral kinetics parameters and epidemic growth rate conditional on the observed data. The user can adjust the priors assumed for the model. Finally, the script will read in the MCMC chains and plot some posterior diagnostics/results. Note that with these chain lengths and the use of parallel tempering, convergence has been consistently reliable. However, users are encouraged to check MCMC trace plots, Gelman-Rubin diagnostics and effective sample sizes.

Rerunning all of the analyses for this section requires the use of a computing cluster due to run time of re-fitting the model for 100 simulated
datasets. Although tailored to the Harvard School of Public Health computing cluster, these scripts can be found in the `scripts` folder:

- `scripts/virosolver_different_grs.sh`, which calls `code/virosolver_different_grs_cluster.R`. Note that this R script is almost identical to `main_infer_dynamics.R`. This underpins the analyses where simultaneously seedded variants with different growth rates are compared early int he epidemic.
- `scripts/virosolver_late_grs.sh`, which calls `code/virosolver_late_grs_cluster.R`. As above, but assuming that samples are obtained late in the epidemic (day 270) and the new variant is seeded late.

Finally, the script `virosolver_gr_comparison_figures.R` combines the outputs of the cluster runs to produce the figures in the paper.