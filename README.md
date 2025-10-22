# Simulations/Empirical Scripts for Honours Research

This repository includes scripts for simulations and empirical analyses, including scripts running on local machines and high-performance computing clusters, scripts for generating results tables and figures.

This work is part of the Honours research project:

**Enhancing Forecasting Reconciliation: A Study of Alternative Covariance Estimators**

For more details, please refer to the accompanying thesis document [here](link_to_thesis_document).

## Quick Access to Key Scripts

Below are links to key scripts producing results for the thesis:

**Simulations**

-   [Simulation scripts (parallel locally)](sim/thesis_sim/main.R)
-   [Computing scores (probabilistic forecasts) for simulations](sim/thesis_sim/main_score.R)
-   [Wrangling simulation results and generating plots for thesis](sim/thesis_sim/read_results_thesis_sim.R)

**Australian Tourism Forecasting**

-   [Tourism forecasting scripts (on computing cluster)](tourism/job/tourism_run.R)
-   [Tourism forecasting scripts for multi-step covariance (on computing cluster)](tourism/job/tourism_run_newcov.R)
-   [Computing scores (probabilistic forecasts) for tourism](tourism/job/results/read_results_cluster.Rmd)
-   [Wrangling tourism results and generating plots for thesis (point reconciliation)](tourism/job/results_newcov/read_results_newcov.Rmd)
-   [Wrangling tourism results and generating plots for thesis (probabilistic reconciliation)](tourism/job/results/read_results_cluster.Rmd)

## Repository Structure

**Main directories:**

-   `sim/`: Contains scripts for running simulation locally and analysing results.
    -   Results files and scripts that generate plots for the thesis are in subfolder `thesis_sim/`.
    -   Others are old scripts
-   `tourism/`: Contains scripts for Australian tourism forecasting and results analysis.
    -   Scripts for analysing and plotting results for the thesis are in subfolders `job/results/`, `job/results_newcov/`.
    -   Scripts for computing clusters are `tourism_run_*.R` and `job.sh`.
    -   Results data files are stored separately due to large sizes.
    -   `forcetrend` files/folder are old scripts.

**Other directories (in archived):**

-   `archived/job/`: Contains initial scripts for running old simulations on computing clusters.
-   `archived/swiss_tourism/`: Contains scripts for Swiss tourism forecasting.