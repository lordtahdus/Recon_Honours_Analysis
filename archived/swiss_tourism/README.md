# Modeling the uncertainty on the covariance matrix for probabilistic forecast reconciliation

This repository accompanies the paper *Modeling the uncertainty on the covariance matrix for probabilistic forecast reconciliation*. The paper introduces **t-Reconciliation**, a Bayesian approach to hierarchical time series forecasting that explicitly accounts for uncertainty in the residuals covariance matrix.

By placing an Inverse-Wishart prior on the covariance matrix, the method incorporates parameter uncertainty and leads to a closed-form solution for the reconciled forecasts. The resulting predictive distribution follows a multivariate t-distribution, ensuring probabilistic coherence across the hierarchy.

The repository includes:

- The **data** used in the paper.
- Code to reproduce all **plots and results** from the paper.
- A **simple illustrative example** of t-Reconciliation in action.

## Main repository elements

- `main.R`: Reproduces empirical results on real-world datasets.
- `main_sim.R`: Reproduces simulation study results.
- `t-Reconciliation.Rmd`: Interactive simple example for t-Reconciliation using the Swiss tourism dataset.


## Datasets used 

The paper utilizes three datasets: *Swiss tourism*, *Australian tourism – M*, and *Australian tourism – Q*. 
Below is a brief description of each dataset along with the associated data files.

*Swiss tourism*: Swiss monthly overnight stays (2005–2025), with a one-level hierarchy:
- National total
- 26 Cantons (cross-sectional)

Data files:
- `data_Swiss_tourism.RDS`
- `A_Swiss_tourism.RDS` (aggregation matrix)

Source:
- derived from `SwissTourism.csv`

*Australian tourism - M*: Australian monthly overnight trips (1998–2016), with a three-levels hierarchy:
- National total
- 7 States
- 27 Zones
- 76 Regions

Data files:
- `data_Australian_tourism_zone.RDS`
- `A_Australian_tourism_zone.RDS` (aggregation matrix)

Source:
- derived from `Regions.csv`

*Australian tourism - Q*: Australian quarterly overnight trips (1998–2017), with a two-levels hierarchy:
- National total
- 8 States 
- 76 Regions

> *Note: In this dataset, ACT is included as a state. In contrast, for Australian tourism – M, ACT is part of the zones.*

Data files:
- `data_Australian_tourism_no_zone.RDS`
- `A_Australian_tourism_no_zone.RDS` (aggregation matrix)

Source:
- obtained from the `tsibble` R package.

> *Note: for reproducing all the results the first step is to clone the repository.*

## How to Run the simple example for the Swiss tourism dataset

1. Open `t-Reconciliation.Rmd` in RStudio.

2. Run the code chunks interactively or click the Knit button to render the full analysis.

3. The following R packages will be installed by running the vignette:

```r
  install.packages(c(
  "forecast",
  "bayesRecon",
  "nloptr",
  "scoringRules",
  "mvtnorm",
  "ggplot2",
  "dplyr",
  "tidyr"
  ))
```


## How to reproduce the paper plots for the real datasets

1. Open `main.R` in RStudio.

2. Install the required R packages (if not already installed):

```r
   install.packages(c(
  "forecast",
  "bayesRecon",
  "nloptr",
  "scoringRules",
  "mvtnorm",
  "ggplot2",
  "dplyr",
  "tidyr",
  "tsibble",
  "tseries",
  "fpp3",
  "future.apply",
  "progressr",
  "parallel",
  "patchwork",
  "tsutils"
))
```

3. Run the script and select which dataset(s) you want to process. To reproduce **all** plots from the paper, make sure to select **all** datasets when prompted.

4. Read the table(s) results in the Console and check the `Paper_plots/` folder for generated figures.


## How to reproduce the paper plots for the simulations

1. Open `main_sim.R` in RStudio.

2. Install the required R packages (if not already installed):

```r
   install.packages(c(
  "forecast",
  "bayesRecon",
  "nloptr",
  "scoringRules",
  "mvtnorm",
  "ggplot2",
  "dplyr",
  "tidyr",
  "tsibble",
  "tseries",
  "fpp3",
  "future.apply",
  "progressr",
  "parallel",
  "patchwork",
  "tsutils"
))
```

3. Run the script.

4. Read the table(s) results in the Console and check the `Paper_plots/` folder for generated figures.


