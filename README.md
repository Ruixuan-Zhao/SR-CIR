# Simulation and Data Analysis Code for “Towards the Revival of Survival Ratio”

#### Code Contributors: Xiaochuan Shi, Jingxin Yan, and Ruixuan Zhao

This supplementary file contains R code for the simulation studies and real data analyses for the manuscript.

## Software Requirements

The code uses R and the following R packages:

- `brm`
- `numDeriv`
- `extraDistr`
- `survival` 
- `survminer` 
- `ggplot2` 
- `ggpubr` 
- `parallel`

## R Files

### Method Files for Survival Ratio Model:

- `SR/getProb_sequent.R`: computes the sequential survival probabilities in forward fashion under the survival ratio and sequential odds-product models.
- `SR/MLE_point.R`: defines the likelihood and computes the maximum likelihood estimator.
- `SR/CallMLE.R`: wrapper for maximum likelihood estimation and inference.
- `SR/WrapResults.R`: organizes point estimates, standard errors, confidence intervals and p-values.
- `SR/DR_point_var.R`: defines the doubly robust estimating equation and computes the doubly robust estimator.
- `SR/CallDR.R`: wrapper for the doubly robust estimator and propensity-score model.
- `SR/DR_var_sandwich.R`: computes the sandwich variance for the doubly robust estimator.

### Method Files for Cumulative Incidence Ratio Model:

- `CIR/helpers_cir.R`: shared helper functions for matrix operations, time-grid handling, and command-line arguments.
- `CIR/getProb_sequent_cir.R`: computes sequential CDFs in backward fashion under the cumulative incidence ratio and sequential odds-product models.
- `CIR/MLE_point_cir.R`: defines the likelihood and computes the maximum likelihood estimator for cumulative incidence ratio model.
- `CIR/CallMLE_cir.R`: wrapper for maximum likelihood estimation and inference.
- `CIR/nuisance_cir_new.R`: fits the propensity-score and censoring nuisance models.
- `CIR/DR_point_cir.R`: defines the doubly robust estimating equation and point estimator.
- `CIR/DR_var_sandwich_cir.R`: computes the numerical Jacobian sandwich variance.
- `CIR/CIR_analytic_patch.R`: computes analytic nuisance derivatives and the full sandwich variance.

### Simulation Files

- `SR/simul_ex.R`: generates simulated data under the survival ratio model.
- `SR/run_simul_DR.R`: runs one simulation configuration and summarizes the results for survival ratio model.
- `SR/run_all_DR.R`: runs both simulation examples for sample sizes 500 and 1,000 for survival ratio model.
- `CIR/simul_ex_cir.R`: generates simulated data under the cumulative incidence ratio model.
- `CIR/run_simul_cir.R`: runs the simulations under the nuisance-model configurations for cumulative incidence ratio model.
- `CIR/cir_censoring_plot.R`: produces the boxplots from the raw simulation results.

### Real Data Files

The prostate cancer trial dataset (`prostate.rda`) is publicly available from the [Vanderbilt Biostatistics data repository](https://hbiostat.org/data/). Download `prostate.rda` and place it in the `SR/` directory before running the real data analysis.

- `SR/RealData_KM.R`: preprocesses the prostate cancer trial data and produces the Kaplan–Meier survival-curve figure.
- `SR/RealData_t_t2.R`: obtains the maximum likelihood and doubly robust estimators, and produces the tables and survival ratio figures.
- `SR/RealData_Addplots.R`: draws the figures in the Supplementary Materials.

### Run the Prostate Cancer Trial Analysis

Run the following script:

```r
setwd("SR")
dir.create("Figures", showWarnings = FALSE)
source("RealData_KM.R")
source("RealData_t_t2.R")
source("RealData_Addplots.R")
```

### Run Simulations for Survival Ratio Model

Run the following script from inside the `SR` directory to reproduce both simulation examples for sample sizes 500 and 1,000:

```bash
cd SR
Rscript run_all_DR.R
```

To run a single simulation configuration, use:

```bash
Rscript run_simul_DR.R --n=500 --example=1 --filename=sr_example1_n500
```

The current script uses 500 Monte Carlo replications. The simulation results are saved in `SR/Results/<filename>/` as:

- `Simu.RData`: replicate-level estimates and simulation objects.
- `simu.out`: simulation settings and summary results.

### Run Simulations for Cumulative Incidence Ratio Model

Run the following commands from inside the `CIR` directory to reproduce both simulation examples for sample sizes 500 and 1,000:

```bash
cd CIR

Rscript run_simul_cir.R --n=500  --nsim=500 --example=1 --seed=3047 --cores=4 --filename=cir_example1_n500
Rscript run_simul_cir.R --n=1000 --nsim=500 --example=1 --seed=3047 --cores=4 --filename=cir_example1_n1000
Rscript run_simul_cir.R --n=500  --nsim=500 --example=2 --seed=3047 --cores=4 --filename=cir_example2_n500
Rscript run_simul_cir.R --n=1000 --nsim=500 --example=2 --seed=3047 --cores=4 --filename=cir_example2_n1000
```

The `--cores` argument can be adjusted to match the number of available local CPU cores. The simulation results are saved in `CIR/Results/<filename>/` as:

- `raw_results.csv`: replicate-level estimates, standard errors, confidence intervals, and biases.
- `Simu.RData`: simulation objects and results.
- `simu.out`: simulation settings and summary results.



After completing all four simulation runs, generate the figure of boxplots using:

```bash
Rscript cir_censoring_plot.R
```



#### Use of AI Tools

The code was developed and reviewed by the authors. OpenAI ChatGPT was used to assist with code organization and documentation.
