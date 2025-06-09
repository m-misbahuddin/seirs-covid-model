# SEIRS MIF2 Parallel Fit

This repository contains R code for fitting a stochastic SEIRS (Susceptible-Exposed-Infectious-Recovered-Susceptible) model to time-series incidence data using the `pomp` package and iterated filtering (MIF2). The model is simulated and calibrated using parallelized Latin Hypercube Sampling (LHS) and likelihood maximization.

## Filename

**`SEIRS_MIF2_parallel_fit.R`**

## Purpose

* Simulate an SEIRS epidemiological model.
* Use Latin Hypercube Sampling to generate initial parameter sets.
* Perform parameter inference via MIF2 over multiple cores.
* Select the best-fit model based on log-likelihood.
* Visualize both pre-fit and post-fit model behavior.

## Dependencies

This script requires the following R packages:

* `pomp`
* `lhs`
* `dplyr`
* `ggplot2`
* `tidyverse`
* `doParallel`
* `foreach`

Ensure your system has:

* R >= 4.1
* `pomp` >= 4.6
* GCC compiler (required by `pomp` for compiling C snippets)

## Input Data

* A CSV file with COVID-19 incidence data.
* Required columns:

  * `Diagnosis.Date` (MM/DD/YYYY format)
  * `Count.of.Covid.Infected`

Update the file path in the script to point to your data:

```R
data <- read.csv("D:/BU/RA/IFH/IFH_pomp_data.csv")
```

## Model Description

* SEIRS structure with the following transitions:

  * S -> E (infection)
  * E -> I (latency)
  * I -> R (recovery)
  * R -> S (waning immunity)
* Stochastic transitions using binomial approximations.
* Negative Binomial measurement model.

## Procedure

1. Load and preprocess data.
2. Define SEIRS model using C snippets.
3. Generate parameter sets with Latin Hypercube Sampling.
4. Run parallel `mif2` inference over multiple parameter sets.
5. Select best-fit model based on log-likelihood.
6. Simulate model with best-fit parameters.
7. Plot observed vs simulated cases and compartment dynamics.

## Outputs

* Pre- and post-MIF2 simulation plots.
* Trajectories for S, E, I, R compartments.
* Simulated vs. observed case comparisons.

## Notes

* If all MIF2 runs fail, detailed error messages are printed.
* You can adjust the number of particles (`Np`), MIF iterations (`Nmif`), and random walk standard deviations (`rw.sd`) as needed.

## License

Licensed under the MIT License. See the LICENSE file for details.

For academic or research use, please cite:

> King AA, Nguyen D, Ionides EL. Statistical inference for partially observed Markov processes via the R package pomp. J Stat Softw. 2016;69(12):1-43.
