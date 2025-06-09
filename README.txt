# SEIRS COVID-19 Modeling using `pomp` in R

This project implements a **Stochastic SEIRS (Susceptible-Exposed-Infectious-Recovered-Susceptible)** model for simulating and fitting COVID-19 transmission dynamics using the **`pomp`** package in R. The workflow includes parameter sampling, model simulation, parameter estimation via iterated filtering (`mif2`), and visualization.

---

## 📦 Dependencies

Ensure the following R packages are installed:

```r
install.packages(c("pomp", "lhs", "doParallel", "doFuture", "dplyr", "ggplot2", "tidyverse"))
```

Additional tools:
- A working C compiler (e.g., GCC) is required for compiling `Csnippet` code.
- R version >= 4.1
- `pomp` version >= 4.6

---

## 📁 Data

Your dataset should be a CSV file with at least:
- `Diagnosis.Date`: the date of infection diagnosis (MM/DD/YYYY format)
- `Count.of.Covid.Infected`: daily case counts

Update the file path in this line:

```r
data <- read.csv("D:/BU/RA/IFH/IFH_pomp_data.csv")
```

---

## 🧠 Model Summary

- **Compartments**: S (Susceptible), E (Exposed), I (Infectious), R (Recovered)
- **Transitions**: Implemented using stochastic Euler approximations with `Csnippet`.
- **Birth/death processes**: Optional and currently commented.
- **Noise**: Includes white noise via Gamma-distributed random variables.

---

## ⚙️ Workflow

### 1. **Preprocessing**
- Converts date to numeric time.
- Initializes `cases` as numeric.

### 2. **Model Definition**
- Uses `pomp()` to define:
  - `rprocess`: SEIRS process model
  - `rinit`: Initial condition generator
  - `rmeasure/dmeasure`: Negative binomial observation model
  - `parameter_trans`: Log and logit transformations

### 3. **Sampling Parameters**
- Uses **Latin Hypercube Sampling (LHS)** to explore parameter space.
- Adjust `n_samples` for how many parameter sets you want.

### 4. **Simulation**
- Simulates one or more realizations of the model.
- Plots SEIRS compartment trajectories.

### 5. **Iterated Filtering (mif2)**
- **Initial testing**: Ensures model works with a single iteration.
- **Parallel execution**: Uses `foreach` and `doParallel` to fit multiple parameter sets in parallel.
- **Sequential execution**: Alternative fallback for environments where parallelization fails.

### 6. **Model Evaluation**
- Selects the best parameter set based on log-likelihood.
- Simulates model using best parameters and compares to observed cases.
- Uses `pfilter()` to estimate final log-likelihood.

---

## 📈 Outputs

- Plots of compartmental dynamics
- Comparison between observed and simulated case trajectories
- Log-likelihoods from `mif2` and `pfilter`
- Best-fit parameter estimates

---

## 🛠 Notes

- You can customize the number of particles (`Np`), iterations (`Nmif`), and cooling settings in the `mif2` routine.
- Consider increasing `Np` and `Nmif` for better estimation accuracy at the cost of computation time.
- Check system environment for C compiler setup using:

```r
Sys.getenv("PATH")
system("gcc --version")
```

---

## ✅ Example: Run Basic Model Fit

```r
# Run single mif2 iteration
mif_fit <- mif2(pomp_model, Np = 200, Nmif = 10, params = initial_params, ...)
logLik(mif_fit)
coef(mif_fit)
```

---

## 📚 References

- King, A. A., Nguyen, D., & Ionides, E. L. (2016). A *pomp* package for R.
- [pomp documentation](https://kingaa.github.io/pomp/)
