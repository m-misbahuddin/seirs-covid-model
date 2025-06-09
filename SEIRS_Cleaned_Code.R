# Check system environment
Sys.getenv("PATH")
system("gcc --version")
stopifnot(getRversion() >= "4.1")
stopifnot(packageVersion("pomp") >= "4.6")

# Load required libraries
library(pomp)        # For building and analyzing POMP models
library(lhs)         # For Latin Hypercube Sampling
library(dplyr)       # For data manipulation
library(ggplot2)     # For plotting
library(tidyverse)   # For data science tools
library(doParallel)  # For parallel processing
library(foreach)     # For parallel for-loop

# Set pomp compilation directory
options(pomp_cdir = "./tmp")

# Read data and preprocess
data <- read.csv("D:/BU/RA/IFH/IFH_pomp_data.csv")
data$Date <- as.Date(data$Diagnosis.Date, format = "%m/%d/%Y")
data$time <- as.numeric(data$Date - min(data$Date))
data$cases <- as.numeric(data$Count.of.Covid.Infected)
sum(is.na(data$cases))

# SEIR model step
seir_step <- Csnippet("
 double dN_SE = rbinom(S, 1 - exp(-Beta * I/N * dt));
 double dN_EI = rbinom(E, 1 - exp(-sigma * dt));
 double dN_IR = rbinom(I, 1 - exp(-gamma * dt));
 double dN_RS = rbinom(R, 1 - exp(-xi * dt));
 S = S + dN_RS - dN_SE;
 E = E + dN_SE - dN_EI;
 I = I + dN_EI - dN_IR;
 R = R + dN_IR - dN_RS;
 N = S + E + I + R;
")

# Initializer
rinit <- Csnippet("
 S = (int)round(S_init * N0);
 E = (int)round(E_init * (N0 - S));
 I = (int)round(I_init * (N0 - S - E));
 R = N0 - S - E - I;
 N = S + E + I + R;
")

# Parameter transformations
SEIRS_untrans <- Csnippet("
 T_Beta = log(Beta);
 T_sigma = log(sigma);
 T_sigmaSE = log(sigmaSE);
 T_gamma = log(gamma);
 T_xi = log(xi);
 T_rho = logit(rho);
 T_theta = log(theta);
 T_mu = log(mu);
 T_S_init = logit(S_init);
 T_E_init = logit(E_init);
 T_I_init = logit(I_init);
")

SEIRS_trans <- Csnippet("
 Beta = exp(T_Beta);
 sigma = exp(T_sigma);
 sigmaSE = exp(T_sigmaSE);
 gamma = exp(T_gamma);
 xi = exp(T_xi);
 rho = expit(T_rho);
 theta = exp(T_theta);
 mu = exp(T_mu);
 S_init = expit(T_S_init);
 E_init = expit(T_E_init);
 I_init = expit(T_I_init);
")

# Initial parameter guess
params <- c(
  Beta = 2.2/2.3, sigma = 1/5.2, gamma = 1/2.3, xi = 0.005,
  rho = 0.1, theta = 0.03, sigmaSE = 0.05, mu = 3.9e-05,
  N0 = 891705, S_init = 0.99, E_init = 0.009, I_init = 0.001
)

# Define parameter ranges for Latin Hypercube Sampling (LHS)
param_ranges <- matrix(c(
  0.1, 3, #1.0, 2.0,   # Beta range
  0.01, 2, #0.1, 0.9,   # sigma range
  0.01, 2, #0.1, 0.9,   # sigmaSE range
  0.05, 0.6,  # gamma range
  0.05, 0.1,  # xi range (waning immunity)
  0.7, 1.0,   # rho range
  5, 20,       # theta range
  1e-6, 1e-4,     # mu range (allowing some variation around the initial value)
  0.1, 0.999,       # S_init range
  0.0001, 0.1,       # E_init range
  0.0001, 0.1       # I_init range
), nrow = 11, byrow = TRUE)

n_samples <- 30
lhs_samples <- lhs::randomLHS(n_samples, 11)
lhs_params <- apply(lhs_samples, 2, function(x, range) {
  range[1] + x * (range[2] - range[1])
}, range = param_ranges)
colnames(lhs_params) <- c("Beta", "sigma", "sigmaSE", "gamma", "xi", "rho",
                          "theta", "mu", "S_init", "E_init", "I_init")

#----------------------------------------
# Parallelizing the 'mif2' Process Over Multiple LHS Samples
#----------------------------------------

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# 2) export all globals your workers need
clusterExport(cl, varlist = c(
  "data",
  "lhs_params",
  "seir_step",
  "rinit",
  "SEIRS_untrans",
  "SEIRS_trans",
  "n_samples"
))

# 3) benchmark start
start_time <- Sys.time()
cat("Start time:", start_time, "\n")

# 4) parallel MIF2 with tryCatch
mif_results <- foreach(
  i = 1:n_samples,
  .packages = "pomp"
) %dopar% {
  
  # build the full named vector of initial params
  initial_params <- c(
    lhs_params[i, "Beta"],
    lhs_params[i, "sigma"],
    lhs_params[i, "sigmaSE"],
    lhs_params[i, "gamma"],
    lhs_params[i, "xi"],
    lhs_params[i, "rho"],
    lhs_params[i, "theta"],
    lhs_params[i, "mu"],
    lhs_params[i, "S_init"],
    lhs_params[i, "E_init"],
    lhs_params[i, "I_init"]
  )
  print(names(initial_params))
  tryCatch({
    # rebuild the pomp model inside each worker
    pomp_model <- pomp(
      data      = data.frame(cases = data$Count.of.Covid.Infected, time = data$time),
      times     = "time",
      t0        = 0,
      rprocess  = euler(seir_step, delta.t = 1/10),
      rinit     = rinit,
      globals   = Csnippet("const double N0 = 891705;"),
      rmeasure  = Csnippet("cases = rnbinom_mu(theta, rho * I);"),
      dmeasure  = Csnippet("
        lik = dnbinom_mu(cases, theta, rho * I + 1e-6, give_log);
        if (!isfinite(lik)) lik = log(1e-100);
      "),
      statenames = c("S","E","I","R","N"),
      obsnames   = "cases",
      paramnames = c("Beta","sigma","sigmaSE","gamma","xi","rho","theta","mu",
                     "S_init","E_init","I_init"),
      params     = initial_params,
      partrans   = parameter_trans(toEst = SEIRS_untrans,
                                   fromEst = SEIRS_trans)
    )
    
    # run mif2
    set.seed(1234 + i)
    mif_fit <- mif2(
      pomp_model,
      Np    = 200,
      Nmif  = 500,
      params = initial_params,
      rw.sd = rw_sd(
        Beta    = 0.02, sigma   = 0.02, sigmaSE = 0.02, gamma = 0.02,
        xi      = 0.02, rho     = 0.02, theta   = 0.02, mu    = 0.02,
        S_init  = ivp(0.1), E_init = ivp(0.1), I_init = ivp(0.1)
      ),
      cooling.type       = "geometric",
      cooling.fraction.50 = 0.5
    )
    
    # return success
    list(
      iteration = i,
      logLik    = logLik(mif_fit),
      params    = coef(mif_fit)
    )
    
  }, error = function(e) {
    # return the error message
    list(
      iteration = i,
      error     = e$message
    )
  })
}

# 5) benchmark end
end_time <- Sys.time()
cat("End time:", end_time, "\nDuration:", end_time - start_time, "\n")

# 6) clean up
stopCluster(cl)

# End timing
end_time <- Sys.time()
total_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("End time:", format(end_time, "%H:%M:%S"), "\n")
cat("Duration:", floor(total_secs / 60), "minutes", round(total_secs %% 60, 1), "seconds\n")

# Extract successful runs
success <- vapply(mif_results, function(x) "logLik" %in% names(x), logical(1))

if (!any(success)) {
  errs <- vapply(mif_results, function(x) x$error %||% "<no error>", character(1))
  stop("All mif2 runs failed:\n", paste0("  iter ", seq_along(errs), ": ", errs, collapse = "\n"))
}

# Extract best result
best_idx <- which.max(vapply(mif_results[success], `[[`, numeric(1), "logLik"))
best_run <- mif_results[which(success)[best_idx]]
cat("\nBest MIF2 run (iteration", best_run$iteration, "):\n")
cat("Log-likelihood:", best_run$logLik, "\n")
print(best_run$params)

