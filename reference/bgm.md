# Bayesian Estimation or Edge Selection for Markov Random Fields

The `bgm` function estimates the pseudoposterior distribution of
category thresholds (main effects) and pairwise interaction parameters
of a Markov Random Field (MRF) model for binary and/or ordinal
variables. Optionally, it performs Bayesian edge selection using
spike-and-slab priors to infer the network structure.

## Usage

``` r
bgm(
  x,
  variable_type = "ordinal",
  baseline_category,
  iter = 1000,
  warmup = 1000,
  pairwise_scale = 2.5,
  main_alpha = 0.5,
  main_beta = 0.5,
  edge_selection = TRUE,
  edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
  inclusion_probability = 0.5,
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  beta_bernoulli_alpha_between = 1,
  beta_bernoulli_beta_between = 1,
  dirichlet_alpha = 1,
  lambda = 1,
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis", "hamiltonian-mc"),
  target_accept,
  hmc_num_leapfrogs = 100,
  nuts_max_depth = 10,
  learn_mass_matrix = FALSE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  interaction_scale,
  burnin,
  save,
  threshold_alpha,
  threshold_beta
)
```

## Arguments

- x:

  A data frame or matrix with `n` rows and `p` columns containing binary
  and ordinal responses. Variables are automatically recoded to
  non-negative integers (`0, 1, ..., m`). For regular ordinal variables,
  unobserved categories are collapsed; for Blume–Capel variables, all
  categories are retained.

- variable_type:

  Character or character vector. Specifies the type of each variable in
  `x`. Allowed values: `"ordinal"` or `"blume-capel"`. Binary variables
  are automatically treated as `"ordinal"`. Default: `"ordinal"`.

- baseline_category:

  Integer or vector. Baseline category used in Blume–Capel variables.
  Can be a single integer (applied to all) or a vector of length `p`.
  Required if at least one variable is of type `"blume-capel"`.

- iter:

  Integer. Number of post–burn-in iterations (per chain). Default:
  `1e3`.

- warmup:

  Integer. Number of warmup iterations before collecting samples. A
  minimum of 1000 iterations is enforced, with a warning if a smaller
  value is requested. Default: `1e3`.

- pairwise_scale:

  Double. Scale of the Cauchy prior for pairwise interaction parameters.
  Default: `2.5`.

- main_alpha, main_beta:

  Double. Shape parameters of the beta-prime prior for threshold
  parameters. Must be positive. If equal, the prior is symmetric.
  Defaults: `main_alpha = 0.5` and `main_beta = 0.5`.

- edge_selection:

  Logical. Whether to perform Bayesian edge selection. If `FALSE`, the
  model estimates all edges. Default: `TRUE`.

- edge_prior:

  Character. Specifies the prior for edge inclusion. Options:
  `"Bernoulli"`, `"Beta-Bernoulli"`, or `"Stochastic-Block"`. Default:
  `"Bernoulli"`.

- inclusion_probability:

  Numeric scalar. Prior inclusion probability of each edge (used with
  the Bernoulli prior). Default: `0.5`.

- beta_bernoulli_alpha, beta_bernoulli_beta:

  Double. Shape parameters for the beta distribution in the
  Beta–Bernoulli and the Stochastic-Block priors. Must be positive. For
  the Stochastic-Block prior these are the shape parameters for the
  within-cluster edge inclusion probabilities. Defaults:
  `beta_bernoulli_alpha = 1` and `beta_bernoulli_beta = 1`.

- beta_bernoulli_alpha_between, beta_bernoulli_beta_between:

  Double. Shape parameters for the between-cluster edge inclusion
  probabilities in the Stochastic-Block prior. Must be positive.
  Default: `beta_bernoulli_alpha_between = 1` and
  `beta_bernoulli_beta_between = 1`

- dirichlet_alpha:

  Double. Concentration parameter of the Dirichlet prior on block
  assignments (used with the Stochastic Block model). Default: `1`.

- lambda:

  Double. Rate of the zero-truncated Poisson prior on the number of
  clusters in the Stochastic Block Model. Default: `1`.

- na_action:

  Character. Specifies missing data handling. Either `"listwise"` (drop
  rows with missing values) or `"impute"` (perform single imputation
  during sampling). Default: `"listwise"`.

- update_method:

  Character. Specifies how the MCMC sampler updates the model
  parameters:

  "adaptive-metropolis"

  :   Componentwise adaptive Metropolis–Hastings with Robbins–Monro
      proposal adaptation.

  "hamiltonian-mc"

  :   Hamiltonian Monte Carlo with fixed path length (number of leapfrog
      steps set by `hmc_num_leapfrogs`).

  "nuts"

  :   The No-U-Turn Sampler, an adaptive form of HMC with dynamically
      chosen trajectory lengths.

  Default: `"nuts"`.

- target_accept:

  Numeric between 0 and 1. Target acceptance rate for the sampler.
  Defaults are set automatically if not supplied: `0.44` for adaptive
  Metropolis, `0.65` for HMC, and `0.80` for NUTS.

- hmc_num_leapfrogs:

  Integer. Number of leapfrog steps for Hamiltonian Monte Carlo. Must be
  positive. Default: `100`.

- nuts_max_depth:

  Integer. Maximum tree depth in NUTS. Must be positive. Default: `10`.

- learn_mass_matrix:

  Logical. If `TRUE`, adapt a diagonal mass matrix during warmup
  (HMC/NUTS only). If `FALSE`, use the identity matrix. Default:
  `FALSE`.

- chains:

  Integer. Number of parallel chains to run. Default: `4`.

- cores:

  Integer. Number of CPU cores for parallel execution. Default:
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).

- display_progress:

  Logical. Whether to show a progress bar during sampling. Default:
  `TRUE`.

- seed:

  Optional integer. Random seed for reproducibility. Must be a single
  non-negative integer.

- interaction_scale, burnin, save, threshold_alpha, threshold_beta:

  \`r lifecycle::badge("deprecated")\` Deprecated arguments as of
  \*\*bgms 0.1.6.0\*\*. Use \`pairwise_scale\`, \`warmup\`,
  \`main_alpha\`, and \`main_beta\` instead.

## Value

A list of class `"bgms"` with posterior summaries, posterior mean
matrices, and access to raw MCMC draws. The object can be passed to
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`coef()`](https://rdrr.io/r/stats/coef.html).

Main components include:

- `posterior_summary_main`: Data frame with posterior summaries (mean,
  sd, MCSE, ESS, Rhat) for category threshold parameters.

- `posterior_summary_pairwise`: Data frame with posterior summaries for
  pairwise interaction parameters.

- `posterior_summary_indicator`: Data frame with posterior summaries for
  edge inclusion indicators (if `edge_selection = TRUE`).

- `posterior_mean_main`: Matrix of posterior mean thresholds (rows =
  variables, cols = categories or parameters).

- `posterior_mean_pairwise`: Symmetric matrix of posterior mean pairwise
  interaction strengths.

- `posterior_mean_indicator`: Symmetric matrix of posterior mean
  inclusion probabilities (if edge selection was enabled).

- Additional summaries returned when `edge_prior = "Stochastic-Block"`.
  For more details about this prior see Sekulovski et al. (2025) .

  - `posterior_summary_pairwise_allocations`: Data frame with posterior
    summaries (mean, sd, MCSE, ESS, Rhat) for the pairwise cluster
    co-occurrence of the nodes. This serves to indicate whether the
    estimated posterior allocations,co-clustering matrix and posterior
    cluster probabilities (see blow) have converged.

  - `posterior_coclustering_matrix`: a symmetric matrix of pairwise
    proportions of occurrence of every variable. This matrix can be
    plotted to visually inspect the estimated number of clusters and
    visually inspect nodes that tend to switch clusters.

  - `posterior_mean_allocations`: A vector with the posterior mean of
    the cluster allocations of the nodes. This is calculated using the
    method proposed in Dahl (2009) .

  - `posterior_mode_allocations`: A vector with the posterior mode of
    the cluster allocations of the nodes.

  - `posterior_num_blocks`: A data frame with the estimated posterior
    inclusion probabilities for all the possible number of clusters.

- `raw_samples`: A list of raw MCMC draws per chain:

  - `main`:

    List of main effect samples.

  - `pairwise`:

    List of pairwise effect samples.

  - `indicator`:

    List of indicator samples (if edge selection enabled).

  - `allocations`:

    List of cluster allocations (if SBM prior used).

  - `nchains`:

    Number of chains.

  - `niter`:

    Number of post–warmup iterations per chain.

  - `parameter_names`:

    Named lists of parameter labels.

- `arguments`: A list of function call arguments and metadata (e.g.,
  number of variables, warmup, sampler settings, package version).

The [`summary()`](https://rdrr.io/r/base/summary.html) method prints
formatted posterior summaries, and
[`coef()`](https://rdrr.io/r/stats/coef.html) extracts posterior mean
matrices.

NUTS diagnostics (tree depth, divergences, energy, E-BFMI) are included
in `fit$nuts_diag` if `update_method = "nuts"`.

## Details

This function models the joint distribution of binary and ordinal
variables using a Markov Random Field, with support for edge selection
through Bayesian variable selection. The statistical foundation of the
model is described in Marsman et al. (2025) , where the ordinal MRF
model and its Bayesian estimation procedure were first introduced. While
the implementation in bgms has since been extended and updated (e.g.,
alternative priors, parallel chains, HMC/NUTS warmup), it builds on that
original framework.

Key components of the model are described in the sections below.

## Ordinal Variables

The function supports two types of ordinal variables:

**Regular ordinal variables**: Assigns a category threshold parameter to
each response category except the lowest. The model imposes no
additional constraints on the distribution of category responses.

**Blume-Capel ordinal variables**: Assume a baseline category (e.g., a
“neutral” response) and score responses by distance from this baseline.
Category thresholds are modeled as:

\$\$\mu\_{c} = \alpha \cdot c + \beta \cdot (c - b)^2\$\$

where:

- \\\mu\_{c}\\: category threshold for category \\c\\

- \\\alpha\\: linear trend across categories

- \\\beta\\: preference toward or away from the baseline

  - If \\\beta \< 0\\, the model favors responses near the baseline
    category;

  - if \\\beta \> 0\\, it favors responses farther away (i.e.,
    extremes).

- \\b\\: baseline category

## Edge Selection

When `edge_selection = TRUE`, the function performs Bayesian variable
selection on the pairwise interactions (edges) in the MRF using
spike-and-slab priors.

Supported priors for edge inclusion:

- **Bernoulli**: Fixed inclusion probability across edges.

- **Beta-Bernoulli**: Inclusion probability is assigned a Beta prior
  distribution.

- **Stochastic-Block**: Cluster-based edge priors with Beta, Dirichlet,
  and Poisson hyperpriors.

All priors operate via binary indicator variables controlling the
inclusion or exclusion of each edge in the MRF.

## Prior Distributions

- **Pairwise effects**: Modeled with a Cauchy (slab) prior.

- **Main effects**: Modeled using a beta-prime distribution.

- **Edge indicators**: Use either a Bernoulli, Beta-Bernoulli, or
  Stochastic-Block prior (as above).

## Sampling Algorithms and Warmup

Parameters are updated within a Gibbs framework, but the conditional
updates can be carried out using different algorithms:

- **Adaptive Metropolis–Hastings**: Componentwise random–walk updates
  for main effects and pairwise effects. Proposal standard deviations
  are adapted during burn–in via Robbins–Monro updates toward a target
  acceptance rate.

- **Hamiltonian Monte Carlo (HMC)**: Joint updates of all parameters
  using fixed–length leapfrog trajectories. Step size is tuned during
  warmup via dual–averaging; the diagonal mass matrix can also be
  adapted if `learn_mass_matrix = TRUE`.

- **No–U–Turn Sampler (NUTS)**: An adaptive extension of HMC that
  dynamically chooses trajectory lengths. Warmup uses a staged
  adaptation schedule (fast–slow–fast) to stabilize step size and, if
  enabled, the mass matrix.

When `edge_selection = TRUE`, updates of edge–inclusion indicators are
carried out with Metropolis–Hastings steps. These are switched on after
the core warmup phase, ensuring that graph updates occur only once the
samplers’ tuning parameters (step size, mass matrix, proposal SDs) have
stabilized.

After warmup, adaptation is disabled. Step size and mass matrix are
fixed at their learned values, and proposal SDs remain constant.

## Warmup and Adaptation

The warmup procedure in `bgm` is based on the multi–stage adaptation
schedule used in Stan (Stan Development Team 2023) . Warmup iterations
are split into several phases:

- **Stage 1 (fast adaptation)**: A short initial interval where only
  step size (for HMC/NUTS) is adapted, allowing the chain to move
  quickly toward the typical set.

- **Stage 2 (slow windows)**: A sequence of expanding, memoryless
  windows where both step size and, if `learn_mass_matrix = TRUE`, the
  diagonal mass matrix are adapted. Each window ends with a reset of the
  dual–averaging scheme for improved stability.

- **Stage 3a (final fast interval)**: A short interval at the end of the
  core warmup where the step size is adapted one final time.

- **Stage 3b (proposal–SD tuning)**: Only active when
  `edge_selection = TRUE` under HMC/NUTS. In this phase, Robbins–Monro
  adaptation of proposal standard deviations is performed for the
  Metropolis steps used in edge–selection moves.

- **Stage 3c (graph selection warmup)**: Also only relevant when
  `edge_selection = TRUE`. At the start of this phase, a random graph
  structure is initialized, and Metropolis–Hastings updates for edge
  inclusion indicators are switched on.

When `edge_selection = FALSE`, the total number of warmup iterations
equals the user–specified `burnin`. When `edge_selection = TRUE` and
`update_method` is `"nuts"` or `"hamiltonian-mc"`, the schedule
automatically appends additional Stage-3b and Stage-3c intervals, so the
total warmup is strictly greater than the requested `burnin`.

After all warmup phases, the sampler transitions to the sampling phase
with adaptation disabled. Step size and mass matrix (for HMC/NUTS) are
fixed at their learned values, and proposal SDs remain constant.

This staged design improves stability of proposals and ensures that both
local parameters (step size) and global parameters (mass matrix,
proposal SDs) are tuned before collecting posterior samples.

For adaptive Metropolis–Hastings runs, step size and mass matrix
adaptation are not relevant. Proposal SDs are tuned continuously during
burn–in using Robbins–Monro updates, without staged fast/slow intervals.

## Missing Data

If `na_action = "listwise"`, observations with missing values are
removed. If `na_action = "impute"`, missing values are imputed during
Gibbs sampling.

## References

Dahl DB (2009). “Modal clustering in a class of product partition
models.” *Bayesian Analysis*, **4**(2), 243–264.
[doi:10.1214/09-BA409](https://doi.org/10.1214/09-BA409) .  
  
Marsman M, van den Bergh D, Haslbeck JMB (2025). “Bayesian analysis of
the ordinal Markov random field.” *Psychometrika*, **90**, 146–-182.  
  
Sekulovski N, Arena G, Haslbeck JMB, Huth KBS, Friel N, Marsman M
(2025). “A Stochastic Block Prior for Clustering in Graphical Models.”
*Retrieved from <https://osf.io/preprints/psyarxiv/29p3m_v1>*. OSF
preprint.  
  
Stan Development Team (2023). *Stan Modeling Language Users Guide and
Reference Manual*. Version 2.33, <https://mc-stan.org/docs/>.

## See also

[`vignette("intro", package = "bgms")`](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/intro.md)
for a worked example.

## Examples

``` r
# \donttest{
# Run bgm on subset of the Wenchuan dataset
fit = bgm(x = Wenchuan[, 1:5])
#> Warning: There were 7 rows with missing observations in the input matrix x.
#> Since na_action = listwise these rows were excluded from the analysis.
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2200 (2.3%)
#> Chain 2 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 56/2200 (2.5%)
#> Chain 3 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 60/2200 (2.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 43/2200 (2.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 209/8800 (2.4%)
#> Elapsed: 2s | ETA: 1m 22s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2200 (4.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 104/2200 (4.7%)
#> Chain 3 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 114/2200 (5.2%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 85/2200 (3.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 403/8800 (4.6%)
#> Elapsed: 4s | ETA: 1m 23s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 150/2200 (6.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 154/2200 (7.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 164/2200 (7.5%)
#> Chain 4 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 127/2200 (5.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 595/8800 (6.8%)
#> Elapsed: 6s | ETA: 1m 22s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 200/2200 (9.1%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 196/2200 (8.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 215/2200 (9.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 161/2200 (7.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 772/8800 (8.8%)
#> Elapsed: 8s | ETA: 1m 23s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 250/2200 (11.4%)
#> Chain 2 (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 245/2200 (11.1%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 262/2200 (11.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 208/2200 (9.5%)
#> Total   (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 965/8800 (11.0%)
#> Elapsed: 10s | ETA: 1m 21s
#> Chain 1 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2200 (13.6%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 282/2200 (12.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 306/2200 (13.9%)
#> Chain 4 (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 247/2200 (11.2%)
#> Total   (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1135/8800 (12.9%)
#> Elapsed: 12s | ETA: 1m 21s
#> Chain 1 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 350/2200 (15.9%)
#> Chain 2 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 336/2200 (15.3%)
#> Chain 3 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 352/2200 (16.0%)
#> Chain 4 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 293/2200 (13.3%)
#> Total   (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1331/8800 (15.1%)
#> Elapsed: 14s | ETA: 1m 18s
#> Chain 1 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/2200 (18.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 389/2200 (17.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 394/2200 (17.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 328/2200 (14.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1511/8800 (17.2%)
#> Elapsed: 16s | ETA: 1m 17s
#> Chain 1 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 450/2200 (20.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 446/2200 (20.3%)
#> Chain 3 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 451/2200 (20.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 384/2200 (17.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1731/8800 (19.7%)
#> Elapsed: 18s | ETA: 1m 13s
#> Chain 1 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/2200 (22.7%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 505/2200 (23.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 503/2200 (22.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 434/2200 (19.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1942/8800 (22.1%)
#> Elapsed: 20s | ETA: 1m 10s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 550/2200 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 554/2200 (25.2%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 553/2200 (25.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 488/2200 (22.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2145/8800 (24.4%)
#> Elapsed: 22s | ETA: 1m 8s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 600/2200 (27.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 603/2200 (27.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 617/2200 (28.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 545/2200 (24.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2365/8800 (26.9%)
#> Elapsed: 24s | ETA: 1m 5s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 650/2200 (29.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 664/2200 (30.2%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 667/2200 (30.3%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 599/2200 (27.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2580/8800 (29.3%)
#> Elapsed: 26s | ETA: 1m 2s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/2200 (31.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 716/2200 (32.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 721/2200 (32.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 662/2200 (30.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2799/8800 (31.8%)
#> Elapsed: 28s | ETA: 1m
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 750/2200 (34.1%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 785/2200 (35.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 770/2200 (35.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 719/2200 (32.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3024/8800 (34.4%)
#> Elapsed: 30s | ETA: 57s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2200 (36.4%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 826/2200 (37.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 811/2200 (36.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 769/2200 (35.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3206/8800 (36.4%)
#> Elapsed: 31s | ETA: 54s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 850/2200 (38.6%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 865/2200 (39.3%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 851/2200 (38.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 810/2200 (36.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3376/8800 (38.4%)
#> Elapsed: 33s | ETA: 53s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/2200 (40.9%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 918/2200 (41.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 894/2200 (40.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 850/2200 (38.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 3562/8800 (40.5%)
#> Elapsed: 35s | ETA: 51s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 950/2200 (43.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 968/2200 (44.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 947/2200 (43.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 898/2200 (40.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 3763/8800 (42.8%)
#> Elapsed: 37s | ETA: 50s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1000/2200 (45.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1024/2200 (46.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1004/2200 (45.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 955/2200 (43.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 3983/8800 (45.3%)
#> Elapsed: 39s | ETA: 47s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1050/2200 (47.7%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1072/2200 (48.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1047/2200 (47.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1005/2200 (45.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4174/8800 (47.4%)
#> Elapsed: 41s | ETA: 45s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/2200 (50.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1139/2200 (51.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 1104/2200 (50.2%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1064/2200 (48.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 4407/8800 (50.1%)
#> Elapsed: 43s | ETA: 43s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2200 (52.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 1175/2200 (53.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1137/2200 (51.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1093/2200 (49.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4555/8800 (51.8%)
#> Elapsed: 44s | ETA: 41s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1200/2200 (54.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1211/2200 (55.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 1175/2200 (53.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1147/2200 (52.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4733/8800 (53.8%)
#> Elapsed: 45s | ETA: 39s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1250/2200 (56.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1254/2200 (57.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1222/2200 (55.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1199/2200 (54.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 4925/8800 (56.0%)
#> Elapsed: 47s | ETA: 37s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1300/2200 (59.1%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1309/2200 (59.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 1267/2200 (57.6%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1245/2200 (56.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 5121/8800 (58.2%)
#> Elapsed: 49s | ETA: 35s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1350/2200 (61.4%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1356/2200 (61.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1322/2200 (60.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1298/2200 (59.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 5326/8800 (60.5%)
#> Elapsed: 50s | ETA: 33s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1400/2200 (63.6%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1412/2200 (64.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1360/2200 (61.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1345/2200 (61.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 5517/8800 (62.7%)
#> Elapsed: 52s | ETA: 31s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1450/2200 (65.9%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1460/2200 (66.4%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1413/2200 (64.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1386/2200 (63.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5709/8800 (64.9%)
#> Elapsed: 53s | ETA: 29s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1500/2200 (68.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1510/2200 (68.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1456/2200 (66.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1432/2200 (65.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5898/8800 (67.0%)
#> Elapsed: 55s | ETA: 27s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1550/2200 (70.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1555/2200 (70.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1497/2200 (68.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1478/2200 (67.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6080/8800 (69.1%)
#> Elapsed: 56s | ETA: 25s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1600/2200 (72.7%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1613/2200 (73.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1557/2200 (70.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1534/2200 (69.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6304/8800 (71.6%)
#> Elapsed: 58s | ETA: 23s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2200 (75.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1666/2200 (75.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1596/2200 (72.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1576/2200 (71.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 6488/8800 (73.7%)
#> Elapsed: 1m | ETA: 21s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/2200 (77.3%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1716/2200 (78.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1636/2200 (74.4%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1621/2200 (73.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 6673/8800 (75.8%)
#> Elapsed: 1m 1s | ETA: 19s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1750/2200 (79.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1776/2200 (80.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1679/2200 (76.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1671/2200 (76.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 6876/8800 (78.1%)
#> Elapsed: 1m 3s | ETA: 18s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/2200 (81.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 1824/2200 (82.9%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1729/2200 (78.6%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1722/2200 (78.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 7075/8800 (80.4%)
#> Elapsed: 1m 4s | ETA: 16s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1850/2200 (84.1%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 1876/2200 (85.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1774/2200 (80.6%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1768/2200 (80.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 7268/8800 (82.6%)
#> Elapsed: 1m 6s | ETA: 14s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/2200 (86.4%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1916/2200 (87.1%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 1823/2200 (82.9%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1803/2200 (82.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7442/8800 (84.6%)
#> Elapsed: 1m 7s | ETA: 12s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1950/2200 (88.6%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1964/2200 (89.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1868/2200 (84.9%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1852/2200 (84.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7634/8800 (86.8%)
#> Elapsed: 1m 9s | ETA: 11s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 2000/2200 (90.9%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2022/2200 (91.9%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1923/2200 (87.4%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1908/2200 (86.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7853/8800 (89.2%)
#> Elapsed: 1m 11s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2050/2200 (93.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2079/2200 (94.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1963/2200 (89.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1951/2200 (88.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8043/8800 (91.4%)
#> Elapsed: 1m 12s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2100/2200 (95.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2125/2200 (96.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2012/2200 (91.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 1990/2200 (90.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 8227/8800 (93.5%)
#> Elapsed: 1m 13s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 2150/2200 (97.7%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 2163/2200 (98.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2056/2200 (93.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2040/2200 (92.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 8409/8800 (95.6%)
#> Elapsed: 1m 15s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2115/2200 (96.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2091/2200 (95.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 8606/8800 (97.8%)
#> Elapsed: 1m 17s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8800/8800 (100.0%)
#> Elapsed: 1m 19s | ETA: 0s
#> NUTS Diagnostics Summary:
#>   Total divergences:         0 
#>   Max tree depth hits:       0 
#>   Min E-BFMI across chains:  1.392 

# Posterior inclusion probabilities
summary(fit)$indicator
#>                      mean         sd        mcse n0->0 n0->1 n1->0
#> intrusion-dreams  1.00000 0.00000000          NA     0     0     0
#> intrusion-flash   1.00000 0.00000000          NA     0     0     0
#> intrusion-upset   0.88675 0.31689815 0.033097148   435    18    18
#> intrusion-physior 0.94625 0.22552370 0.023709256   206     9     9
#> dreams-flash      1.00000 0.00000000          NA     0     0     0
#> dreams-upset      0.99850 0.03870078 0.001367047     4     2     2
#> dreams-physior    0.06700 0.25002200 0.008926280  3649    82    82
#> flash-upset       0.08275 0.27550397 0.010723912  3582    86    86
#> flash-physior     1.00000 0.00000000          NA     0     0     0
#> upset-physior     1.00000 0.00000000          NA     0     0     0
#>                   n1->1     n_eff     Rhat
#> intrusion-dreams   3999        NA       NA
#> intrusion-flash    3999        NA       NA
#> intrusion-upset    3528  91.67655 1.024333
#> intrusion-physior  3775  90.47915 1.145775
#> dreams-flash       3999        NA       NA
#> dreams-upset       3991 801.44296 1.293917
#> dreams-physior      186 784.54068 1.015340
#> flash-upset         245 660.00800 1.016598
#> flash-physior      3999        NA       NA
#> upset-physior      3999        NA       NA

# Posterior pairwise effects
summary(fit)$pairwise
#>                          mean          sd        mcse     n_eff
#> intrusion-dreams  0.630011559 0.001608864 0.066123967 1689.1933
#> intrusion-flash   0.338858635 0.001426958 0.061194884 1839.1102
#> intrusion-upset   0.182650623 0.082863156 0.006962362  141.6478
#> intrusion-physior 0.197786720 0.072082734 0.005180593  193.5992
#> dreams-flash      0.499269400 0.001282161 0.059751778 2171.7848
#> dreams-upset      0.231409344 0.055056569 0.002092234  692.4647
#> dreams-physior    0.006397386 0.024184302 0.000927634  679.6946
#> flash-upset       0.008209138 0.027710743 0.001128422  603.0501
#> flash-physior     0.306644868 0.001209639 0.053094538 1926.5846
#> upset-physior     0.710281824 0.001426547 0.059325674 1729.4699
#>                        Rhat
#> intrusion-dreams  1.0031150
#> intrusion-flash   1.0009726
#> intrusion-upset   1.0076115
#> intrusion-physior 1.0122796
#> dreams-flash      1.0009953
#> dreams-upset      1.0028343
#> dreams-physior    1.0053788
#> flash-upset       1.0048742
#> flash-physior     1.0064246
#> upset-physior     0.9999516
# }
```
