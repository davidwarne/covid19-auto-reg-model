# Bayesian analysis of the COVID-19 pandemic

Bayesian analysis of COVID-19 case time-series for individual regions using a simulation based model and Approximate Bayesian Computation (ABC) using adaptive sequential Monte Carlo. 

To accompany the paper DJ Warne, A Ebert, D Drovandi, W Hu, A Mira, K Mengersen (2020) ``Hindsight is 2020 vision: a characterisation of the global response to the COVID-19 pandemic'' (link TBA)

## Authors

1. Christopher Drovandi (c.drovandi@qut.edu.au),
                School of Mathematical Sciences, 
                Science and Engineering Faculty, 
                Queensland University of Technology 
Google Scholar: (https://scholar.google.com.au/citations?user=jc0284YAAAAJ&hl=en)

2. David J. Warne (david.warne@qut.edu.au),
                School of Mathematical Sciences, 
                Science and Engineering Faculty, 
                Queensland University of Technology 
Google Scholar: (https://scholar.google.com.au/citations?user=t8l-kuoAAAAJ&hl=en)

## Model summary
The model is a discrete state continuous-time Markov Process model that is based on an SIR model with important extentions. The form of the interactions are:
$$S \overset{g(A_t,R_t,D_t)I/P}{\rightarrow} I,\quad A \overset{\beta}{\rightarrow} R,\quad I\overset{\gamma}{\rightarrow} A, \quad A\overset{\delta}{\rightarrow} D,\quad\text{and}\quad I\overset{\eta\beta}{\rightarrow} R^u$$
where $S_t$ is susceptible (latent), $I_t$ is infected (latent), $A_t$ are active confirmed cases, $R_t$ are case recoveries, $D_t$ are case fatalities, $C_t = A_t + R_t + D_t$ are the cumulative confirmed cases, and $R_t^u$ are removed (recoveries of deaths of infected latent variable) at time $t > 0$. Parameters are the infection rates $\alpha_0 > 0$ and $\alpha > 0$, case recovery rate $\beta > 0$, case fatality rate $\delta > 0$, case detection rate $\gamma > 0$, and overal removal rate relative to case recovery $\eta > 0$. The function $g(A_t,R_t,D_t)$ represents a regulatory response of the community (e.g., social distancing etc...) we assume the form $g(A_t,R_t,D_t) = \alpha_0 + \dfrac{\alpha}{1+U(A_t,R_t,D_t)^{n}}$ where $n \geq 0$ and $U(A_t,R_t,D_t)$ is a utility function of observables. Inital contion for $I$ is given by $I_0 = \kappa A_0$ where $\kappa \geq 0$. 

## Functions and scripts

The following files are provided:

* `run_smc_intensity_reg.m` runs analysis pipeline (computes posteriors and samples prior predictive) given a country ID (ISO-3166 alpha3 code).

* `run_smc_intensity_reg_cluster.m` and `submit_all.sh` runs analysis pipeline in parallel on a HPC cluster for  all 252 countries.

* `plot_results.m` plot posterior distributions and point estimate scatter plots using results from HPC cluster run.

* `simuldata_reg.m` forwards stochastic simulation of model given a vector of parameters. Uses the utility function $U(A_t,R_t,D_t) = A_t + R_t + D_t = C_t$. 

* `simuldata_reg_fD.m` forwards stochastic simulation of model given a vector of parameters. Uses the utility function $U(A_t,R_t,D_t) = D_t$. 

* `smry.m` summary statistic function for usage in the ABC method.

* `smc_abc_rw.m` Adaptive sequential Monte Carlo sampler for ABC. This should not need to be modified even if the model completely changes.

* `smc_abc_rw_par.m` multithreaded version Adaptive sequential Monte Carlo sampler for ABC. This should not need to be modified even if the model completely changes.

* `TauLeapingMethod.m` routine for approximate stochastic simulation of discrete-state continuous-time Markov process (used in model simulation).


## Usage

1. Edit `run_smc_intensity_reg.m` to ensure  the variable `DATA_DIR` is pointing to a directory that contains clones of the two github repos `COVID19data` and `COVID-19_ISO-3166`.
2. Define a country/region id in the MATLAB workspace (i.e., the row number in the population table for the region of interest), e.g., for China
use `iso_id = 'CHN'` or for Italy use `iso_id = 'ITA'`.
3. If required, edit the SMC parameters, defaults are quite reasonable for the moment. Please contact the authors if there are any difficulties editing these.
4. Run the main script `run_smc_intensity_reg`. This may take a while to run. Currently, the code is not parallelised to ensure reproducibility (RNG in parallel are not reproducible). If Parallel SMC is desired, then edit `run_smc_intensity_reg.m` to utilise `smc_abc_rw_par.m`.
5. To run on a HPC cluster use the script `submit_all.sh` shell script. WARNING: due to differences in HPC cluster management across institutions, this script will require modification.
