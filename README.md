# DynamicGGM
Code used for the paper "Change point detection in dynamic Gaussian graphical models: the impact of COVID-19 pandemic on the US stock market"  by B. Franzolini, A. Beskos, M. De Iorio, W. Poklewski Koziell and K. Grzeszkiewicz  arXiv:2208.00952v2 , 2023

## main_to_run.R 
--> R code to run to reproduce results and simulate data.

## DGGM.R 
--> R functions to estimate the model (needed to run main_to_run.R).

details on functions in DGGM: 
- compute_prior: returns the inverse of the normalizing constant of uniform prior for each possible number of change points
- detect_free_points: given a sequence of change points find available change points
- detect_free_points_btw2: given a sequence of change points find available change points between two change points
- ess: compute the effective sample size
- lik_block: compute the likelihood of a block given the corresponding graph
- mutate_all_G: perform a MH step for the whole sequence of graphs and all particles to sample for their posterior (not used in the main code)
- mutate_G: perform MH step for one graph and all particles to sample for their temperated-posterior
- particle_filter: inner component of the algorithm, compute the marg lik and sample the graphs
- sample_ChangePoints: outer component, MH for change points
- sim.data: simulate data from the model (not used in the main code)
- simulate_data: simulate data from scenarios described in Franzolini et al. (2023)
- sim.G: sample a graph at time 0 from the prior
- sim.N.G: sample N graphs at time 0 from the prior
- sim.GG: sample a graph from the prior given the previous graph
- sim.N.GG: sample N graphs from the prior given N previous graphs
- temperatures_tuning: compute temperatures adaptively

## comparisonGFGL.m
--> MatLab script to perform comparisons with group-fused graphical lasso (GFGL) of Gibberd and Nelson (2017)

Ref: 
Gibberd, A. J. and J. D. Nelson (2017). Regularized estimation of piecewise constant Gaussian graphical models:
The group-fused graphical lasso. Journal of Computational and Graphical Statistics 26, 623–634.
