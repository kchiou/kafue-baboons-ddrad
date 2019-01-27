#!/bin/bash

# ------------------------------------------------------------------------------
# --- Run bgc
# ------------------------------------------------------------------------------

# Arguments set
# -O 0      : use HDF5 output (requires program estpost)
# -x 400000 : 400,000 MCMC generations
# -n 200000 : discard 50% (200,000 steps) as burn-in
# -t 40     : sample MCMC chain every 40 generations
# -N 1      : use read counts
# -E 0.0001 : 99.99% accuracy (error rate corresponds to Phred score >= 40)
# -q 1      : calculate cline parameter quantiles
# -I 1      : do not initialize hybrid index/ancestry based on data
# -p 1      : print default parameters + precision parameters

# Defaults kept
# -d 1      : all loci are diploid
# -o 1      : do not assume a constant population-level cline parameter variance for all loci
# -i 0      : do not print interspecific heterozygosity
# -s 1      : sum-to-zero constraint on locus cline parameters
# -T 0      : use a full gamma prior
# -u 0.1    : maximum deviation from uniform for proposed hybrid index
# -g 0.05   : standard deviation for Gaussian proposal of cline parameter gamma
# -z 0.05   : standard deviation for Gaussian proposal of cline parameter zeta
# -e 0.02   : standard deviation for Gaussian proposal of cline parameters eta and kappa

bgc -a data/bgc_p0in.txt -b data/bgc_p1in.txt -h data/bgc_admixedin.txt \
	-O 0 -x 400000 -n 200000 -t 40 -N 1 -E 0.0001 -q 1 -I 1 -p 1 \
	-F results/bgc_mcmcout_${1}

# ------------------------------------------------------------------------------
# --- Run estpost
# ------------------------------------------------------------------------------

# Arguments set
# -p LnL    : summarize log-likelihood; other options: alpha beta eta eta-quantile
#           :      gamma-quantile gamma-quantile-local hi interspecific-het kappa
#           :      kappa-quantile rho tau-alpha tau-beta zeta-quantile zeta-quantile-local
# -s 2      : Convert summary to plain text; other options: 0 = posterior estimates and 
#           :      credibility intervals, 1 = histogram of posterior samples
# -w 1      : Write parameter identification and headers to file

# Defaults kept
# -c 0.95   : 95% credibility interval
# -b 0      : Discard 0 generations as burn-in (based on number of thinned samples)
# -h 20     : number of bins for posterior sample histogram

# Get log-likelihood
estpost -i results/bgc_mcmcout_${1}.hdf5 -p LnL -o results/bgc_stat_ln1_${1} -s 2 -w 0

# Get alpha and beta
estpost -i results/bgc_mcmcout_${1}.hdf5 -p alpha -o results/bgc_stat_a0_${1} -s 2 -w 0
estpost -i results/bgc_mcmcout_${1}.hdf5 -p beta -o results/bgc_stat_b0_${1} -s 2 -w 0

# Get hybrid index
estpost -i results/bgc_mcmcout_${1}.hdf5 -p hi -o results/bgc_stat_hi_${1} -s 2 -w 0

exit