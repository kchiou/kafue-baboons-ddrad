#!/bin/bash

# ------------------------------------------------------------------------------
# --- Aggregate bgc results
# ------------------------------------------------------------------------------

# Get combined log-likelihood
paste -d ',' results/bgc_stat_ln1_* > results/bgc_stat_ln1

# Get combined alpha and beta
paste -d ',' results/bgc_stat_a0_* > results/bgc_stat_a0
paste -d ',' results/bgc_stat_b0_* > results/bgc_stat_b0

# Get combined hybrid index
paste -d ',' results/bgc_stat_hi_* > results/bgc_stat_hi

exit