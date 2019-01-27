#!/bin/bash

# ------------------------------------------------------------------------------
# --- Run ADMIXTURE in cross-validation mode
# ------------------------------------------------------------------------------

for K in {1..10}; do 
	admixture --cv results/all.cleaned.LDpruned.bed ${K} | tee reports/ADMIXTURE_log${K}.out
done

grep -h CV `ls -v reports/ADMIXTURE_log[0-9]*` > reports/ADMIXTURE_log_cv.out

# Get rid of these results, since we later run ADMIXTURE with chosen K
rm all.cleaned.LDpruned.[0-9]*.P
rm all.cleaned.LDpruned.[0-9]*.Q

exit;