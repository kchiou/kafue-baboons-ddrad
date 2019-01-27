#/bin/bash

# ------------------------------------------------------------------------------
# --- Run ADMIXTURE for a specific k
# ------------------------------------------------------------------------------

ADMIX_K=$1

admixture -B200 results/all.cleaned.LDpruned.bed ${ADMIX_K}
mv all.cleaned.LDpruned.${ADMIX_K}.Q results/
mv all.cleaned.LDpruned.${ADMIX_K}.Q_se results/
mv all.cleaned.LDpruned.${ADMIX_K}.Q_bias results/
mv all.cleaned.LDpruned.${ADMIX_K}.P results/

exit