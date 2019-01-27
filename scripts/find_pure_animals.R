#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

mds = read.table('results/mds.txt', header = T)

# Read ancestry results from ADMIXTURE analysis with K = 2
adm = read.table("results/all.cleaned.LDpruned.2.Q")

# ========================================================================================
# === Make lists of pure Kinda and chacma animals
# ========================================================================================

# Give the ADMIXTURE results descriptive names
names(adm) = c("ADMIXTURE_1", "ADMIXTURE_2")

# Add the ADMIXTURE results to the MDS results (with extra info)
mds.adm = cbind(mds, adm)

# Find extremes
pure_high_adm = mds.adm[mds.adm$ADMIXTURE_1 < 0.001,]
pure_low_adm  = mds.adm[mds.adm$ADMIXTURE_1 > 0.999,]

# Find hybrids
hybrids_adm = mds.adm[(mds.adm$ADMIXTURE_1 >= 0.001) & (mds.adm$ADMIXTURE_1 <= 0.999),]

# Figure out which represents the pure chacma and which the pure Kinda
# since the ADMIXTURE results could be inverted

num.high.kind = length(pure_high_adm$taxon[pure_high_adm$taxon == "Kinda"])
num.high.chac = length(pure_high_adm$taxon[pure_high_adm$taxon == "Grayfoot"])

if (num.high.kind > num.high.chac) {
	# Subset with high ADMIXTURE values is predominantly Hamadryas
	pure.kind.subset = pure_high_adm
	pure.chac.subset = pure_low_adm
} else {
	# Subset with high ADMIXTURE values is predominantly Anubis
	pure.kind.subset = pure_low_adm
	pure.chac.subset = pure_high_adm
}

# Save these datasets as text files
write.table(pure.kind.subset, file="results/pure_kind_data.txt")
write.table(pure.chac.subset, file="results/pure_chac_data.txt")
write.table(hybrids_adm,      file="results/hybrid_data.txt")

# Also output first two columns (FID and IID, which are currently the same thing)
# for plink to read.

write.table(pure.kind.subset[,1:2], file="results/pure_kind_list.txt", 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(pure.chac.subset[,1:2], file="results/pure_chac_list.txt", 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(hybrids_adm[,1:2],      file="results/hybrid_list.txt", 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)