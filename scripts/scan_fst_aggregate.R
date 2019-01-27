#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

system('cat checkpoints/p.booster.[0-9]*.txt > checkpoints/p.booster.txt')
p.booster = read.table('checkpoints/p.booster.txt')

checkpoint.files = list.files('checkpoints')
checkpoint.files = checkpoint.files[grep('^pval.perm',checkpoint.files)]

# Calculate checkpoint file with max permutations
max.permutations = checkpoint.files[which.max(gsub('^pval\\.perm\\.1e([0-9]+)\\.RData','\\1',checkpoint.files))]

load(max.permutations)

names(p.booster) = c('i','nam','fst','pval','nreps')

for (gene in unique(p.booster$nam)) {
	p.pre = real.fst$pval[real.fst$nam %in% gene]
	p.post = p.booster$pval[p.booster$nam %in% gene]
	reps.pre = real.fst$nreps[real.fst$nam %in% gene]
	real.fst$nreps[real.fst$nam %in% gene] = reps.pre * 10
	real.fst$pval[real.fst$nam %in% gene] = mean(c(p.pre,p.post))
}

# Write final scan results
write.csv(real.fst,file='results/scan_selection_fst.csv',row.names=FALSE)
