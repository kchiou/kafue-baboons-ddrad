#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(parallel)

options(stringsAsFactors=FALSE)

# Source: http://goo.gl/K4yh
lm_eqn = function(df){
    m = lm(y ~ x, df);
    eq = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
            list(   a = format(coef(m)[1], digits = 2), 
                    b = format(coef(m)[2], digits = 2), 
                    r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

# Plot log-likelihoods, reported every 1000 generations

lnl = as.numeric(system('cat reports/bgc_err.txt | grep "^mcmc iteration" | sed -e "s/.*LnL: //g"',intern=TRUE))

lnl = data.frame(gen=(1:length(lnl) - 1) * 1000,lnl=lnl)

ggplot(lnl,aes(gen/1e6,lnl)) + geom_point(pch=21) + theme_classic() + xlab('MCMC generation (millions)') + ylab('log likelihood')
ggsave(file='results/bgc_convergence.pdf',width=6,height=4,useDingbats=FALSE)

ln1 = round(scan(file='results/bgc_stat_ln1',sep=','))

library(coda)
mcmc.out = mcmc(ln1)
h.out = heidel.diag(mcmc.out,0.1)

if (h.out['var1','pvalue'] < 0.05) stop('Failed stationarity test')
if (abs(h.out['var1','halfwidth'] / h.out['var1','mean']) >= 0.1) stop('Failed half-width test')

ln1 = data.frame(samp=1:length(ln1),ln1=ln1)
ggplot(ln1,aes(samp,ln1)) + geom_point(pch=21) + theme_classic() + xlab('MCMC sample') + ylab('log likelihood')
ggsave(file='results/bgc_log_likelihood.pdf',width=6,height=4,useDingbats=FALSE)

# Plot hybrid indices

hi = read.csv('results/bgc_stat_hi')
hi.mean = rowMeans(as.matrix(hi))

mds = read.table('results/mds.txt', header = T)

ind_info = read.csv('data/individual_info.csv')
ind_info = ind_info[ind_info$Individual.ID %in% mds$IID,]
ind_info$pop = NA

pure1 = gsub('\\.PE$','',system(paste('cut -f 1','results/pure_kind_list.txt'),intern=TRUE))
pure2 = gsub('\\.PE$','',system(paste('cut -f 1','results/pure_chac_list.txt'),intern=TRUE))

ind_info$pop[ind_info$Individual.ID %in% pure1] = 'pure 1'
ind_info$pop[ind_info$Individual.ID %in% pure2] = 'pure 2'
# ind_info$pop[is.na(ind_info$pop)] = paste('pop',as.numeric(factor(ind_info$Group[is.na(ind_info$pop)])) - 1)
ind_info$pop[is.na(ind_info$pop)] = 'pop 0'

ind_info = ind_info[order(ind_info$Group,ind_info$Individual.ID),c('Individual.ID','pop','Group','Locality')]
names(ind_info) = c('id','pop','group','locality')
row.names(ind_info) = NULL

# Bring in admixture results
adm = read.table('results/all.cleaned.LDpruned.2.Q')

# Swap ADM1 and ADM2 if they are not the order I arbitrarily like them in.
if (adm[which(mds$ind_id == "BZ11-001"),]$V1 > 0.5) {
	adm = data.frame(V1=adm$V2, V2=adm$V1)
}

names(adm) = paste0("ADM", 1:2)
mds.adm = cbind(mds, adm)
rownames(mds.adm) = mds.adm$IID

hybrids = ind_info[ind_info$pop %in% 'pop 0',]$id

hybrid.indices = mds.adm[hybrids,]
rownames(hybrid.indices) = NULL

hybrid.indices$BGC = hi.mean

data.label = data.frame(x=1/3, y=0.9, label = c(lm_eqn(data.frame(x=hybrid.indices$ADM1, y=hybrid.indices$BGC))))

p = ggplot(aes(ADM1, BGC), data=hybrid.indices, pch=16) + 
	xlim(0,1) + ylim(0,1) + 
	xlab("ADMIXTURE hybrid index") + ylab("bgc hybrid index") + 
	geom_point() + 
	geom_smooth(method=lm, se=FALSE) +
	geom_text(data = data.label, aes(x=x, y=y, label=label), size=3, parse=TRUE) +
	coord_fixed() +
	theme_classic()

lm.fit = lm(hybrid.indices$ADM1 ~ hybrid.indices$BGC)
summary(lm.fit)

ggsave(p, file='results/ancestries_adm_vs_bgc.pdf',width=4,height=4,useDingbats=FALSE)

a.out = read.csv('results/bgc_stat_a.out')
b.out = read.csv('results/bgc_stat_b.out')
qa.out = read.csv('results/bgc_stat_qa.out')
qb.out = read.csv('results/bgc_stat_qb.out')

names(a.out)[4:5] = c('lb','ub')
names(b.out)[4:5] = c('lb','ub')
names(qa.out)[4:5] = c('lb','ub')
names(qb.out)[4:5] = c('lb','ub')

# Excess ancestry compared to genome-wide average introgression

a.out$excess = NA
a.out$excess[a.out$lb > 0] = 'pos'
a.out$excess[a.out$ub < 0] = 'neg'

b.out$excess = NA
b.out$excess[b.out$lb > 0] = 'pos'
b.out$excess[b.out$ub < 0] = 'neg'

# Add in quantiles
a.out$q = qa.out$mean
b.out$q = qb.out$mean

# qnorm takes SD, so take the square root of the quantile estimate
a.out$qlb = qnorm(0.025,0,sqrt(a.out$q))
a.out$qub = qnorm(0.975,0,sqrt(a.out$q))

b.out$qlb = qnorm(0.025,0,sqrt(b.out$q))
b.out$qub = qnorm(0.975,0,sqrt(b.out$q))

a.out$outlier = NA
a.out$outlier[a.out$median > a.out$qub] = 'pos'
a.out$outlier[a.out$median < a.out$qlb] = 'neg'

b.out$outlier = NA
b.out$outlier[b.out$median > b.out$qub] = 'pos'
b.out$outlier[b.out$median < b.out$qlb] = 'neg'

snps = read.delim('results/bgc_snp_info.txt')

snps$alpha = a.out$mean
snps$beta = b.out$mean

snps$alpha.excess = a.out$excess
snps$beta.excess = b.out$excess
snps$alpha.outlier = a.out$outlier
snps$beta.outlier = b.out$outlier

snps$crazy.a = !is.na(snps$alpha.excess) & !is.na(snps$alpha.outlier)
snps$crazy.b = !is.na(snps$beta.excess) & !is.na(snps$beta.outlier)

# Excess Kinda ancestry and fast introgression rate
sum((snps$alpha.excess %in% 'pos' | snps$alpha.outlier %in% 'pos') & (snps$beta.excess %in% 'pos' | snps$beta.outlier %in% 'pos'))

# Excess Kinda ancestry and slow introgression rate
sum((snps$alpha.excess %in% 'pos' | snps$alpha.outlier %in% 'pos') & (snps$beta.excess %in% 'neg' | snps$beta.outlier %in% 'neg'))

# Excess grayfoot ancestry and fast introgression rate
sum((snps$alpha.excess %in% 'neg' | snps$alpha.outlier %in% 'neg') & (snps$beta.excess %in% 'pos' | snps$beta.outlier %in% 'pos'))

# Excess grayfoot ancestry and slow introgression rate
sum((snps$alpha.excess %in% 'neg' | snps$alpha.outlier %in% 'neg') & (snps$beta.excess %in% 'neg' | snps$beta.outlier %in% 'neg'))

write.csv(snps,file='results/snps_bgc.csv')
