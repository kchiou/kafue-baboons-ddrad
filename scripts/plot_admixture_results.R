#!/usr/bin/env Rscript

library (ggplot2)
library(car)

options(stringsAsFactors = FALSE)

mds = read.table('results/mds.txt', header = T)

# ========================================================================================
# === Analyze and plot ADMIXTURE results
# ========================================================================================

adm.files = list.files(path = "results/", pattern = "all\\.cleaned\\.LDpruned\\.[0-9]*\\.Q$")
pattern = "[0-9]+.Q"
k.vals = gsub(".Q", "", regmatches(adm.files, regexpr(pattern, adm.files)))

adm = list()
se  = list()

# Hard code colors
cols = list()
cols[["1"]] = c("#999999")
cols[["2"]] = c("#377eb8", "#ff7f00")
cols[["3"]] = c("#984ea3", "#377eb8", "#ff7f00")
cols[['4']] = c("#4daf4a", '#984ea3', '#ff7f00', '#377eb8')
cols[['5']] = c("#4daf4a", "#e41a1c", "#984ea3", "#377eb8", "#ff7f00")

for (k in k.vals) {

	k.ch = as.character(k)
	adm[[k.ch]] = read.table(paste0("results/all.cleaned.LDpruned.", k, ".Q"))

	# Swap in consistent order, based on total amount of ancestry ascribed to component
	col.order = order(colSums(adm[[k.ch]]))
	adm[[k.ch]] = data.frame(adm[[k.ch]][,col.order])
	names(adm[[k.ch]]) = paste0("V", 1:ncol(adm[[k.ch]]))
	
	group = mds$group
	ind_id = mds$ind_id
	
	num.ind = nrow(mds)
	
	new.order = order(adm[[k.ch]]$V1,
		if ("V2" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V2 } else { rep(NA, num.ind) },
		if ("V3" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V3 } else { rep(NA, num.ind) },
		if ("V4" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V4 } else { rep(NA, num.ind) },
		if ("V5" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V5 } else { rep(NA, num.ind) },
		if ("V6" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V6 } else { rep(NA, num.ind) },
		mds$group,
		mds$ind_id)
	
	# Plot animals by ancestry, largely ignoring group
	if (k == 2) {
		pdf(file='results/admixture_plot.pdf', width=15, height=5,useDingbats=FALSE)
	
			newplot = barplot(t(as.matrix(adm[[k.ch]][new.order,])), 
				col=cols[[k.ch]], 
				xlab="Individual", ylab="Ancestry", border=NA, xaxt='n', main='ADMIXTURE ancestry (k = 2)')
			text(newplot,-0.02,labels=ind_id[new.order],srt=-90,adj=c(0,0.5),xpd=NA,cex=0.5)
		
		dev.off()
	}
	
	# Plot animals by group first, then ancestry within
	new.order.group = order(mds$group,
		adm[[k.ch]]$V1,
		if ("V2" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V2 } else { rep(NA, num.ind) },
		if ("V3" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V3 } else { rep(NA, num.ind) },
		if ("V4" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V4 } else { rep(NA, num.ind) },
		if ("V5" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V5 } else { rep(NA, num.ind) },
		if ("V6" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V6 } else { rep(NA, num.ind) },
		mds$ind_id)

	groups = mds[new.order.group,]$group
	group.starts = which(groups != c(groups[-1], NA))

	# Read in file with SE estimates
	se[[k.ch]] = read.table(paste0("results/all.cleaned.LDpruned.", k, ".Q_se"))
	
	pdf(file=paste0('results/admixture_plot_K', k, '_by_group.pdf'), width=15, height=5,useDingbats=FALSE)
	
		mp = barplot(t(as.matrix(adm[[k.ch]][new.order.group,])), 
			col=cols[[k.ch]], 
			xlab="Group", ylab="Ancestry", border=NA, xaxt='n', main=paste0('ADMIXTURE ancestry (k = ',k,')'))
	
		bar.width = mp[2] - mp[1]
		border.locs = mp[group.starts] + (0.5 * bar.width)
		abline(v=border.locs, lwd=3)

		group.mps = colMeans(rbind(	c(0, border.locs), 
									c(border.locs, max(mp) + 0.5*bar.width)))

		uniq.groups = unique(groups)
	
		mtext(uniq.groups, side=1, at=group.mps, padj=1)

		if (k == 2) {
			arrows(	mp, adm[[k.ch]][new.order.group,]$V1+se[[k.ch]][new.order.group,]$V1, 
					mp, adm[[k.ch]][new.order.group,]$V1, 
					angle=90, code=1, length=0)
			arrows(	mp, adm[[k.ch]][new.order.group,]$V1-se[[k.ch]][new.order.group,]$V1, 
					mp, adm[[k.ch]][new.order.group,]$V1, 
					angle=90, code=1, length=0)
		}
	dev.off()

}

# ========================================================================================
# === Group results from various runs (k: 2-5)
# ========================================================================================

pdf(file=paste0('results/admixture_plot_K2-',k.vals[length(k.vals)],'.pdf'), width=15, height=10,useDingbats=FALSE)

	layout(matrix(1:(length(k.vals)-1),ncol=1))

	# Order individuals by group, then ancestry (k=2), then individual
	fixed.order = order(mds$group,adm[[2]]$V1,adm[[2]]$V2,mds$ind_id)

	for (k in k.vals[2:length(k.vals)]) {

		k.ch = as.character(k)

		lab1 = if (k == k.vals[length(k.vals)]) 'Group' else NULL
		lab2 = NULL

		groups = mds[fixed.order,]$group
		group.starts = which(groups != c(groups[-1], NA))	
			mp = barplot(t(as.matrix(adm[[k.ch]][fixed.order,])), 
				col=cols[[k.ch]], 
				xlab=lab1, ylab=lab2, border=NA, yaxt='n', xaxt='n',cex.lab=2,main=paste0('ADMIXTURE ancestry (k = ',k,')'),cex.main=2)
	
			bar.width = mp[2] - mp[1]
			border.locs = mp[group.starts] + (0.5 * bar.width)
			abline(v=border.locs, lwd=3)

			group.mps = colMeans(rbind(	c(0, border.locs), 
										c(border.locs, max(mp) + 0.5*bar.width)))

			uniq.groups = unique(groups)
	
			if (k == k.vals[length(k.vals)]) mtext(uniq.groups, side=1, at=group.mps, padj=1)

	}

dev.off()

# ========================================================================================
# === Plot various ancestry inferences against one another
# ========================================================================================

# Set adm equal to k=2 run, just for simplicity's sake
adm = read.table(paste0("results/all.cleaned.LDpruned.2.Q"))

# Swap ADM1 and ADM2 if they are not the order I arbitrarily like them in.
if (adm[which(mds$ind_id == "BZ11-001"),]$V1 > 0.5) {
	adm = data.frame(V1=adm$V2, V2=adm$V1)
}

names(adm) = paste0("ADM", 1:2)
mds.adm = cbind(mds, adm)

# Source: http://goo.gl/K4yh
lm_eqn = function(df){
    m = lm(y ~ x, df);
    eq = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
            list(   a = format(coef(m)[1], digits = 2), 
                    b = format(coef(m)[2], digits = 2), 
                    r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

# Create negligible variation in ancestry indices (to avoid errors in quantile calculation)
mds.adm$adm = mds.adm$ADM1 + runif(nrow(mds.adm),-1e-10,1e-10)

# ----------------------------------------------------------------------------------------
# --- Create violin / boxplot of ADMIXTURE inferred ancestry by taxon
# ----------------------------------------------------------------------------------------

mds.adm$taxon.f = factor(mds.adm$taxon, levels = c("Kinda","Hybrid","Grayfoot"))

p = ggplot(aes(y=ADM1, x=taxon.f), data=mds.adm, pch=16) + 
	geom_boxplot() + geom_jitter(height=0,alpha=0.25) + theme_bw() + 
	xlab("Phenotypic Taxon") + ylab("ADMIXTURE Inferred Ancestry") +
	ggtitle("Ancestry by taxon (ADMIXTURE)")

ggsave(p, file='results/boxplot_ancestry_by_taxon_adm.pdf',width=6,height=6,useDingbats=FALSE)

p = ggplot(aes(y=adm, x=taxon.f), data=mds.adm, pch=16) + 
	geom_violin(scale='width',draw_quantiles=0.5,trim=TRUE) +
	geom_jitter(height=0,alpha=0.25) + theme_bw() + 
	xlab("Phenotypic Taxon") + ylab("ADMIXTURE Inferred Ancestry") + ggtitle("Ancestry by taxon (ADMIXTURE)")

ggsave(p, file='results/violin_ancestry_by_taxon_adm.pdf',width=6,height=6,useDingbats=FALSE)

# ----------------------------------------------------------------------------------------
# --- Create violin / boxplot of ADMIXTURE inferred ancestry by group
# ----------------------------------------------------------------------------------------

p = ggplot(aes(y=ADM1, x=group), data=mds.adm, pch=16) + 
	geom_boxplot() + geom_jitter(height=0,alpha=0.25) + theme_bw() + 
	xlab("Group") + ylab("ADMIXTURE Inferred Ancestry") +
	ggtitle("Ancestry by group (ADMIXTURE)")

ggsave(p, file='results/boxplot_ancestry_by_group_adm.pdf',width=6,height=3)

p = ggplot(aes(y=adm, x=group), data=mds.adm, pch=16) + 
	geom_violin(scale='width',draw_quantiles=0.5,trim=TRUE) +
	geom_jitter(height=0,alpha=0.25) + theme_bw() +
	xlab("Group") + ylab("ADMIXTURE inferred ancestry") + # + ggtitle("Ancestry by group (ADMIXTURE)")
	theme_classic()

ggsave(p, file='results/violin_ancestry_by_group_adm.pdf',width=6,height=3,useDingbats=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot first dimension of IBS MDS as function of ADMIXTURE inferred ancestry
# ----------------------------------------------------------------------------------------

data.label = data.frame(x=1/3, y=0.9, label = c(lm_eqn(data.frame(x=mds.adm$ADM1, y=mds.adm$C1.norm))))

p = ggplot(aes(ADM1, C1.norm), data=mds.adm, pch=16) + 
	xlim(0,1) + ylim(0,1) + 
	xlab("ADMIXTURE Inferred Ancestry") + ylab("MDS Dimension 1 (Scaled)") + 
	geom_point() + 
	geom_smooth(method=lm, se=FALSE) + theme_bw() +
	geom_text(data = data.label, aes(x=x, y=y, label=label), size=3, parse=TRUE) +
	coord_fixed() +
	theme_classic()

lm.fit = lm(mds.adm$ADM1 ~ mds.adm$C1.norm)
#summary(lm.fit)

ggsave(p, file='results/ancestries_mds_vs_adm.pdf',width=4,height=4,useDingbats=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot spatial map of ADMIXTURE ancestry
# ----------------------------------------------------------------------------------------

# Bring in Flora Jay's POPS script (http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R)
source('scripts/POPSutilities.R')

# show.key function (http://membres-timc.imag.fr/Olivier.Francois/TESS_Plot.html)
show.key = function(cluster=1,colorGradientsList=lColorGradients){
        ncolors=length(colorGradientsList[[cluster]])
        barplot(matrix(rep(1/10,10)),col=colorGradientsList[[cluster]][(ncolors-9):ncolors],main=paste("Cluster",cluster))}

coord = read.table('data/all.cleaned.LDpruned.coord')

pdf('results/adm_geo_plot_zambia.pdf',width=6,height=5,useDingbats=FALSE)
	layout(matrix(c(rep(1,4),2,3),2,3,byrow=FALSE), widths=c(3.8,0.2), respect=FALSE)
	par(ps=22,mar=c(5,5,5,2))
	asc.raster = 'data/zambia.asc'
	grid=createGridFromAsciiRaster(asc.raster)
	constraints=getConstraintsFromAsciiRaster(asc.raster,cell_value_min=0)
	maps(matrix=adm, coord, grid, constraints, method="max", xlab="Longitude", ylab="Latitude")
	par(ps=16,mar=c(5,2,5,5))
	for(k in 1:2) show.key(k)
dev.off()

taxon=c('Grayfoot','Kinda')

show.key=function(cluster=1,colorGradientsList=lColorGradients){
        ncolors=length(colorGradientsList[[cluster]])
        barplot(matrix(rep(1/10,10)),col=colorGradientsList[[cluster]][(ncolors-9):ncolors],main=paste0(taxon[cluster],'\nancestry'),cex.main=12/13,cex.axis=10/11)}

pdf('results/adm_geo_plot_studyarea.pdf',width=6,height=4,useDingbats=FALSE)
	layout(matrix(c(rep(1,4),2,3),2,3,byrow=FALSE), widths=c(3.8,0.2), respect=FALSE)
	par(ps=22,mar=c(5,5,2,2),las=1)
	asc.raster = 'data/studyarea.asc'
	grid=createGridFromAsciiRaster(asc.raster)
	constraints=getConstraintsFromAsciiRaster(asc.raster,cell_value_min=0)
	maps(matrix=adm, coord, grid, constraints, method="max", xlab="Longitude", ylab="Latitude",axes=FALSE)
	axis(side=1,at=25:29,cex.axis=0.75)
	axis(side=2,at=-17:-14,cex.axis=0.75)
	box()
	par(ps=16,mar=c(5,2,5,5))
	for(k in 1:2) show.key(k)
dev.off()

pdf('results/adm_geo_plot_kafue.pdf',width=5,height=5,useDingbats=FALSE)
	layout(matrix(c(rep(1,4),2,3),2,3,byrow=FALSE), widths=c(1.2,0.8), respect=FALSE)
	par(ps=22,mar=c(5,5,5,2))
	asc.raster = 'data/kafue.asc'
	grid=createGridFromAsciiRaster(asc.raster)
	constraints=getConstraintsFromAsciiRaster(asc.raster,cell_value_min=0)
	maps(matrix=adm, coord, grid, constraints, method="max", xlab="Longitude", ylab="Latitude")
	par(ps=16,mar=c(5,2,5,5))
	for(k in 1:2) show.key(k)
dev.off()

# ----------------------------------------------------------------------------------------
# --- Plot IBS MDS plot again, this time coloring by ADMIXTURE-inferrred ancestry
# ----------------------------------------------------------------------------------------

mds.adm$taxon.adm = 'Hybrid'
mds.adm$taxon.adm[mds.adm$ADM1 > 0.999] = 'Grayfoot'
mds.adm$taxon.adm[mds.adm$ADM1 < 0.001] = 'Kinda'

key = c("Chunga HQ", "Mwengwa Rp", "Namiyezhi R", "Malala Cmp", "Musa Br", "Top Musa", "N Nkala Rd", "Ngoma Airstr", "Dendro P", "Nanzhila Pl", "Choma", "L Zambezi NP")
names(key) = LETTERS[1:12]
key = key[names(key) %in% mds$group]

names(key) = LETTERS[1:length(key)]
mds.adm$group = factor(mds.adm$group)
levels(mds.adm$group) = LETTERS[1:nlevels(mds.adm$group)]

ggplot(mds.adm,aes(C1,C2,shape=group,color=taxon.adm)) +
	geom_point(cex=3) +
	scale_shape_manual(name='Group',values=names(key),labels=key) +
	scale_color_manual(name='Taxon',values=c('#e41a1c','#4daf4a','#377eb8')) +
	coord_fixed() +
	theme_classic() +
	theme(legend.position='bottom',panel.grid=element_blank()) +
	xlab('Dimension 1') + ylab('Dimension 2')
ggsave(filename='results/mds_plot.pdf',width=6,height=6,useDingbats=FALSE)
