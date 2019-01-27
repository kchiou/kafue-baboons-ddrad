#!/usr/bin/env Rscript

library(fields)

library(maptools)
library(raster)

library(ggplot2)
library(ggtern)

source('scripts/POPSutilities.R')

# Read in data

asc.raster = 'data/studyarea.asc'

study.area=createGridFromAsciiRaster(asc.raster)

cliff = read.delim('data/sex_chr_genotypes.txt')
mt = cliff[!is.na(cliff$mt),]
ychr = cliff[!is.na(cliff$ychr),]
mt.coord = mt[,c('x','y')]
ychr.coord = ychr[,c('x','y')]

mt = data.frame(mt1=mt$mt,mt2=1-mt$mt)
ychr = data.frame(ychr1=ychr$ychr,ychr2=1-ychr$ychr)

# Set adm equal to k=2 run, just for simplicity's sake
adm = read.table(paste0("results/all.cleaned.LDpruned.2.Q"))
mds = read.table('results/mds.txt', header = T)


# Swap ADM1 and ADM2 if they are not the order I arbitrarily like them in.
if (adm[which(mds$ind_id == "BZ11-001"),]$V1 > 0.5) {
	adm = data.frame(V1=adm$V2, V2=adm$V1)
}

names(adm) = paste0("ADM", 1:2)

coord = read.table('data/all.cleaned.LDpruned.coord')
names(coord) = c('x','y')

adm.predictions = predict(Krig(coord, adm[[1]], theta = 10),study.area)
mt.predictions = predict(Krig(mt.coord, mt[[1]], theta = 10),study.area)
ychr.predictions = predict(Krig(ychr.coord, ychr[[1]], theta = 10),study.area)

# Now use sample localities to select values

# Because cells are defined by their center points, converting coordinates is a matter of scanning for the nearest x and y
# Write a function that does this
# Function takes a matrix with two columns of coordinates and figures out which x and y indices of the grid matrix correspond

get.indices = function(coord,predictions) {
	as.numeric(apply(coord,1,function(r) {
		x = r[1]
		y = r[2]
		xlen = length(predictions$x)
		ylen = length(predictions$y)
		x.ind = which.min(abs(x - predictions$x))
		y.ind = which.min(abs(y - predictions$y))
		(y.ind-1) * xlen + x.ind
#		c(x=x.ind,y=y.ind)
	}))
}

# It works!!!

unique.coords = unique(rbind(coord,mt.coord,ychr.coord))
rownames(unique.coords) = NULL

# The second argument to get.indices does not matter which prediction set it comes from as only the geographic coordinates are needed
adm.interpolations = adm.predictions[get.indices(unique.coords,as.surface(study.area,adm.predictions)),]
mt.interpolations = mt.predictions[get.indices(unique.coords,as.surface(study.area,adm.predictions)),]
ychr.interpolations = ychr.predictions[get.indices(unique.coords,as.surface(study.area,adm.predictions)),]

real.interpolations = data.frame(n=adm.interpolations,m=mt.interpolations,y=ychr.interpolations,type='real')

num.cells = nrow(adm.predictions)

# Cut off the eastern half to avoid extrapolation
# This arranges the indices into a matrix in which it is the eastern half we wish to crop off
# Technically, the matrix is flipped heightwise but that is not important here

fake.matrix = t(matrix(1:num.cells,nrow=504))
kafue.indices = sort(as.numeric(fake.matrix[,1:round(ncol(fake.matrix)/2)]))

set.seed(141251)
sim.indices = sample(kafue.indices,40,replace=FALSE)

adm.simulations = adm.predictions[sim.indices,]
mt.simulations = mt.predictions[sim.indices,]
ychr.simulations = ychr.predictions[sim.indices,]

sim.interpolations = data.frame(n=adm.simulations,m=mt.simulations,y=ychr.simulations,type='sim')

interpolations = rbind(real.interpolations,sim.interpolations)

normalize = function(x) (x - min(x)) / (max(x) - min(x))
clip.extremes = function(x) {
	x[x>1] = 1
	x[x<0] = 0
	x
}

interpolations$n = clip.extremes(interpolations$n)
interpolations$m = clip.extremes(interpolations$m)
interpolations$y = clip.extremes(interpolations$y)

# interpolations = interpolations[order(interpolations$n),]
xlimits = c(25.0041673762705,25.0041673762705+0.008333334*504)
ylimits = c(-17.1958329017171,-17.1958329017171+0.008333334*384)

p = ggtern(interpolations[-25,],aes(m,n,y,col=type)) +
	geom_point() +
	xlab('mtDNA') +
	ylab('Autosomes') +
	zlab('Y') +
	theme_classic() +
	scale_color_manual(values=c('#e41a1c','#377eb8'),guide=FALSE)
ggsave(p,file='results/ancestries_interpolated.pdf',width=3,height=3,useDingbats=FALSE)

p = ggplot(interpolations[-25,],aes(m,n,col=type)) +
	geom_point() + geom_smooth(method=lm,se=FALSE) +
	xlab('mtDNA') +
	ylab('Autosomes') +
	theme_classic() +
	scale_color_manual(values=c('#e41a1c','#377eb8'),guide=FALSE) +
	coord_cartesian(xlim=c(0,1),ylim=c(0,1))
ggsave(p,file='results/ancestries_interpolated_mt.pdf',width=3,height=3,useDingbats=FALSE)


p = ggplot(interpolations[-25,],aes(y,n,col=type)) +
	geom_point() + geom_smooth(method=lm,se=FALSE) +
	xlab('Y') +
	ylab('Autosomes') +
	theme_classic() +
	scale_color_manual(values=c('#e41a1c','#377eb8'),guide=FALSE) +
	coord_cartesian(xlim=c(0,1),ylim=c(0,1))
ggsave(p,file='results/ancestries_interpolated_y.pdf',width=3,height=3,useDingbats=FALSE)

p = ggplot(interpolations[-25,],aes(y,m,col=type)) +
	geom_point() + geom_smooth(method=lm,se=FALSE) +
	xlab('Y') +
	ylab('mtDNA') +
	theme_classic() +
	scale_color_manual(values=c('#e41a1c','#377eb8'),guide=FALSE) +
	coord_cartesian(xlim=c(0,1),ylim=c(0,1))
ggsave(p,file='results/ancestries_interpolated_y_mt.pdf',width=3,height=3,useDingbats=FALSE)

# Regressions

st.interpolations = interpolations
st.interpolations = within(st.interpolations,{
	n = (n-mean(n))/sd(n)
	m = (m-mean(m))/sd(m)
	y = (y-mean(y))/sd(y)
})

m.n = lm(m~n,interpolations[interpolations$type %in% c('real','sim'),])
y.n = lm(y~n,interpolations[interpolations$type %in% c('real','sim'),])
y.m = lm(y~m,interpolations[interpolations$type %in% c('real','sim'),])
m_n = lm(m~n,interpolations[interpolations$type %in% c('real'),])
y_n = lm(y~n,interpolations[interpolations$type %in% c('real'),])
y_m = lm(y~m,interpolations[interpolations$type %in% c('real'),])

library(car)

linearHypothesis(m.n,c(0,1),1)
linearHypothesis(y.n,c(0,1),1)
linearHypothesis(y.m,c(0,1),1)
linearHypothesis(m_n,c(0,1),1)
linearHypothesis(y_n,c(0,1),1)
linearHypothesis(y_m,c(0,1),1)

# Get t- and p-value for the null hypothesis that slope = h0 (http://stats.stackexchange.com/a/111587)
# test.slope = function(model,h0=1) {
# 	if (class(model) == 'lm') model = summary(model)
# 	d.f = 1 # max(model$df)
# 	t.val = (coef(model)[2,'Estimate'] - h0) / coef(model)[2,'Std. Error']
# 	p.val = pt(t.val,d.f)
# 	c(slope=coef(model)[2,'Estimate'],intercept=coef(model)[1,'Estimate'],tval=t.val,pval=p.val,df=d.f)
# }

real.coord.matrix = t(matrix(as.numeric(adm.predictions),nrow=504))[nrow(t(matrix(as.numeric(adm.predictions),nrow=504))):1,]
real.lon = seq(25.0041673762705,25.0041673762705 + (504 - 1) * 0.008333334,0.008333334)
real.lat = seq(-17.1958329017171,-17.1958329017171 + (384 - 1) * 0.008333334,0.008333334)[384:1]
rownames(real.coord.matrix) = real.lat
colnames(real.coord.matrix) = real.lon

# Go from real matrix to fake matrix to match coordinates

fake.coord.matrix = t(real.coord.matrix[nrow(real.coord.matrix):1,])

library(reshape2)

sim.coords = melt(fake.coord.matrix)[sim.indices,1:2]
names(sim.coords) = c('x','y')
rownames(sim.coords) = NULL
all.coords = rbind(unique.coords,sim.coords)
names(all.coords) = c('lon','lat')

interpolations = data.frame(interpolations,all.coords)
# data.frame(interpolations,all.coords)

inn = interpolations[,c('n','type','lon','lat')]
inm = interpolations[,c('m','type','lon','lat')]
iny = interpolations[,c('y','type','lon','lat')]
inn$marker = 'autosomes'
inm$marker = 'mitochondrial DNA'
iny$marker = 'Y chromosome'
names(inn) = c('index','type','lon','lat','marker')
names(inm) = c('index','type','lon','lat','marker')
names(iny) = c('index','type','lon','lat','marker')
int = rbind(inn,inm,iny)

int = int[int$lon >= xlimits[1] & int$lon <= xlimits[2] & int$lat >= ylimits[1] & int$lat <= ylimits[2],]

wilcox.test(index~marker,int[int$marker %in% c('autosomes','mitochondrial DNA'),],paired=TRUE)
wilcox.test(index~marker,int[int$marker %in% c('autosomes','Y chromosome'),],paired=TRUE)
wilcox.test(index~marker,int[int$marker %in% c('mitochondrial DNA','Y chromosome'),],paired=TRUE)

wilcox.test(index~marker,int[int$type %in% 'real' & int$marker %in% c('autosomes','mitochondrial DNA'),],paired=TRUE)
wilcox.test(index~marker,int[int$type %in% 'real' & int$marker %in% c('autosomes','Y chromosome'),],paired=TRUE)
wilcox.test(index~marker,int[int$type %in% 'real' & int$marker %in% c('mitochondrial DNA','Y chromosome'),],paired=TRUE)

p = ggplot(int,aes(marker,index,color=marker)) +
	geom_violin(draw_quantiles=0.5) +
	geom_jitter(alpha=0.5) +
	theme_minimal() +
	theme(axis.title.x=element_blank(),legend.position='none') +
	ylab('Hybrid index')
ggsave(p,file='results/ancestries_interpolated_boxplot.pdf',width=2.4,height=2.4,scale=1.8,useDingbats=FALSE)

p = ggplot(int,aes(abs(lat),index,color=marker)) +
	geom_point() +
	geom_smooth(method='lm',se=FALSE) +
	theme_minimal() + coord_cartesian(ylim=c(0,1)) +
	scale_color_discrete(name='Marker set') +
	xlab(expression(paste('Latitude (',degree,'S)'))) + ylab('Hybrid index')
ggsave(p,file='results/ancestries_interpolated_latitude.pdf',width=3.6,height=2.4,scale=1.8,useDingbats=FALSE)

