#!/usr/bin/env Rscript

# Bring in Flora Jay's POPS script (http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R)
source('scripts/POPSutilities.R')

library(maptools)
library(raster)

maps = function(matrix,coord,grid,constraints=NULL,method="treshold",colorGradientsList=lColorGradients,onemap=T,onepage=T,...)

# project membership/admixture coefficients on a grid using Krig
# gradients in coefficients are represented by gradients in colors 
# if onemap=T & method="treshold" only coefficients > 0.5 are plotted
# if onemap=T & method="max" at each point the cluster for which the coefficient is maximal is plotted (even if the value is less than 0.5) 
# if onemap=F all values are plotted (since there is one cluster represented on each map, there is no overlap problem)

{

	require(fields)
	if ( (method != "treshold") & (method != "max")) {stop(paste("Unknown method",method))}
	if (class(constraints)!= "NULL") {
	   if ( nrow(grid) != nrow(constraints)*ncol(constraints) ) {
	      stop(paste("Argument grid assumes", nrow(grid), "pixels, but argument constaints assumes", nrow(constraints)*ncol(constraints),"pixels"))
	   }
	}

	if (onemap & method=="max") {
		mapsMethodMax(matrix=matrix,coord=coord,grid=grid,constraints=constraints,colorGradientsList=colorGradientsList,...)

	} else {
	K=ncol(matrix)
	if (length(colorGradientsList)<K) 
	{
		stop(paste(K,"clusters detected but only",length(colorGradientsList),"color gradient(s) defined.", 
				"You should complete colorGradientsList to have as many gradients as clusters."))
	}

	if (!onemap & onepage) {KK = as.integer(K/2); par(mfrow = c(2,KK+1)) }

	for (k in 1:K)
	{
		clust=NULL
		clust= Krig(coord, matrix[,k], theta = 10)  
		look<- predict(clust,grid) # evaluate on a grid of points
		out<- as.surface( grid, look)

		if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }

		ncolors=length(colorGradientsList[[k]])
		if (onemap) 
		{
			out[[8]][ out[[8]] < .5 ] = NA
			image(out,add=(k>1),col=colorGradientsList[[k]][(ncolors-4):ncolors],breaks=c(seq(.5,.9,.1),+200),...)
		} else {
			image(out,col=colorGradientsList[[k]][(ncolors-9):ncolors],breaks=c(-200,.1,seq(.2,.9,.1),+200),...)
		}
	}

}
}

taxon=c('Grayfoot','Kinda')

show.key=function(cluster=1,colorGradientsList=lColorGradients){
        ncolors=length(colorGradientsList[[cluster]])
        barplot(matrix(rep(1/10,10)),col=colorGradientsList[[cluster]][(ncolors-9):ncolors],main=paste0(taxon[cluster],'\nancestry'),cex.main=12/13,cex.axis=10/11)}

cliff = read.delim('data/sex_chr_genotypes.txt')
mt = cliff[!is.na(cliff$mt),]
ychr = cliff[!is.na(cliff$ychr),]
mt.coord = mt[,c('x','y')]
ychr.coord = ychr[,c('x','y')]

mt = data.frame(mt1=mt$mt,mt2=1-mt$mt)
ychr = data.frame(ychr1=ychr$ychr,ychr2=1-ychr$ychr)


pdf('results/mt_geo_plot_studyarea.pdf',width=6,height=4,useDingbats=FALSE)
layout(matrix(c(rep(1,4),2,3),2,3,byrow=FALSE), widths=c(3.8,0.2), respect=FALSE)
par(ps=22,mar=c(5,5,2,2),las=1)
asc.raster = 'data/studyarea.asc'
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster,cell_value_min=0)
maps(matrix=mt, mt.coord, grid, constraints, method="max", xlab="Longitude", ylab="Latitude",axes=FALSE,cex.lab=0.8)
# axis(side=1,at=25:29,cex.axis=0.6818182)
# axis(side=2,at=-17:-14,cex.axis=0.6818182)
# box()
points(unique(mt.coord),pch=19,col='#ffffff',cex=1.25)
par(ps=16,mar=c(5,2,5,5))
# for(k in 2:1) show.key(k)
dev.off()


pdf('results/ychr_geo_plot_studyarea.pdf',width=6,height=4,useDingbats=FALSE)
layout(matrix(c(rep(1,4),2,3),2,3,byrow=FALSE), widths=c(3.8,0.2), respect=FALSE)
par(ps=22,mar=c(5,5,2,2),las=1)
asc.raster = 'data/studyarea.asc'
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster,cell_value_min=0)
maps(matrix=ychr, ychr.coord, grid, constraints, method="max", xlab="Longitude", ylab="Latitude",axes=FALSE,cex.lab=0.8)
# axis(side=1,at=25:29,cex.axis=0.6818182)
# axis(side=2,at=-17:-14,cex.axis=0.6818182)
# box()
points(unique(ychr.coord),pch=19,col='#ffffff',cex=1.25)
par(ps=16,mar=c(5,2,5,5))
# for(k in 2:1) show.key(k)
dev.off()