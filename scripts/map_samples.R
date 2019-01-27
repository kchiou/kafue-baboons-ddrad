#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

library(maptools)
library(mapplots)
library(methods)

zambia = readShapePoly('data/zambia_borders.shp')
proj4string(zambia) = '+proj=longlat +datum=WGS84'

map.kafue = function(pts=FALSE,xlim=NULL,ylim=NULL,add=FALSE,extras=TRUE) {
	require(maptools)

	if (!is.null(lake)) {
		xlim = c(25.3,26.4)
		ylim = c(-16.3,-15.2)
	}

	riverpoly = readShapePoly('data/kafue_riverpoly.shp')
	parkboundary = readShapePoly('data/kafue_parkboundary.shp')

	proj4string(riverpoly) = "+proj=longlat +datum=WGS84"
	proj4string(parkboundary) = "+proj=longlat +datum=WGS84"

	plot(parkboundary,col='forestgreen',xlim=xlim,ylim=ylim,add=add)
	if (extras) {
		plot(riverpoly,col='blue',add=TRUE)
	}
}

# Get samples that passed QC
mds = read.table('results/mds.txt',header=TRUE)

ind_info = read.csv('data/individual_info.csv')
ind_info = ind_info[ind_info$Individual.ID %in% mds$FID,]
ind_info$Sample.type[!ind_info$Sample.type %in% 'feces'] = 'blood'
locales = unique(ind_info[,c('Group','Taxon','Longitude','Latitude','Sample.type','Locality')])
locales$n = table(mds$group)[locales$Group]
locales = locales[order(locales$Group),]
locales$Locality = factor(locales$Locality,levels=locales$Locality)
locales$Group = factor(locales$Group,levels=locales$Group)
rownames(locales) = NULL

# Calculate diameter for circles (pch = 21)
locales$size = sqrt(locales$n/pi) * 2
locales=locales[order(locales$size,decreasing=TRUE),]

# Calculate diagonal for square diamonds (pch = 23)
# locales$size = sqrt(locales$n) * sqrt(2)

taxa = c('Kinda' = '#fdb863', 'Hybrid' = '#f7f7f7', 'Grayfoot' = '#b2abd2')

# Get range of sample size values
n.range = min(locales$n):max(locales$n)
n.intervals = c(n.range[1],n.range[round(length(n.range) / 16)],n.range[round(length(n.range) / 8)],n.range[round(length(n.range))])

pdf(file='results/samples_map_kafue.pdf',useDingbats=FALSE)
	plot(zambia,xlim=c(25.5, 28.7),ylim=c(-17.2,-14),lwd=2,main='Samples')
	map.kafue(add=TRUE)
	points(locales[,c('Longitude','Latitude')],col='black',pch=21,bg=taxa[locales$Taxon],cex=locales$size*2/3)
	text(locales[,c('Longitude','Latitude')],labels=locales$Group,cex=1/2)
	legend("bottomleft",inset=0.025, names(taxa),fill=taxa,horiz=FALSE,cex=0.8)
	legend("topright",inset=0.025, paste(': ',levels(locales$Locality)),pch=levels(locales$Group),horiz=FALSE,cex=0.8)
	legend("bottomright",inset=0.025, c(paste0('= ',min(locales$n),' samples'),'= 12 samples'),pch=21,pt.cex=sqrt(c(min(locales$n),12)/pi) * 2 * 2/3,horiz=FALSE,cex=0.8,bg='#ffffff',y.intersp=1.5)
dev.off()

pdf(file='results/samples_map_zambia.pdf',useDingbats=FALSE)
	plot(zambia,main='Samples')
	map.kafue(add=TRUE,extras=FALSE)
	points(locales[,c('Longitude','Latitude')],col='black',pch=21,bg=taxa[locales$Taxon],cex=1)
	legend("topleft",inset=0.025, names(taxa),fill=taxa,horiz=FALSE,cex=0.8)
dev.off()

adm.2 = read.table('results/all.cleaned.LDpruned.2.Q')
# Swap in consistent order, based on total amount of ancestry ascribed to component
col.order = order(colSums(adm.2))
adm.2 = data.frame(adm.2[,col.order])
names(adm.2) = paste0("V", 1:2)

if (adm.2[which(mds$ind_id == "BZ11-001"),]$V1 > 0.5) {
	adm.2 = data.frame(adm.grayness=adm.2$V2, adm.kindness=adm.2$V1)
} else {
	adm.2 = data.frame(adm.grayness=adm.2$V1, adm.kindness=adm.2$V2)
}

sample.info = data.frame(mds,adm.2)
locales$adm.grayness = as.numeric(tapply(sample.info$adm.grayness,sample.info$group,mean)[locales$Group])
locales$adm.kindness = 1 - locales$adm.grayness

pdf(file='results/samples_map_adm_ancestry.pdf',useDingbats=FALSE,height=6,width=6)
	plot(zambia,xlim=c(25.5, 28.7),ylim=c(-17.2,-14),lwd=2)
	map.kafue(add=TRUE)
	for (i in 1:nrow(locales)) {
		r = locales[i,]
		add.pie(z=c(r$adm.grayness,r$adm.kindness),x=r$Longitude,y=r$Latitude,radius=r$size/(60),col=c("#b2abd2","#fdb863"),labels=NA,clockwise=FALSE)
	}
	legend("bottomleft",inset=0.025, names(taxa)[-2],fill=taxa[-2],horiz=FALSE,cex=0.8)
	legend("bottomright",inset=0.025, c(paste0('= ',min(locales$n),' samples'),'= 12 samples'),pch=21,pt.cex=sqrt(c(min(locales$n),12)/pi) * 2 * 2/3,horiz=FALSE,cex=0.8,bg='#ffffff',y.intersp=1.5)
dev.off()
