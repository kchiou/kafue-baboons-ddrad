#!/usr/bin/env Rscript

library(RSAGA)
library(raster)
library(maptools)

africa = read.ascii.grid('data/bio01.asc') # BIO1 variable from WorldClim dataset, clipped to Africa
africa$data[!is.na(africa$data)] = 1

x1 = africa$header$xllcorner
x2 = africa$header$xllcorner + africa$header$ncols * africa$header$cellsize
y1 = africa$header$yllcorner
y2 = africa$header$yllcorner + africa$header$nrows * africa$header$cellsize
CRS = "+proj=longlat +datum=WGS84"

africa = raster(africa$data, xmn=x1, xmx=x2, ymn=y1, ymx=y2, crs=CRS)

####################################################################################
# Crop Africa to Zambian extent
####################################################################################

zambia = readShapePoly('data/zambia_borders.shp') # AKA "ZMB_adm0.shp" from DIVA-GIS
proj4string(zambia) = CRS

ascii = crop(mask(africa,zambia),extent(zambia))

x1 = slot(slot(ascii,'extent'),'xmin')
x2 = slot(slot(ascii,'extent'),'xmax')
y1 = slot(slot(ascii,'extent'),'ymin')
y2 = slot(slot(ascii,'extent'),'ymax')
ncols = slot(ascii,'ncols')
nrows = slot(ascii,'nrows')
cellsize = (x2-x1)/ncols

# # Create a list version expected by RSAGA::write.ascii.grid
# asc = list(header=list())
# asc$header$ncols = as.integer(ncols)
# asc$header$nrows = as.integer(nrows)
# asc$header$xllcorner = x1
# asc$header$yllcorner = y1
# asc$header$cellsize = cellsize
# asc$header$nodata_value = as.integer(-9999)
# asc$header$xllcenter = x1 + cellsize / 2
# asc$header$yllcenter = y1 + cellsize / 2
# asc$data = as.matrix(ascii,nrow=nrows)

header = list(
    NCOLS=as.integer(ncols),
    NROWS=as.integer(nrows),
    XLLCENTER=x1 + cellsize / 2,
    YLLCENTER=y1 + cellsize / 2,
    CELLSIZE=round(cellsize,9),
    NODATA_VALUE=as.integer(-9999)
)

write.table(t(t(header)),file='data/zambia.asc',quote=FALSE,col.names=FALSE)
write.table(as.matrix(ascii,nrow=nrows),file='data/zambia.asc',row.names=FALSE,col.names=FALSE,append=TRUE,na='-9999')

####################################################################################
# Crop above to study extent
####################################################################################

ascii = crop(ascii,extent(c(25.5,28.7,-17.2,-14)))

x1 = slot(slot(ascii,'extent'),'xmin')
x2 = slot(slot(ascii,'extent'),'xmax')
y1 = slot(slot(ascii,'extent'),'ymin')
y2 = slot(slot(ascii,'extent'),'ymax')
ncols = slot(ascii,'ncols')
nrows = slot(ascii,'nrows')
cellsize = (x2-x1)/ncols

header = list(
    NCOLS=as.integer(ncols),
    NROWS=as.integer(nrows),
    XLLCENTER=x1 + cellsize / 2,
    YLLCENTER=y1 + cellsize / 2,
    CELLSIZE=round(cellsize,9),
    NODATA_VALUE=as.integer(-9999)
)

write.table(t(t(header)),file='data/studyarea.asc',quote=FALSE,col.names=FALSE)
write.table(as.matrix(ascii,nrow=nrows),file='data/studyarea.asc',row.names=FALSE,col.names=FALSE,append=TRUE,na='-9999')

####################################################################################
# Crop Africa to Kafue extent
####################################################################################

kafue = readShapePoly('data/kafue_parkboundary.shp') # courtesy of ZAWA
proj4string(zambia) = CRS

ascii = crop(mask(africa,kafue),extent(kafue))

x1 = slot(slot(ascii,'extent'),'xmin')
x2 = slot(slot(ascii,'extent'),'xmax')
y1 = slot(slot(ascii,'extent'),'ymin')
y2 = slot(slot(ascii,'extent'),'ymax')
ncols = slot(ascii,'ncols')
nrows = slot(ascii,'nrows')
cellsize = (x2-x1)/ncols

header = list(
    NCOLS=as.integer(ncols),
    NROWS=as.integer(nrows),
    XLLCENTER=x1 + cellsize / 2,
    YLLCENTER=y1 + cellsize / 2,
    CELLSIZE=round(cellsize,9),
    NODATA_VALUE=as.integer(-9999)
)

write.table(t(t(header)),file='data/kafue.asc',quote=FALSE,col.names=FALSE)
write.table(as.matrix(ascii,nrow=nrows),file='data/kafue.asc',row.names=FALSE,col.names=FALSE,append=TRUE,na='-9999')