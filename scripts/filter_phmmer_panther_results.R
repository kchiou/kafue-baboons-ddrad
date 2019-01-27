#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

uniprot = read.table('data/genes_to_uniprot.txt')

names(uniprot) = c('target','query','e.value','description')

# Parse target to get only macaque accessions

if (!identical(nrow(uniprot),nrow(uniprot[grep('MACMU$',uniprot$target),]))) stop('macaque?')

uniprot$target = gsub('.+?\\|([A-Z0-9]+)_MACMU$','\\1',uniprot$target)

uniprot$type = 'UniProt'
uniprot$type[nchar(uniprot$target) != 6] = 'Not UniProt'

# Split by query and return best match

uniprot.best.match = lapply(split(uniprot,uniprot$query),function(x) {
	x = x[order(x$e.value),]
	x[1,]
})

uniprot.best.match = do.call(rbind,uniprot.best.match)

write.csv(uniprot.best.match, file ='results/genes_uniprot_best_match.csv',row.names=FALSE)

rownames(uniprot.best.match) = NULL

# Filter for matches with E value < 0.0001
uniprot.matches = uniprot.best.match[uniprot.best.match$e.value < 1e-4,]


########### PREPARE PANTHER macaque data to match on UniprotKB ID

macaque.data = read.delim('data/panther_pathways_short_macaque.txt',header=FALSE)
names(macaque.data) = c('uniprot.accession','panther.family','go.molecular.function','go.biological.process','panther.protein.class','panther.pathway')

macaque.data$uniprot.accession = gsub('^MACMU.*?([A-Z0-9]{6})$','\\1',macaque.data$uniprot.accession)

macaque.pathways = macaque.data[!macaque.data$panther.pathway %in% '',]

# Find and attach "PWY" prefix to PANTHER pathways
macaque.pathways$panther.pathway = gsub('#(P[0-9]{5})>','#PWY:\\1#',macaque.pathways$panther.pathway)

# Find and attach "CMP" prefix to PANTHER components
macaque.pathways$panther.pathway = gsub('#([GP][0-9]{5});*','#CMP:\\1#',macaque.pathways$panther.pathway)

macaque.split = strsplit(macaque.pathways$panther.pathway,'#')

# Let's deal with pathways first
pathways = lapply(macaque.split,function(x) {
	x[grep('PWY:[PG][0-9]{5}',x)]
})

num.hits = sapply(pathways,length)

uniprot.list = do.call(c,lapply(1:length(num.hits),function(i) {
	rep(macaque.pathways$uniprot.accession[i],num.hits[i])
}))

pathway.list = do.call(c,pathways)

macaque.pathways = data.frame(uniprot.accession=uniprot.list,pathway.accession=pathway.list)

gene.pathway.annotations = merge(uniprot.matches,macaque.pathways,by.x='target',by.y='uniprot.accession',sort=FALSE)

# Now change the ORF accessions back into protein accessions

gene.pathway.annotations$cds.accession = gsub('_([0-9]+)_[0-9]+','.\\1',gene.pathway.annotations$query)

# Now split by protein accessions and get rid of redundant pathway classifications

gene.pathways.split = split(gene.pathway.annotations,gene.pathway.annotations$cds.accession)

gene.pathways.split = lapply(gene.pathways.split,function(x) {
	unique(x[,c('cds.accession','pathway.accession')])
})

gene.pathway.annotations = do.call(rbind,gene.pathways.split)
rownames(gene.pathway.annotations) = NULL

# Get rid of the PWY prefix
gene.pathway.annotations$pathway.accession = gsub('^PWY:','',gene.pathway.annotations$pathway.accession)

write.csv(gene.pathway.annotations,file='data/panther_pathway_annotations.csv',row.names=FALSE)

# Now figure out the components matched to proteins
components = lapply(macaque.split,function(x) {
	x[grep('CMP:[PG][0-9]{5}',x)]
})

num.hits = sapply(components,length)

uniprot.list = do.call(c,lapply(1:length(num.hits),function(i) {
	rep(macaque.pathways$uniprot.accession[i],num.hits[i])
}))

components.list = do.call(c,components)

macaque.components = data.frame(uniprot.accession=uniprot.list,component.accession=components.list)

# Get rid of prefix

macaque.components$component.accession = gsub('^CMP:','',macaque.components$component.accession)
gene.component.annotations = merge(uniprot.matches,macaque.components,by.x='target',by.y='uniprot.accession',sort=FALSE)

gene.component.annotations = gene.component.annotations[,c('query','target','component.accession')]
names(gene.component.annotations) = c('cds.accession','uniprot.accession','component.accession')

gene.component.annotations$cds.accession = gsub('_([0-9]+)_[0-9]+$','.\\1',gene.component.annotations$cds.accession)
write.csv(gene.component.annotations,file='data/panther_component_annotations.csv',row.names=FALSE)

