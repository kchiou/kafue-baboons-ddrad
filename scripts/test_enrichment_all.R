#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

selected.genes = read.csv('results/full_scans.csv')

# ------------------------------------------------------------------------------
# --- Prep baboon and GO annotations
# ------------------------------------------------------------------------------

library(biomaRt)
library(topGO)

panu = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='panubis_gene_ensembl')

panu.go = getBM(attributes=c('ensembl_gene_id', 'go_id','name_1006','definition_1006'), filters = 'ensembl_gene_id', values = selected.genes$nam, mart = panu)

gene2go = lapply(unique(panu.go$ensembl_gene_id),function(x){sort(panu.go[panu.go$ensembl_gene_id==x,'go_id'])})

names(gene2go) = unique(panu.go$ensembl_gene_id)

go.def = unique(subset(panu.go,select=c('go_id','name_1006','definition_1006')))

# ------------------------------------------------------------------------------
# --- Prep baboon and PANTHER annotations
# ------------------------------------------------------------------------------

# Get gene info for later
gene.info = read.csv('results/baboon_genes_info.csv')
gene.match = gene.info[gene.info$nam %in% selected.genes$nam,]
gene.match = gene.match[,c('nam','gene_name','chr','sta','sto')]

cds.info = read.csv('results/baboon_cds_info.csv')

panther.pathways = read.csv('data/panther_pathway_annotations.csv')

genes.with.pthr = merge(selected.genes,gene.match,by='nam')

names(genes.with.pthr)[names(genes.with.pthr) %in% 'nam'] = 'gene.accession'
names(genes.with.pthr)[names(genes.with.pthr) %in% 'sta'] = 'gene.start'
names(genes.with.pthr)[names(genes.with.pthr) %in% 'sto'] = 'gene.stop'

genes.with.pthr = genes.with.pthr[order(genes.with.pthr$chr,genes.with.pthr$gene.start),]

# Create subsets that are significant with alpha=0.01 and alpha=0.05
# genes.with.go$sig.fst.01 = genes.with.go$pval.fst <= 0.01
# genes.with.go$sig.fst = genes.with.go$pval.fst <= 0.05
genes.with.pthr$sig.fst.fdr = genes.with.pthr$pval.fst.fdr <= 0.05

# Merge in the CDS info to match PTHR accessions
cds.minimal = cds.info
cds.minimal = cds.minimal[,c('gene','nam')]
names(cds.minimal) = c('gene.accession','cds.accession')

genes.with.pthr = merge(genes.with.pthr,cds.minimal,by='gene.accession',all.x=TRUE,sort=FALSE)

# Merge datasets together
genes.with.pthr = merge(genes.with.pthr,panther.pathways,by='cds.accession',all.x=FALSE,all.y=FALSE,sort=FALSE)
rownames(genes.with.pthr) = NULL
names(genes.with.pthr)[names(genes.with.pthr) %in% 'pathway.accession'] = 'pthr.pathway'

# Now that pathways are in, we no longer need the CDS accession
genes.with.pthr$cds.accession = NULL

# Get rid of duplicates
genes.with.pthr = unique(genes.with.pthr)

# Naming convention to match the previous GO analyses
pthr.with.genes = genes.with.pthr

# Get rid of NAs, as NAs will screw up category counts later
pthr.with.genes.fst = pthr.with.genes[!is.na(pthr.with.genes$fst.fst),]
pthr.with.genes.bgc = pthr.with.genes[!is.na(pthr.with.genes$bgc.a),]

pthr.definitions = read.delim('data/panther_pathway_definitions.txt',header=FALSE)
names(pthr.definitions) = c('pthr.accession','pthr.definition')

# ------------------------------------------------------------------------------
# --- Fst GO enrichment
# ------------------------------------------------------------------------------

fst.genes = subset(selected.genes,!is.na(fst.fst))

pval.fst = fst.genes$pval.fst
names(pval.fst) = fst.genes$nam

fst.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=pval.fst,geneSel=function(x) x < 0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)

fst.topgo.fet.pc = runTest(fst.topgo,algorithm='parentchild',statistic='fisher')
fst.topgo.ks.w1 = runTest(fst.topgo,algorithm='weight01',statistic='ks')

fst.topgo.results = data.frame(
	go.id = names(fst.topgo.ks.w1@score),
	p.ks = fst.topgo.ks.w1@score,
	p.fet = fst.topgo.fet.pc@score[names(fst.topgo.ks.w1@score)],
	fdr.ks = p.adjust(fst.topgo.ks.w1@score,'fdr'),
	fdr.fet = p.adjust(fst.topgo.fet.pc@score[names(fst.topgo.ks.w1@score)],'fdr')
)

# fst.topgo.results = merge(fst.topgo.results,go.def,by.x='go.id',by.y='go_id',all.x=TRUE,all.y=FALSE)

fst.topgo.results = fst.topgo.results[order(fst.topgo.results$p.ks),]

library(XML)
library(tidyr)

go.join = go.def$name_1006
names(go.join) = go.def$go_id

fst.topgo.results$go.name = go.join[fst.topgo.results$go.id]

# Fetch missing go names from the AMIGO site
to.do = which(is.na(fst.topgo.results$go.name))

for (i in to.do) {
	cat(i,'\n')
	search.url = paste0('http://amigo.geneontology.org/amigo/term/',fst.topgo.results$go.id[i])
	fst.topgo.results$go.name[i] = search.url %>% htmlParse %>% xmlChildren %>% `[[`(3) %>% xmlChildren %>% `[[`(4) %>% xmlChildren %>% `[[`(24) %>% xmlChildren %>% `[[`(6) %>% xmlChildren %>% `[[`(2) %>% xmlChildren %>% `[[`('text') %>% xmlValue
}

# ------------------------------------------------------------------------------
# --- BGC GO enrichment
# ------------------------------------------------------------------------------

bgc.genes = subset(selected.genes,!is.na(bgc.b))

a.bgc = bgc.genes$bgc.a
b.bgc = bgc.genes$bgc.b
p.a.pos = 1 - bgc.genes$a.pos
p.a.neg = 1 - bgc.genes$a.neg
p.b.pos = 1 - bgc.genes$b.pos
p.b.neg = 1 - bgc.genes$b.neg

names(a.bgc) = names(b.bgc) = names(p.a.pos) = names(p.a.neg) = names(p.b.pos) = names(p.b.neg) = bgc.genes$nam

# enrichment for negative alpha
#bgc.a1.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=a.bgc,geneSel=function(x) x,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)
bgc.a1.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=p.a.neg,geneSel=function(x) x,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)
bgc.a1.topgo.ks.w1 = runTest(bgc.a1.topgo,algorithm='weight01',statistic='ks')

# enrichment for positive alpha
#bgc.a2.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=-a.bgc,geneSel=function(x) x,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)
bgc.a2.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=p.a.pos,geneSel=function(x) x,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)
bgc.a2.topgo.ks.w1 = runTest(bgc.a2.topgo,algorithm='weight01',statistic='ks')

# enrichment for negative beta
#bgc.b1.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=b.bgc,geneSel=function(x) x,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)
bgc.b1.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=p.b.neg,geneSel=function(x) x,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)
bgc.b1.topgo.ks.w1 = runTest(bgc.b1.topgo,algorithm='weight01',statistic='ks')

# enrichment for positive beta
#bgc.b2.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=-b.bgc,geneSel=function(x) x,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)
bgc.b2.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=p.b.pos,geneSel=function(x) x,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=gene2go)
bgc.b2.topgo.ks.w1 = runTest(bgc.b2.topgo,algorithm='weight01',statistic='ks')

bgc.topgo.results = data.frame(
	go.id = names(bgc.b2.topgo.ks.w1@score),
	p.a1 = bgc.a1.topgo.ks.w1@score[names(bgc.b2.topgo.ks.w1@score)],
	p.a2 = bgc.a2.topgo.ks.w1@score[names(bgc.b2.topgo.ks.w1@score)],
	p.b1 = bgc.b1.topgo.ks.w1@score[names(bgc.b2.topgo.ks.w1@score)],
	p.b2 = bgc.b2.topgo.ks.w1@score[names(bgc.b2.topgo.ks.w1@score)]
)

go.join = go.def$name_1006
names(go.join) = go.def$go_id

bgc.topgo.results$go.name = go.join[bgc.topgo.results$go.id]

# Fetch missing go names from the AMIGO site
to.do = which(is.na(bgc.topgo.results$go.name))

for (i in to.do) {
	cat(i,'\n')
	search.url = paste0('http://amigo.geneontology.org/amigo/term/',bgc.topgo.results$go.id[i])
	bgc.topgo.results$go.name[i] = search.url %>% htmlParse %>% xmlChildren %>% `[[`(3) %>% xmlChildren %>% `[[`(4) %>% xmlChildren %>% `[[`(24) %>% xmlChildren %>% `[[`(6) %>% xmlChildren %>% `[[`(2) %>% xmlChildren %>% `[[`('text') %>% xmlValue
}

write.table(bgc.topgo.results,file='results/kafue_bgc_topgo.tsv',sep='\t',row.names=FALSE,quote=FALSE)

# ------------------------------------------------------------------------------
# --- Generic function for conducting enrichment analysis (not hierarchical like topGO)
# ------------------------------------------------------------------------------

# Enrichment test. Compare stats associated with the functional term to stats not associated with it
enrich.test = function(annotation.dataset,gene.dataset,stat.column,alt='greater',accession.column='go.accession',annotation.column='go.annotations') {
	do.call(c,lapply(annotation.dataset[[accession.column]],function(x) {
		enrich.dataframe = data.frame(
			in.group = factor(gene.dataset[[annotation.column]] %in% x,levels=c('TRUE','FALSE')),
			stat = gene.dataset[[stat.column]]
		)
		
		wilcox.test(stat ~ in.group,data=enrich.dataframe,alternative=alt)$p.value
	}))
}


# ------------------------------------------------------------------------------
# --- Fst PANTHER enrichment
# ------------------------------------------------------------------------------

pthr.categories.fst = data.frame(pthr.accession=names(table(pthr.with.genes.fst$pthr.pathway)))

pthr.categories.fst$enrichment.fst = enrich.test(pthr.categories.fst,pthr.with.genes.fst,stat.column='fst.fst',alt='greater',accession.column='pthr.accession',annotation.column='pthr.pathway')
pthr.categories.fst$enrichment.p = enrich.test(pthr.categories.fst,pthr.with.genes.fst,stat.column='pval.fst',alt='less',accession.column='pthr.accession',annotation.column='pthr.pathway')

pthr.categories.fst = merge(pthr.categories.fst,pthr.definitions,by='pthr.accession',all.x=TRUE,all.y=FALSE,sort=FALSE)

write.table(pthr.categories.fst,file='results/kafue_fst_panther.tsv',sep='\t',row.names=FALSE,quote=FALSE)


# ------------------------------------------------------------------------------
# --- BGC PANTHER enrichment
# ------------------------------------------------------------------------------

pthr.categories.bgc = data.frame(pthr.accession=names(table(pthr.with.genes.bgc$pthr.pathway)))

pthr.categories.bgc$enrichment.a1 = enrich.test(pthr.categories.bgc,pthr.with.genes.bgc,stat.column='bgc.a',alt='greater',accession.column='pthr.accession',annotation.column='pthr.pathway')
pthr.categories.bgc$enrichment.a2 = enrich.test(pthr.categories.bgc,pthr.with.genes.bgc,stat.column='bgc.a',alt='less',accession.column='pthr.accession',annotation.column='pthr.pathway')
pthr.categories.bgc$enrichment.b1 = enrich.test(pthr.categories.bgc,pthr.with.genes.bgc,stat.column='bgc.b',alt='greater',accession.column='pthr.accession',annotation.column='pthr.pathway')
pthr.categories.bgc$enrichment.b2 = enrich.test(pthr.categories.bgc,pthr.with.genes.bgc,stat.column='bgc.b',alt='less',accession.column='pthr.accession',annotation.column='pthr.pathway')

pthr.categories.bgc$enrichment.p.a1 = enrich.test(pthr.categories.bgc,pthr.with.genes.bgc,stat.column='a.neg',alt='less',accession.column='pthr.accession',annotation.column='pthr.pathway')
pthr.categories.bgc$enrichment.p.a2 = enrich.test(pthr.categories.bgc,pthr.with.genes.bgc,stat.column='a.pos',alt='less',accession.column='pthr.accession',annotation.column='pthr.pathway')
pthr.categories.bgc$enrichment.p.b1 = enrich.test(pthr.categories.bgc,pthr.with.genes.bgc,stat.column='b.neg',alt='less',accession.column='pthr.accession',annotation.column='pthr.pathway')
pthr.categories.bgc$enrichment.p.b2 = enrich.test(pthr.categories.bgc,pthr.with.genes.bgc,stat.column='b.pos',alt='less',accession.column='pthr.accession',annotation.column='pthr.pathway')

pthr.categories.bgc = merge(pthr.categories.bgc,pthr.definitions,by='pthr.accession',all.x=TRUE,all.y=FALSE,sort=FALSE)

write.table(pthr.categories.bgc,file='results/kafue_bgc_panther.tsv',sep='\t',row.names=FALSE,quote=FALSE)








