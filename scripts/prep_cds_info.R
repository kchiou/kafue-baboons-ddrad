#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

arguments = commandArgs(trailingOnly=TRUE)

cds.file = arguments[1]
gtf.file = arguments[2]

headers = system(paste('cat',cds.file,'| grep "^>"'),intern=TRUE)

gtf = system(paste('cat ',gtf.file,' | grep -v "^#"'),intern=TRUE)
gtf.split = strsplit(gtf,'\\t')
gtf = as.data.frame.matrix(do.call(rbind,gtf.split))
names(gtf) = c('seqname','source','feature','start','end','score','strand','frame','attribute')
gtf$gene_id = NA
gtf$gene_id[grep('gene_id',gtf$attribute)] = gsub('.*?gene_id "(.+?)".+','\\1',gtf$attribute[grep('gene_id',gtf$attribute)])
gtf$gene_name = NA
gtf$gene_name[grep('gene_name',gtf$attribute)] = gsub('.*?gene_name "(.+?)".+','\\1',gtf$attribute[grep('gene_name',gtf$attribute)])

gtf.genes = gtf[gtf$feature %in% 'gene',]
gtf.genes = gtf.genes[,c('seqname','source','gene_id','gene_name','start','end','score','strand','frame')]

# Pull out CDS info with columns: ID, chromosome/scaffold, chromosome number, start, stop
ids = gsub('>(.+?) cds (.+?):.+?:(.+?):(.+?):(.+?):(.+?) (.+?) .+','\\1 \\2 \\3 \\4 \\5 \\6 \\7',headers)

ids = paste(1:length(ids),ids)

genes = strsplit(ids,' ')
genes = do.call(rbind,genes)
genes = as.data.frame.matrix(genes)

names(genes) = c('id','nam','src','chr','sta','sto','str','gene')
genes$gene = gsub('gene:(.+?)\\.[0-9]+','\\1',genes$gene)

genes = within(genes,{
	id = as.numeric(id)
	src = factor(src)
	chr[src == 'scaffold'] = NA
	chr = factor(chr,levels=c(1:20,'X','MT'))
	sta = as.numeric(sta)
	sto = as.numeric(sto)
	str = factor(str)
	gene = as.character(gene)
})

genes$id2 = gsub('\\.','_',genes$id)

# Start of window is 50 bp before the start of the CDS (negative is okay)
genes$sta.50 = genes$sta - 50

# End of window is 50 bp after the end of the CDS (past end of chromosome is okay)
genes$sto.50 = genes$sto + 50

# Midpoint of window is the average of the start and stop position (rounded to whole number)
genes$mid = floor((genes$sta + genes$sto) / 2)

gtf.genes = within(gtf.genes,{
	nam = as.character(gene_id)
	src = factor(source)
	chr = seqname
	gene_name = as.character(gene_name)
	sta = as.numeric(start)
	sto = as.numeric(end)
	score = factor(score)
	strand = factor(strand)
	frame = factor(frame)
})

gtf.genes = gtf.genes[,c('nam','gene_name','chr','sta','sto','score','strand','frame','src')]

gtf.genes$sta.50 = gtf.genes$sta - 50
gtf.genes$sto.50 = gtf.genes$sto - 50
gtf.genes$mid = floor((gtf.genes$sta + gtf.genes$sto) / 2)


write.csv(genes,file='results/baboon_cds_info.csv',row.names=FALSE)
write.csv(gtf.genes,file='results/baboon_genes_info.csv',row.names=FALSE)