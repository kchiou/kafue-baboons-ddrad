#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

arguments = commandArgs(trailingOnly=TRUE)

if (length(arguments) < 2) {
	stop("Two arguments required.\n", call.=FALSE)
}

# Set input and output filenames
FULL = arguments[1]
PATH = arguments[2]
OUT = paste0(gsub('\\.[A-z0-9]*$','',PATH),'.bgc')

# ========================================================================================
# === Pull out genotype info and metadata from input file
# ========================================================================================

# Use sed and awk to identify start and end rows / columns for the genotype matrix
LENGTH = system(paste('wc -l',PATH,'| sed "s/^ *//g" | tr -s " " | cut -d " " -f 1'),intern=TRUE)
START = system(paste('grep -nr "#CHROM"',PATH,' | sed -E "s/.*:([0-9]+):#CHROM/\\1:#CHROM/g" | cut -d ":" -f 1'),intern=TRUE)
END = system(paste('wc -l',PATH,'| sed "s/^ *//g" | tr -s " " | cut -d " " -f 1'),intern=TRUE)
SNP_COL = system(paste('cat',PATH,'|',paste0('sed -n ',START,'p'),'| awk -v RS="\\t" "/ID/ { print NR }"'),intern=TRUE)
FIRST_COL = as.numeric(system(paste('cat',PATH,'| awk "NR ==',START,'{ print }" | awk -v RS="\\t" "/FORMAT/ { print NR }"'),intern=TRUE)) + 1
LAST_COL = system(paste('cat',PATH,'| awk "NR ==',START,'{ print }" | awk "{ print NF }"'),intern=TRUE)

# Pull out individual IDs and parse them
ind_id = system(paste0('cat ',PATH,' | sed -n ',START,'p | cut -f ',FIRST_COL,'-',LAST_COL),intern=TRUE)
ind_id = strsplit(ind_id,'\\t')[[1]]
ind_id = gsub('\\.PE.*','',ind_id)             # annoying formatting issue

snp_id = system(paste0('cat ',PATH,' | sed -n ',as.numeric(START)+1,',',END,'p | cut -f ',SNP_COL),intern=TRUE)

# Pull a table that has group info (post filtering)
mds = read.table('results/mds.txt',header=TRUE)
groups = factor(mds$group)
names(groups) = mds$FID

# Pull out SNP matrix
snp.matrix = system(paste0('cat ',PATH,' | sed -n ',as.numeric(START)+1,',',END,'p | cut -f ',FIRST_COL,'-',LAST_COL),intern=TRUE)

# Substitute missing data with string "NA"
snp.matrix = gsub('\\./\\.','NA',snp.matrix)

# For each row, explode the string by the tab character
snp.matrix = strsplit(snp.matrix,'\\t')

# bind rows together
snp.matrix = do.call(rbind,snp.matrix)
snp.matrix[snp.matrix %in% 'NA'] = NA

# ========================================================================================
# === Pull out genotype info and metadata from full file containing read counts
# ========================================================================================

LENGTH.full = system(paste('wc -l',FULL,'| sed "s/^ *//g" | tr -s " " | cut -d " " -f 1'),intern=TRUE)
START.full = system(paste('grep -nr "#CHROM"',FULL,' | sed -E "s/.*:([0-9]+):#CHROM/\\1:#CHROM/g" | cut -d ":" -f 1'),intern=TRUE)
END.full = system(paste('wc -l',FULL,'| sed "s/^ *//g" | tr -s " " | cut -d " " -f 1'),intern=TRUE)
CHR_COL.full = system(paste('cat',FULL,'|',paste0('sed -n ',START.full,'p'),'| awk -v RS="\\t" "/#CHROM/ { print NR }"'),intern=TRUE)
POS_COL.full = system(paste('cat',FULL,'|',paste0('sed -n ',START.full,'p'),'| awk -v RS="\\t" "/POS/ { print NR }"'),intern=TRUE)
FIRST_COL.full = as.numeric(system(paste('cat',FULL,'|',paste0('sed -n ',START.full,'p'),'| awk -v RS="\\t" "/FORMAT/ { print NR }"'),intern=TRUE)) + 1
LAST_COL.full = system(paste('cat',FULL,'|',paste0('sed -n ',START.full,'p'),'| awk "{ print NF }"'),intern=TRUE)
FORMAT_COL.full = FIRST_COL.full - 1

# Pull out individual IDs and parse them
ind_id.full = system(paste0('cat ',FULL,' | sed -n ',START.full,'p | cut -f ',FIRST_COL.full,'-',LAST_COL.full),intern=TRUE)
ind_id.full = strsplit(ind_id.full,'\\t')[[1]]
ind_id.full = gsub('\\.PE.*','',ind_id.full)             # annoying formatting issue

# Concatenate together chromosome IDs and position IDs
snp_id.full = paste(
	system(paste0('cat ',FULL,' | sed -n ',as.numeric(START.full)+1,',',END.full,'p | cut -f ',CHR_COL.full),intern=TRUE),
	system(paste0('cat ',FULL,' | sed -n ',as.numeric(START.full)+1,',',END.full,'p | cut -f ',POS_COL.full),intern=TRUE),
	sep=':'
)

# Pull out genotypic metadata (VCF "FORMAT" column)
gt.format = unique(system(paste0('cat ',FULL,' | sed -n ',as.numeric(START.full)+1,',',END.full,'p | cut -f ',FORMAT_COL.full),intern=TRUE))

# Flip an error if there is more than one genotype format
if (length(gt.format) > 1) stop('Only one call format supported')

# Find where the AD (allelic depth) call is in the format string
AD.position = match('AD',strsplit(gt.format,':')[[1]])

# Flip an error if there are any SNPs in the input file that are not in the full file
if (!all(snp_id %in% snp_id.full)) stop('Missing SNPs\n')

snp.line.numbers = match(snp_id,snp_id.full) + as.numeric(START.full)
ind.col.numbers = match(ind_id,ind_id.full) + as.numeric(FIRST_COL.full) - 1

# Power through and extract all the matching lines and columns

# First cycle through all SNPs and search for runs of consecutives (sed will fail if given too many inputs)
snp.distances = c(1,c(snp.line.numbers[2:length(snp.line.numbers)]) - snp.line.numbers[1:(length(snp.line.numbers) - 1)])
i = 1
current.group = 0
snp.groups = numeric(0)
while (i <= length(snp.distances)) {
	if (snp.distances[i] == 1) {
	} else {
		current.group = current.group + 1
	}
	snp.groups = c(snp.groups,current.group)
	i = i + 1
}

snp.line.groups = as.character(do.call(c,lapply(split(snp.line.numbers,snp.groups),function(x) {
	if (length(x) == 1) x else paste0(x[1],',',x[length(x)])
})))

snp.matrix.full = system(paste0('cat ',FULL,' | sed -n -e ',paste(snp.line.groups,collapse='p -e '),'p | cut -f ',paste(ind.col.numbers,collapse=',')),intern=TRUE)

snp.matrix.full = gsub('\\./\\.','NA',snp.matrix.full)

# For each row, explode the string by the tab character
snp.matrix.full = strsplit(snp.matrix.full,'\\t')

# Some genotypes weirdly looked like this: "./.:.:1" --> this was unforeseen but fixable
snp.matrix.full = lapply(snp.matrix.full,function(x) gsub('^NA.*','NA',x))

rc.matrix = lapply(snp.matrix.full,function(x) {
	gsub(paste0('^',paste(c(rep('.*?',AD.position-1),'(.*?)'),collapse=':'), if (AD.position == length(strsplit(gt.format,':')[[1]])) '$' else ':.*'),'\\1',x)
})

rc.matrix = do.call(rbind,rc.matrix)
rc.matrix[rc.matrix %in% 'NA'] = NA

if (!identical(is.na(snp.matrix),is.na(rc.matrix))) {
	warning('Genotype and read count have different missing values!\nInserting corresponding missing values.')
	rc.matrix[is.na(snp.matrix)] = NA
}

# Read counts for missing data should be 0
rc.matrix[is.na(rc.matrix)] = '0,0'

# Convert the comma delimiter to a space
rc.matrix = gsub(',',' ',rc.matrix)

# Get read count stats
allele1 = gsub('([0-9]+) ([0-9]+)','\\1',rc.matrix)
allele2 = gsub('([0-9]+) ([0-9]+)','\\2',rc.matrix)
mode(allele1) = 'numeric'
mode(allele2) = 'numeric'
snp.allele.counts = allele1 + allele2

snp_info = data.frame(id=0:(length(snp_id)-1),label=paste('locus',0:(length(snp_id)-1)),name=snp_id,chr=gsub(':.*$','',snp_id),pos=gsub('^.*?:','',snp_id))
write.table(snp_info,file='results/bgc_snp_info.txt',quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')

# ========================================================================================
# === Parse individual/group info
# ========================================================================================

ind_info = read.csv('data/individual_info.csv')

# Keep only individuals that passed filters
ind_info = ind_info[ind_info$Individual.ID %in% ind_id,]

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

write.table(ind_info,file='results/bgc_ind_info.txt',quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')

grp_info = unique(ind_info[,c('pop','group','locality')])
row.names(grp_info) = NULL

write.table(grp_info,file='results/bgc_grp_info.txt',quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')

# ========================================================================================
# === Write read counts
# ========================================================================================

# Parental populations, genotype uncertainty: The data for each locus begins with
#      a line that gives the locus number (e.g., locus 32). The locus line is followed by
#      a data line for each individual in the population. Each data line gives the number
#      of sequences (i.e., reads) of each allele for each individual. If no reads are
#      observed for an individual, the data line should simply contain zeros. A separate
#      file is required for each of the two parental populations. A short example with
#      two loci and two alleles is given in Table 3.


parse.bgc = function(mat,snps,pops=NULL) {
	if (is.null(pops)) {
		paste(paste0('locus ',snps,'\n',apply(mat,1,paste,collapse='\n')),collapse='\n')
	} else {

		# Split indices per population, subset the genotype matrix, and bind genotypes together
		indices = split(1:length(pops),pops)

		# Subset genotype matrix for each population
		mat0 = do.call(cbind,lapply(names(indices),function(x) {
			apply(mat[,indices[[x]]],1,function(y) paste0(x,'\n',paste(y,collapse='\n')))
		}))

		paste(paste0('locus ',snps,'\n',apply(mat0,1,paste,collapse='\n')),collapse='\n')

	}
}

colnames(rc.matrix) = ind_id
rownames(rc.matrix) = snp_info$id

pure1 = parse.bgc(
			mat = rc.matrix[,ind_info$id[ind_info$pop %in% 'pure 1']],
			snps = snp_info$id)

pure2 = parse.bgc(
			mat = rc.matrix[,ind_info$id[ind_info$pop %in% 'pure 2']],
			snps = snp_info$id)

hybrids = parse.bgc(
			mat = rc.matrix[,ind_info$id[ind_info$pop %in% 'pop 0']],
			snps = snp_info$id,
			pops = rep('pop 0',length(grep('^pop',ind_info$pop))))

write(pure1,file='data/bgc_p0in.txt')
write(pure2,file='data/bgc_p1in.txt')
write(hybrids,file='data/bgc_admixedin.txt')

