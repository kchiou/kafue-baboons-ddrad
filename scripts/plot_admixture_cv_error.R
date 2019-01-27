#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
prefix = args[1]
if (is.na(prefix)) {
	prefix=""
}

admix_cv_error = read.csv(paste("reports/", prefix, "ADMIXTURE_log_cv.out", sep=""), sep=":", header=FALSE)

lowest.k = which (admix_cv_error$V2 == min(admix_cv_error$V2))


p = ggplot(admix_cv_error,aes(V3,V2)) +
	geom_line() +
	geom_point(size=6,color='white') +
	geom_point(size=2) +
	geom_point(aes(V3,V2),data=admix_cv_error[2,],pch=21,size=4,color='red') +
	geom_text(aes(V3,V2,label=round(V2,3)),data=admix_cv_error[2,],
		size=3,color='red',nudge_x=-2/3) +
	scale_x_continuous(breaks=seq(2,10,2)) +
	scale_y_continuous(breaks=seq(0.36,0.46,0.02)) +
	theme_classic() +
	xlab(expression(italic(K))) +
	ylab('Cross-validation error')
ggsave(p,file='reports/admixture_cv.pdf',width=4.5,height=4.5,useDingbats=FALSE)

pdf(file=paste('reports/', prefix, 'ADMIXTURE_CV_plot.pdf', sep=""),width=6,height=6,useDingbats=FALSE)

	plot(admix_cv_error$V2 ~ seq(nrow(admix_cv_error)), type="b", pch=19,
		xlab="K", ylab="Cross-Validation Error")
	
	ideal.k = 2
	points(admix_cv_error$V2[ideal.k] ~ ideal.k, col="red", cex=2)
	text(ideal.k, admix_cv_error$V2[ideal.k], 
		labels=round(admix_cv_error$V2[ideal.k], digits=3), 
		cex=0.8, pos=2, col="red")

	lowest.k = which (admix_cv_error$V2 == min(admix_cv_error$V2))
	points(admix_cv_error$V2[lowest.k] ~ lowest.k, col="red", cex=2)
	text(lowest.k, admix_cv_error$V2[lowest.k], 
		labels=round(admix_cv_error$V2[lowest.k], digits=3), 
		cex=0.8, pos=3, col="red")
		
dev.off()