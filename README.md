# kafue-baboons-ddrad

This repository contains all analysis code for two related projects on the baboons of Kafue National Park. Baboon DNA was prepared using [double-digest RADseq](https://doi.org/10.1371/journal.pone.0037135) and used for two related projects. The citations for the papers are below (citations will be updated as their statuses change):

> Chiou, K. L., Bergey, C. M., Burrell, A. S., Disotell, T. R., Rogers, J., Jolly, C. J., & Phillips-Conroy, J. E. Genomic signatures of extreme body size divergence in baboons.
> 
> Chiou, K. L., Bergey, C. M., Burrell, A. S., Disotell, T. R., Rogers, J., Jolly, C. J., & Phillips-Conroy, J. E. Genome-wide ancestry and introgression in a Zambian baboon hybrid zone.

Some scripts in this pipeline borrow heavily from the [awash-dopamine](https://github.com/bergeycm/awash-dopamine) repository, written by coauthor [Christina Bergey](http://christinabergey.com).

The analysis pipeline here begins with one pair of FASTQ files for each animal in the study. This should be most convenient for those wishing to replicate this analysis, as the files are the same as the ones available on [SRA](https://www.ncbi.nlm.nih.gov/sra).

For those wishing to recreate the steps for demultiplexing and, in some cases, merging libraries across independent sequencing runs, raw reads are available upon request (output from [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)). Demultiplexing code is available in the [demultiplex](https://github.com/kchiou/demultiplex) repository, which was originally written by [Christina Bergey](http://christinabergey.com).

The code throughout is documented as though all code is run line by line, with all programs available in the `$PATH`. In reality, the majority of the code was run on the [NYU high performance computing](https://wikis.nyu.edu/display/NYUHPC/Clusters) cluster, particularly those steps requiring extra time, memory, or processors. Most of the programs were available as preinstalled modules maintained by the HPC staff. Here, the modules are omitted and instead software links and version numbers are provided to help users find the appropriate versions. Similarly, for [R](https://cran.r-project.org) packages, you will need to check each script and install required libraries using `install.packages` or [Bioconductor](https://bioconductor.org/install).

If you find yourself unable to recreate any of the steps in this pipeline, start by checking version numbers and consulting the documentation for the relevant programs. If you still have questions, feel free to contact me and I will do my best to help.

Contact:
Kenneth Chiou ([website](http://kennychiou.com))
kchiou [at] uw [dot] edu

## Table of Contents

* [Mapping and genotyping pipeline](#mapping-and-genotyping-with-the-ngs-map-pipeline)

* [Analysis pipeline](#analysis-pipeline)

	* [SNP filtering](#snp-filtering)

	* [Exploratory data visualization](#exploratory-data-visualization)

	* [Ancestry inference](#ancestry-inference)

	* [Exploratory population genetic analyses](#exploratory-population-genetic-analyses)

	* [Genome scan for selection](#genome-scan-for-selection)

	* [Genomic cline analysis](#genomic-cline-analysis)

	* [Gene Ontology and PANTHER pathway enrichment](#gene-ontology-and-panther-pathway-enrichment)

## Mapping and genotyping with the [NGS-map](https://github.com/bergeycm/NGS-map) pipeline

To begin, ensure that the [NGS-map](https://github.com/bergeycm/NGS-map) repository is up-to-date (written by [Christina Bergey](http://christinabergey.com)).

```
git submodule update NGS-map
```

Reads downloaded from [SRA](https://www.ncbi.nlm.nih.gov/sra) (BioProject PRJNA486659) will have different file names from those used here. To ensure that the pipeline works properly, you may need to rename the reads to those indicated in the table below. Once finished, place reads in the `NGS-map/data` folder.

| SRA accession | Read 1 filename            | Read 2 filename            |
| ------------- | -------------------------- | -------------------------- |
| SRR7717351    | BZ06\_218.read1.fastq      | BZ06\_218.read2.fastq      |
| SRR7717350    | BZ06\_220.read1.fastq      | BZ06\_220.read2.fastq      |
| SRR7717353    | BZ06\_221.read1.fastq      | BZ06\_221.read2.fastq      |
| SRR7717352    | BZ06\_224.read1.fastq      | BZ06\_224.read2.fastq      |
| SRR7717355    | BZ06\_225.read1.fastq      | BZ06\_225.read2.fastq      |
| SRR7717354    | BZ06\_227.read1.fastq      | BZ06\_227.read2.fastq      |
| SRR7717357    | BZ07\_001.read1.fastq      | BZ07\_001.read2.fastq      |
| SRR7717356    | BZ07\_004.read1.fastq      | BZ07\_004.read2.fastq      |
| SRR7717359    | BZ07\_005.read1.fastq      | BZ07\_005.read2.fastq      |
| SRR7717358    | BZ07\_007.read1.fastq      | BZ07\_007.read2.fastq      |
| SRR7717328    | BZ07\_029.read1.fastq      | BZ07\_029.read2.fastq      |
| SRR7717327    | BZ07\_030.read1.fastq      | BZ07\_030.read2.fastq      |
| SRR7717326    | BZ07\_032.read1.fastq      | BZ07\_032.read2.fastq      |
| SRR7717325    | BZ07\_034.read1.fastq      | BZ07\_034.read2.fastq      |
| SRR7717332    | BZ07\_035.read1.fastq      | BZ07\_035.read2.fastq      |
| SRR7717331    | BZ07\_039.read1.fastq      | BZ07\_039.read2.fastq      |
| SRR7717330    | BZ07\_041.read1.fastq      | BZ07\_041.read2.fastq      |
| SRR7717329    | BZ07\_042.read1.fastq      | BZ07\_042.read2.fastq      |
| SRR7717324    | BZ07\_045.read1.fastq      | BZ07\_045.read2.fastq      |
| SRR7717323    | BZ07\_047.read1.fastq      | BZ07\_047.read2.fastq      |
| SRR7717395    | BZ07\_100.read1.fastq      | BZ07\_100.read2.fastq      |
| SRR7717396    | BZ11\_001.read1.fastq      | BZ11\_001.read2.fastq      |
| SRR7717393    | BZ11\_002.read1.fastq      | BZ11\_002.read2.fastq      |
| SRR7717394    | BZ11\_003.read1.fastq      | BZ11\_003.read2.fastq      |
| SRR7717399    | BZ11\_004.read1.fastq      | BZ11\_004.read2.fastq      |
| SRR7717400    | BZ11\_005.read1.fastq      | BZ11\_005.read2.fastq      |
| SRR7717397    | BZ11\_006.read1.fastq      | BZ11\_006.read2.fastq      |
| SRR7717398    | BZ11\_007.read1.fastq      | BZ11\_007.read2.fastq      |
| SRR7717401    | BZ11\_008.read1.fastq      | BZ11\_008.read2.fastq      |
| SRR7717402    | BZ11\_009.read1.fastq      | BZ11\_009.read2.fastq      |
| SRR7717384    | BZ11\_010.read1.fastq      | BZ11\_010.read2.fastq      |
| SRR7717383    | BZ11\_011.read1.fastq      | BZ11\_011.read2.fastq      |
| SRR7717386    | BZ11\_012.read1.fastq      | BZ11\_012.read2.fastq      |
| SRR7717385    | BZ11\_013.read1.fastq      | BZ11\_013.read2.fastq      |
| SRR7717388    | BZ11\_014.read1.fastq      | BZ11\_014.read2.fastq      |
| SRR7717387    | BZ11\_015.read1.fastq      | BZ11\_015.read2.fastq      |
| SRR7717390    | BZ11\_016.read1.fastq      | BZ11\_016.read2.fastq      |
| SRR7717389    | BZ11\_017.read1.fastq      | BZ11\_017.read2.fastq      |
| SRR7717392    | BZ11\_018.read1.fastq      | BZ11\_018.read2.fastq      |
| SRR7717391    | BZ11\_019.read1.fastq      | BZ11\_019.read2.fastq      |
| SRR7717293    | BZ11\_020.read1.fastq      | BZ11\_020.read2.fastq      |
| SRR7717294    | BZ11\_021.read1.fastq      | BZ11\_021.read2.fastq      |
| SRR7717295    | BZ11\_022.read1.fastq      | BZ11\_022.read2.fastq      |
| SRR7717296    | BZ11\_023.read1.fastq      | BZ11\_023.read2.fastq      |
| SRR7717297    | BZ11\_024.read1.fastq      | BZ11\_024.read2.fastq      |
| SRR7717298    | BZ11\_025.read1.fastq      | BZ11\_025.read2.fastq      |
| SRR7717299    | BZ11\_026.read1.fastq      | BZ11\_026.read2.fastq      |
| SRR7717300    | BZ11\_028.read1.fastq      | BZ11\_028.read2.fastq      |
| SRR7717301    | BZ11\_029.read1.fastq      | BZ11\_029.read2.fastq      |
| SRR7717302    | BZ11\_030.read1.fastq      | BZ11\_030.read2.fastq      |
| SRR7717281    | BZ11\_031.read1.fastq      | BZ11\_031.read2.fastq      |
| SRR7717280    | BZ11\_032.read1.fastq      | BZ11\_032.read2.fastq      |
| SRR7717279    | BZ11\_033.read1.fastq      | BZ11\_033.read2.fastq      |
| SRR7717278    | BZ11\_034.read1.fastq      | BZ11\_034.read2.fastq      |
| SRR7717277    | BZ11\_035.read1.fastq      | BZ11\_035.read2.fastq      |
| SRR7717276    | BZ11\_036.read1.fastq      | BZ11\_036.read2.fastq      |
| SRR7717275    | BZ11\_037.read1.fastq      | BZ11\_037.read2.fastq      |
| SRR7717274    | BZ11\_038.read1.fastq      | BZ11\_038.read2.fastq      |
| SRR7717283    | BZ11\_039.read1.fastq      | BZ11\_039.read2.fastq      |
| SRR7717282    | BZ11\_040.read1.fastq      | BZ11\_040.read2.fastq      |
| SRR7717321    | BZ11\_041.read1.fastq      | BZ11\_041.read2.fastq      |
| SRR7717322    | BZ11\_042.read1.fastq      | BZ11\_042.read2.fastq      |
| SRR7717319    | BZ11\_043.read1.fastq      | BZ11\_043.read2.fastq      |
| SRR7717320    | BZ11\_045.read1.fastq      | BZ11\_045.read2.fastq      |
| SRR7717317    | BZ11\_046.read1.fastq      | BZ11\_046.read2.fastq      |
| SRR7717318    | BZ11\_047.read1.fastq      | BZ11\_047.read2.fastq      |
| SRR7717315    | BZ11\_048.read1.fastq      | BZ11\_048.read2.fastq      |
| SRR7717316    | BZ11\_050.read1.fastq      | BZ11\_050.read2.fastq      |
| SRR7717313    | BZ11\_051.read1.fastq      | BZ11\_051.read2.fastq      |
| SRR7717314    | BZ11\_052.read1.fastq      | BZ11\_052.read2.fastq      |
| SRR7717310    | BZ11\_053.read1.fastq      | BZ11\_053.read2.fastq      |
| SRR7717309    | BZ11\_054.read1.fastq      | BZ11\_054.read2.fastq      |
| SRR7717312    | BZ11\_056.read1.fastq      | BZ11\_056.read2.fastq      |
| SRR7717311    | BZ11\_057.read1.fastq      | BZ11\_057.read2.fastq      |
| SRR7717306    | BZ11\_058.read1.fastq      | BZ11\_058.read2.fastq      |
| SRR7717305    | BZ11\_059.read1.fastq      | BZ11\_059.read2.fastq      |
| SRR7717308    | BZ11\_061.read1.fastq      | BZ11\_061.read2.fastq      |
| SRR7717307    | BZ11\_062.read1.fastq      | BZ11\_062.read2.fastq      |
| SRR7717304    | BZ11\_063.read1.fastq      | BZ11\_063.read2.fastq      |
| SRR7717303    | BZ11\_064.read1.fastq      | BZ11\_064.read2.fastq      |
| SRR7717368    | BZ11\_065.read1.fastq      | BZ11\_065.read2.fastq      |
| SRR7717369    | BZ11\_066.read1.fastq      | BZ11\_066.read2.fastq      |
| SRR7717370    | BZ11\_067.read1.fastq      | BZ11\_067.read2.fastq      |
| SRR7717371    | BZ11\_068.read1.fastq      | BZ11\_068.read2.fastq      |
| SRR7717364    | BZ11\_069.read1.fastq      | BZ11\_069.read2.fastq      |
| SRR7717365    | BZ11\_070.read1.fastq      | BZ11\_070.read2.fastq      |
| SRR7717366    | BZ11\_071.read1.fastq      | BZ11\_071.read2.fastq      |
| SRR7717367    | BZ11\_072.read1.fastq      | BZ11\_072.read2.fastq      |
| SRR7717361    | BZ11\_073.read1.fastq      | BZ11\_073.read2.fastq      |
| SRR7717362    | BZ11\_074.read1.fastq      | BZ11\_074.read2.fastq      |
| SRR7717340    | BZ11\_075.read1.fastq      | BZ11\_075.read2.fastq      |
| SRR7717334    | BZ11\_076.read1.fastq      | BZ11\_076.read2.fastq      |
| SRR7717335    | BZ12\_001.read1.fastq      | BZ12\_001.read2.fastq      |
| SRR7717333    | BZ12\_002.read1.fastq      | BZ12\_002.read2.fastq      |
| SRR7717341    | BZ12\_003.read1.fastq      | BZ12\_003.read2.fastq      |
| SRR7717374    | BZ12\_004.read1.fastq      | BZ12\_004.read2.fastq      |
| SRR7717339    | BZ12\_005.read1.fastq      | BZ12\_005.read2.fastq      |
| SRR7717338    | BZ12\_006.read1.fastq      | BZ12\_006.read2.fastq      |
| SRR7717363    | BZ12\_007.read1.fastq      | BZ12\_007.read2.fastq      |
| SRR7717360    | BZ12\_008.read1.fastq      | BZ12\_008.read2.fastq      |
| SRR7717344    | BZ12\_009.read1.fastq      | BZ12\_009.read2.fastq      |
| SRR7717345    | BZ12\_010.read1.fastq      | BZ12\_010.read2.fastq      |
| SRR7717342    | BZ12\_011.read1.fastq      | BZ12\_011.read2.fastq      |
| SRR7717343    | BZ12\_012.read1.fastq      | BZ12\_012.read2.fastq      |
| SRR7717348    | BZ12\_030.read1.fastq      | BZ12\_030.read2.fastq      |
| SRR7717349    | BZ12\_031.read1.fastq      | BZ12\_031.read2.fastq      |
| SRR7717346    | BZ12\_032.read1.fastq      | BZ12\_032.read2.fastq      |
| SRR7717347    | BZ12\_033.read1.fastq      | BZ12\_033.read2.fastq      |
| SRR7717336    | Chiou\_14\_001.read1.fastq | Chiou\_14\_001.read2.fastq |
| SRR7717337    | Chiou\_14\_003.read1.fastq | Chiou\_14\_003.read2.fastq |
| SRR7717376    | Chiou\_14\_004.read1.fastq | Chiou\_14\_004.read2.fastq |
| SRR7717375    | Chiou\_14\_005.read1.fastq | Chiou\_14\_005.read2.fastq |
| SRR7717378    | Chiou\_14\_030.read1.fastq | Chiou\_14\_030.read2.fastq |
| SRR7717377    | Chiou\_14\_036.read1.fastq | Chiou\_14\_036.read2.fastq |
| SRR7717380    | Chiou\_14\_039.read1.fastq | Chiou\_14\_039.read2.fastq |
| SRR7717379    | Chiou\_14\_041.read1.fastq | Chiou\_14\_041.read2.fastq |
| SRR7717382    | Chiou\_14\_042.read1.fastq | Chiou\_14\_042.read2.fastq |
| SRR7717381    | Chiou\_14\_044.read1.fastq | Chiou\_14\_044.read2.fastq |
| SRR7717373    | Chiou\_14\_050.read1.fastq | Chiou\_14\_050.read2.fastq |
| SRR7717372    | Chiou\_14\_054.read1.fastq | Chiou\_14\_054.read2.fastq |
| SRR7717284    | Chiou\_14\_056.read1.fastq | Chiou\_14\_056.read2.fastq |
| SRR7717285    | Chiou\_14\_057.read1.fastq | Chiou\_14\_057.read2.fastq |
| SRR7717286    | Chiou\_14\_058.read1.fastq | Chiou\_14\_058.read2.fastq |
| SRR7717287    | Chiou\_14\_059.read1.fastq | Chiou\_14\_059.read2.fastq |
| SRR7717288    | Chiou\_14\_065.read1.fastq | Chiou\_14\_065.read2.fastq |
| SRR7717289    | Chiou\_14\_069.read1.fastq | Chiou\_14\_069.read2.fastq |
| SRR7717290    | Chiou\_15\_003.read1.fastq | Chiou\_15\_003.read2.fastq |
| SRR7717291    | Chiou\_15\_004.read1.fastq | Chiou\_15\_004.read2.fastq |
| SRR7717292    | Chiou\_15\_005.read1.fastq | Chiou\_15\_005.read2.fastq |

Proceed to map samples and call and filter variants according to the [NGS-map](https://github.com/bergeycm/NGS-map) pipeline. The following outline may be helpful but will not work as is--some components of [NGS-map](https://github.com/bergeycm/NGS-map) will need to be manually configured (particularly paths to software).

1. Prepare sample list and save to a file.
```
ls -1 NGS-map/data/*.fastq | \
	sed 's/NGS-map\/data\/\([A-Za-z0-9_]*\)\.read[1-2].fastq/\1/' | \
	sort -u > NGS-map/data/individual_list.txt
```

2. Change to `NGS-map` directory.

```
cd NGS-map
```

3. Download and index the baboon reference genome (Panu2.0).

```
genomes/download_papAnu2.sh
scripts/index_genome.sh genomes/papAnu2/papAnu2.fa
```

4. Map samples to the reference genome using [BWA mem](http://bio-bwa.sourceforge.net/) (v0.7.15).

```
qsub -t 1-129 pbs/call_make.pbs
```

5. Call genotypes using [GATK](https://software.broadinstitute.org/gatk/) (v3.5.0).

```
qsub -t 1-20 pbs/call_gatk_genotyper.pbs
```

6. Filter genotype calls using [GATK](https://software.broadinstitute.org/gatk/) (v3.5.0).

```
qsub -t 1-20 pbs/filter_gatk_snps.pbs
```

7. Change to directory containing variant calls.

```
cd baboon_snps
```

8. Concatenate SNP calls from autosomes using [VCFtools](http://vcftools.sourceforge.net/) (v0.1.14).

```
vcf-concat $(echo chr{1..20}.pass.snp.vcf) > kafue.pass.snp.noX.vcf
```

9. Change to base project directory.

```
cd ../..
```

Following the process above, transfer results files to the `data` folder (at the base level of this repository)

```
mv NGS-map/baboon_snps/kafue.pass.snp.noX.vcf data
```

## Analysis pipeline

### SNP filtering

* Download [Repeatmasker](www.repeatmasker.org/) and [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) files for the baboon reference genome (Panu2.0).

```
scripts/download_repeats.sh
```

* Convert VCF files to [PLINK](http://zzz.bwh.harvard.edu/plink) PED and MAP files using [VCFtools](http://vcftools.sourceforge.net/) (v0.1.14).

```
scripts/vcf_to_ped.sh kafue.pass.snp.noX
```

* Calculate stats on missingness using [PLINK](http://zzz.bwh.harvard.edu/plink) (v1.90).

```
plink --noweb --double-id --vcf data/kafue.pass.snp.noX.vcf --missing --out reports/kafue.pass.snp.noX
```

* Find repeats based on the [Repeatmasker](www.repeatmasker.org/) and [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) files and the [PLINK](http://zzz.bwh.harvard.edu/plink) `.map` file generated in the previous step. [BEDtools](https://bedtools.readthedocs.io/en/latest/) (v2.26.0) must be installed and in the `$PATH`.

```
scripts/repeat_finder.pl
```

* Filter dataset to remove repeats identified in the previous step as well as samples and loci with excessive missingness using [PLINK](http://zzz.bwh.harvard.edu/plink) (v1.90). This will produce two datasets: the larger `results/all.cleaned.*` dataset and the more stringent `results/strict.cleaned.*` dataset (the latter is not used in subsequent steps).

```
scripts/make_clean_dataset.sh
```

* Calculate SNPs in linkage disequilibrium and remove them from the dataset using [PLINK](http://zzz.bwh.harvard.edu/plink) (v1.90). This will produce the `results/all.cleaned.LDpruned.*` dataset.

```
scripts/find_and_remove_ld.sh
```

### Exploratory data visualization

* Perform multidimensional scaling and principal components analyses using [PLINK](http://zzz.bwh.harvard.edu/plink) (v1.90). The scripts also plot respective results using [R](https://cran.r-project.org) (v3.4.2).

```
scripts/mds_plot.sh
scripts/pca_plot.sh
```

* Visualize library mapping statistics from individuals passing filters using [R](https://cran.r-project.org) (v3.4.2). This step will read the `results/mds.txt` file generated in the previous step.

```
scripts/plot_mapping.R
```

### Ancestry inference

* Raster and vector imagery for Zambia and Kafue National Park are used for subsequent spatial analyses or plotting. The following datasets were used to prepare a custom [ASCII](http://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/esri-ascii-raster-format.htm)-formatted grid file for this project using [R](https://cran.r-project.org) (v3.4.2):

	* The [bio1](http://www.worldclim.org/current) bioclimatic dataset, obtained from [WorldClim](http://www.worldclim.org/version1) and stored in [ASCII](http://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/esri-ascii-raster-format.htm) format. File was moved to: `data/bio01.asc`.

	* Zambia country administrative borders, obtained from [DIVA-GIS](http://diva-gis.org/gdata). The specific file used is `ZMB_adm0.shp`, renamed to: `data/zambia_borders.shp`.

	* Park boundaries for Kafue National Park, kindly shared by the Zambia Wildlife Authority (now the [Department of National Parks and Wildlife](https://www.zambiatourism.com/about-zambia/conservation/department-of-national-parks-wildlife/) of Zambia. The particular file is not mine to share, but a perfectly suitable [alternative](https://www.protectedplanet.net/1085) is available for download at [Protected Planet](https://www.protectedplanet.net). The file was moved to: `data/kafue_parkboundary.shp`.

```
scripts/make_asc_file.R
```

* Convert VCF file to BED file using [PLINK](http://zzz.bwh.harvard.edu/plink) (v1.90) to prepare dataset for ancestry inference. Will generate `results/all.cleaned.LDpruned.vcf`.

```
plink --bfile results/all.cleaned.LDpruned --recode vcf --out results/all.cleaned.LDpruned
mv results/all.cleaned.LDpruned.log reports/bed_to_vcf.log
```

* Write spatial coordinates into file that is line-matched to the VCF data using [R](https://cran.r-project.org) (v3.4.2). Will generate `data/all.cleaned.LDpruned.coord`.

```
scripts/make_coord_file.R
```

* Infer ancestry using [ADMIXTURE](http://software.genetics.ucla.edu/admixture/) (v1.3.0) in cross-validation mode testing values of *K* from 1 to 10. Plot cross-validated errors using [R](https://cran.r-project.org) (v3.4.2).

```
scripts/run_admixture_cv.sh
scripts/plot_admixture_cv_error.R
```

* Run [ADMIXTURE](http://software.genetics.ucla.edu/admixture/) (v1.3.0) on the optimal value of *K* (*K* = 2).

```
scripts/run_admixture_k.sh 2
```

* Generate a variety of ancestry plots using [R](https://cran.r-project.org) (v3.4.2), including comparisons to multidimensional scaling analyses done earlier. The script uses functions from `POPSutilities.R`, which is not provided here and must be downloaded ([link](http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R)) and moved to the `scripts` folder. 

```
scripts/plot_admixture_results.R
```

* Generate ancestry plots from [Jolly et al. 2011](http://doi.org/10.1002/ajp.20896). The script uses the data set `data/sex_chr_genotypes.txt`, which is not provided here but can be requested from [Cliff Jolly](https://wp.nyu.edu/csho/people/faculty/clifford_jolly/).

```
scripts/plot_jolly_et_al_2011.R
```

* Compare autosomal and sex chromosome ancestries.

```
scripts/test_autosomes_vs_uniparental_markers.R
```

### Exploratory population genetic analyses

* Partition populations based on ancestry results using [PLINK](http://zzz.bwh.harvard.edu/plink) (v1.90) and prepare and convert requisite files for population genomic scans using [PLINK/SEQ](https://atgu.mgh.harvard.edu/plinkseq/) (v0.10), [VCFtools](http://vcftools.sourceforge.net/) (v0.1.14), [htslib](http://www.htslib.org/doc/tabix.html) `tabix` (v1.4.1), and [R](https://cran.r-project.org) (v3.4.2).

```
scripts/prep_popgen.sh
```

* Calculate a variety of within-population and between-population (<em>F<sub>ST</sub></em>) statistics using [VCFtools](http://vcftools.sourceforge.net/) (v0.1.14).

```
scripts/popgen_stats.sh
```

* Generate exploratory plot of inbreeding according to ancestry using [R](https://cran.r-project.org) (v3.4.2).

```
scripts/plot_inbreeding.R
```

* Make maps for spatial visualization of sample locales and geospatial ancestry using [R](https://cran.r-project.org) (v3.4.2).

```
scripts/map_samples.R
```

* Prepare for genome scans by downloading baboon and macaque annotations, finding open reading frames (ORFs), mapping baboon ORFs to macaque proteins, and finalizing baboon PANTHER annotations derived from homology. [sesbio](https://github.com/sestaton/sesbio) (v2016-10-18), [HMMER2GO](https://github.com/sestaton/HMMER2GO) (v0.17.2), [HMMER](http://hmmer.org) (v3.1b2), [SAMtools](http://samtools.sourceforge.net/) (v1.6), and [R](https://cran.r-project.org) (v3.4.2) are used and should be in the `$PATH`, along with the dependency [EMBOSS](http://emboss.sourceforge.net) (v6.6.0).

### Genome scan for selection

```
scripts/prep_genome_scans.sh
```

* Assign and average <em>F<sub>ST</sub></em> values across genes and conduct permutation test to build empirical <em>F<sub>ST</sub></em> null distribution using [R](https://cran.r-project.org) (v3.4.2).

```
scripts/scan_fst.R
```

* Some *p* values from the above procedure might be so close to the threshold of 0.05 that we cannot confidently determine their significance. For these sets of genes, we can run more permutations to augment ("boost") their empirical <em>F<sub>ST</sub></em> null distributions to bring the total up one order of magnitude (e.g., from 1 million to 10 million permutations). The procedure will make use of the `.RData` files in the `checkpoints` folder. To make this more efficient, we can parallelize the process using the following script. The script uses a modified R script `scripts/scan_fst_boost.R` and passes as argument two variables:

	1. `$reps` is automatically calculated and shows how many permutations have already been run (on a log10 scale)

	2. `$1` (may be used for [GNU parallel](https://www.gnu.org/software/parallel/), [PBS](https://www.pbspro.org/) [replace `$1` with `$PBS_ARRAYID`], or [Slurm](https://slurm.schedmd.com) [replace `$1` with `$SLURM_ARRAY_TASK_ID`]) is an arbitrary integer that ensures that result files do not conflict.

```
#!/bin/sh

reps=$(ls checkpoints/pval.perm.1e*.RData | sed -e 's/checkpoints\/pval.perm.1e\([0-9]*\)\.RData/\1/' | sort -rn | head -n 1)

scripts/scan_fst_boost.R $reps $1
```

* If the boost procedure described above was run, consolidate the results using [R](https://cran.r-project.org) (v3.4.2) to generate the final empirical null distribution.

```
scripts/scan_fst_aggregate.R
```

### Genomic cline analysis

* Convert genotypes into the genotype files require for [bgc](https://sites.google.com/site/bgcsoftware/) (v1.03). The script takes two VCF files as arguments because the first includes read count information but also extraneous SNPs (those not passing filters) while the second contains SNPs passing filters but no read counts. The conversion script merges information from the two to produce read counts for SNPs passing filters in the proper format.

```
scripts/vcf_to_bgc.R data/kafue.pass.snp.noX.vcf results/all.cleaned.LDpruned.vcf
```

* Run [bgc](https://sites.google.com/site/bgcsoftware/) (v1.03). The script below can be run in parallel (5 parallel runs were done in this case) by passing an arbitrary argument `$1` (may be used for [GNU parallel](https://www.gnu.org/software/parallel/), [PBS](https://www.pbspro.org/) [replace `$1` with `$PBS_ARRAYID`], or [Slurm](https://slurm.schedmd.com) [replace `$1` with `$SLURM_ARRAY_TASK_ID`]). [GSL](https://www.gnu.org/software/gsl/doc/html) (v1.16) and [HDF5](https://portal.hdfgroup.org/display/support) (v1.8.17) are required as dependencies, and the `bgc` and `estpost` programs must be properly compiled and in the `$PATH`. The stdout of one chain was saved to the file `reports/bgc_err.txt`. Running on a remote machine is essentially required as each run will take a long time to converge and complete.

```
#!/bin/sh

scripts/run_bgc.sh $1
```

* Aggregate posterior MCMC chains across the 5 [bgc](https://sites.google.com/site/bgcsoftware/) (v1.03) runs.

```
scripts/run_bgc_aggregate.sh
```

* Visualize genomic cline results and compare to other ancestry results.

```
scripts/plot_bgc_results.R
```

* Assign genomic cline parameter estimates to genes based on overlap.

```
scripts/run_bgc_genes.R
```

* Plot results from both <em>F<sub>ST</sub></em> and genomic cline scans.

```
scripts/plot_all_scans.R
```

* Test for difference in beta cline parameter estimates between genic and nongenic SNPs.

```
scripts/run_bgc_genic_vs_nongenic.R
```

### Gene Ontology and PANTHER pathway enrichment

* Test for enrichment of [Gene Ontology](http://www.geneontology.org) terms and [PANTHER](http://www.pantherdb.org/) pathways for both <em>F<sub>ST</sub></em> and genomic cline analyses.

```
scripts/test_enrichment_all.R
```
