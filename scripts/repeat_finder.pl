#!/usr/bin/perl

use strict;
use warnings;

use File::Copy qw(copy);

#	$ head data/papAnu2.fa.out 
#	   SW  perc perc perc  query    position in query           matching  repeat      ...
#	score  div. del. ins.  sequence  begin     end    (left)    repeat    class/family...
#	
#	 2891  17.4  5.9  6.3  chr1          1     600 (220367099) +  L1MB7     LINE/L1   ...
#	 1679  17.6  5.3  0.0  chr1        601     884 (220366815) +  AluJo     SINE/Alu  ...
#	 2021  13.7  0.3  5.2  chr1        885    1212 (220366487) +  AluSx1    SINE/Alu  ...
#	
#	$ head data/papAnu2.trf.bed 
#	chr1	44934	45014	trf	7	11.1	7	76	5	79	3	18	51	26	1.63	CTG...
#	chr1	51643	51683	trf	10	3.9	10	84	15	55	20	10	70	1.16	GGGAGGGGCA
#	chr1	62440	62505	trf	8	8.2	8	76	10	53	1	58	0	40	1.07	CTCCCTTC
#	
#	$ head results/all.cleaned.map 
#	1	52883	0	52883
#	1	53024	0	53024
#	1	93012	0	93012
#	1	126058	0	126058
#	1	126069	0	126069

my $INPUT_PREFIX = "kafue.pass.snp.noX";                        # Change this to toggle dataset

my $rm_out_file  = "data/papAnu2.fa.out";                       # $ENV{'RM_OUT'}
my $trf_out_file = "data/papAnu2.trf.bed";                      # $ENV{'TRF_OUT'}

my $map_file_full = "data/${INPUT_PREFIX}.map";

# File to contain excluded regions
my $exclude_bed = "data/excluded_repetitive_regions.bed";
my $exclude_txt = $exclude_bed;
$exclude_txt =~ s/\.bed/.plink.txt/;

# ----------------------------------------------------------------------------------------
# --- Turn map files into BEDs
# ----------------------------------------------------------------------------------------

 open(MAP_FULL, "<$map_file_full")
 	or die "ERROR: Could not open map file [$map_file_full].\n";
 
 open(MAP_FULL_BED, ">${map_file_full}.bed")
 	or die "ERROR: Could not open new BED file for map.\n";
 
 while (<MAP_FULL>) {
 
 	my @map_info = split /\t/;
 	
 	my $chr   = "chr" . $map_info[0];
 	my $start = $map_info[1];
 	$start =~ s/chr[X-Y0-9]*://g;
 	my $end   = $start + 1;
 
 	print MAP_FULL_BED $chr . "\t" . $start . "\t" . $end . "\n";
 }
 
 close MAP_FULL;
 close MAP_FULL_BED;

# ----------------------------------------------------------------------------------------
# --- Turn RM out file into BED
# ----------------------------------------------------------------------------------------

open(RM, "<$rm_out_file")
	or die "ERROR: Could not open RM out file [$rm_out_file].\n";

open(RM_BED, ">${rm_out_file}.bed")
	or die "ERROR: Could not open new BED file for RM.\n";

for (my $i = 0; $i < 3; $i++) {
	my $header = <RM>;
}

while (<RM>) {

	# Remove leading space
	$_ =~ s/^\s+//g;
	
	my @rm_info = split / +/;
	
	my $chr   = $rm_info[4];
	my $start = $rm_info[5];
	my $end   = $rm_info[6] + 1;
	
	if ($start > $end) {
		my $start_backup = $start;
		$start = $end;
		$end = $start_backup;
	}
	
	print RM_BED $chr . "\t" . $start . "\t" . $end . "\n";
}

close RM;
close RM_BED;

# ----------------------------------------------------------------------------------------
# --- Make BED of all regions to exclude; Sort all the things
# ----------------------------------------------------------------------------------------

my $rm_out_file_sort   = "${rm_out_file}.bed";
my $trf_out_file_sort  = $trf_out_file;
my $map_file_full_sort = "${map_file_full}.bed";

$rm_out_file_sort   =~ s/\.bed/.sort.bed/;
$trf_out_file_sort  =~ s/\.bed/.sort.bed/;
$map_file_full_sort =~ s/\.bed/.sort.bed/;

# Sort RM BED and remove non-chromosomal stuff
system("sort -k1,1 -k2,2n ${rm_out_file}.bed | grep '^chr'> $rm_out_file_sort");

# Sort TRF BED and remove non-chromosomal stuff (and extra columns)
system("sort -k1,1 -k2,2n $trf_out_file | cut -f1-3 | grep '^chr'> $trf_out_file_sort");

# Sort MAP BED while we're here
system("sort -k1,1 -k2,2n ${map_file_full}.bed | grep '^chr'> $map_file_full_sort");

# Combine RM and TRF regions into one BED
my $exclusion_regions_cmd = "cat ${trf_out_file_sort} ${rm_out_file_sort} | ";
$exclusion_regions_cmd .= "sort -k1,1 -k2,2n | ";
$exclusion_regions_cmd .= "bedtools merge -i stdin > $exclude_bed";
system($exclusion_regions_cmd);

# ----------------------------------------------------------------------------------------
# --- Intersect repetitive regions and SNPs
# ----------------------------------------------------------------------------------------

# Get count of SNPs before filtering
my $start_count_full = `wc -l $map_file_full_sort`;

# From SNP set subtract RM regions.
my $full_sub_rm_cmd = "bedtools intersect -v -a $map_file_full_sort ";
$full_sub_rm_cmd .= "-b $rm_out_file_sort | wc -l";
my $full_sub_rm = `$full_sub_rm_cmd`;

# From SNP set subtract TRF regions.
my $full_sub_trf_cmd = "bedtools intersect -v -a $map_file_full_sort ";
$full_sub_trf_cmd .= "-b $trf_out_file_sort | wc -l";
my $full_sub_trf = `$full_sub_trf_cmd`;

# From SNP set subtract RM and TRF regions.
my $full_sub_rm_trf_cmd = "bedtools intersect -v -a $map_file_full_sort ";
$full_sub_rm_trf_cmd .= "-b $exclude_bed | wc -l";
my $full_sub_rm_trf = `$full_sub_rm_trf_cmd`;

# ----------------------------------------------------------------------------------------
# --- Write info on repetitive regions excluded to log and STDERR
# ----------------------------------------------------------------------------------------

my $info_string = "Count of SNPs before filtering:\t$start_count_full\n";
$info_string .= "Count of SNPs less RM regions:\t$full_sub_rm\n";
$info_string .= "Count of SNPs less TRF regions:\t$full_sub_trf\n";
$info_string .= "Count of SNPs less RM+TRF regions:\t$full_sub_rm_trf\n";

open (LOG, ">reports/repetitive_region_filtering.txt")
	or die "ERROR: Could not open report file.\n";

print LOG $info_string;
print STDERR $info_string;

close LOG;

# ----------------------------------------------------------------------------------------
# --- Convert excluded BED file to range file for plink's exclude command
# ----------------------------------------------------------------------------------------

# Range file should look like this:
#	2 30000000 35000000  R1
#	2 60000000 62000000  R2
#	X 10000000 20000000  R3

open(EX_BED, "<$exclude_bed")
	or die "ERROR: Could not open BED file of excluded regions [$exclude_bed].\n";
	
open(EX_TXT, ">$exclude_txt")
	or die "ERROR: Could not open new file of excluded ranges.\n";

my $exclude_region_count = 1;

while (<EX_BED>) {

	my @bed_info = split /\s+/;
	my ($chr, $start, $end) = @bed_info;
	$chr =~ s/chr//;
	# BEDtools and plink handle ends of intervals differently
	$end--;
	print EX_TXT "$chr $start $end R${exclude_region_count}\n";
	$exclude_region_count++;
}

close EX_BED;
close EX_TXT;

exit;
