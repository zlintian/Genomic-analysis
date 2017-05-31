# Genomic-analysis

Currently there are five cohorts I have finished imputation:

BSGS, collected from Allan.

QLD_MND, collected from Beben.

CNH_MND, collected from Beben.

PD, collected from Costanza.

LBCW1, collected from Qian.

The genotype data can be found in each directory under /shares/compbio/PCTG/methylation/mQTL_project/#_cohort in Delta cluster.

There's a readme file in each cohort. You can read it and decide whether you would like to start with my imputed data, or QC and impute by yourself.

Imputed data are saved in /shares/compbio/PCTG/methylation/mQTL_project/Sange_imputed, and /shares/compbio/PCTG/methylation/mQTL_project/Michigan_imputed. They are under processing now.

Always check with me before you use the imputed data to make sure I have updated it to the latest I have.


########################Genotype preparation workflow###############

The following steps should fit most SNP data, but there can be special situations come up. 

0. define sample list. 
1. recode to ATGC if original data was coded with 1234
2. liftover to build 37 if in build 36 (this step can be omited because it can be done with the flip step, step6.
3. get the subset of data of selected individuals
4. pick the subset of chr1-23 if in original has 1-26. Only 1-23 can be imputed.
5. quality control. Parameter can be different.
6. flip the minus strand based on the strand information. Define which chip is the data from first.
7. convert the SNPs identifiers if a lot of non rs IDs exist.
8. separate the whole data into chromosomes if use Michigan imputation server.
9. fix the ref allele if use the Sanger imputation server.
10. zip the vcf file into vcf.gz
11. submit to imputation
12. download results, and transfer to plink format
13. combine the chromosomes into one file if neccessary.
14. Choose a proper software for GWAS


########################Quick search of Core Codes##################

	ssh username@delta.imb.uq.edu.au

	cd /shares/compbio/PCTG/methylation/mQTL_project/#_cohort/COHORT_genotype

	qsub -l select=1:mem=8G -l walltime=12:00:00 -I -q interactive

	for tar.gz files, use tar -xvf 

	for .tar files, use tar -xvf 

	for .gz files, use unzip

	./plink --bfile cohort --recode --make-bed --alleleACGT --out cohort_recoded

	python liftMap.py -m LBCW1_recoded.map -p LBCW1_recoded.ped -o LBCW1_recoded_lifted

	./plink --file LBCW1_recoded_lifted --keep genotype.IDs.txt --chr 1-23  --make-bed --out LBCW1_subset

	./plink --bfile cohort --mind 0.05 --geno 0.05 --hwe 0.000001 --maf 0.01  --make-bed --out Tian_cleand_cohort

	./plink --bfile Tian_cleand_bsgs --flip Human610Quadv1_B-b37-strand-v2/Human610-Quadv1_B-b37-v2_minus.strand  --recode  --make-bed  --out flip_Tian_cleand

	for i in $(seq 1 22); do ./plink --bfile idsed.fliped.QLD_MND  --chr $i --recode vcf --out QLD_MND-chr$i; done

	for i in $(seq 1 22); do vcf-sort QLD_MND-chr$i.vcf | bgzip -c > QLD_MND-chr$i.vcf.gz; done

	./plink --bfile idsed.fliped.QLD_MND  --exclude SNPs.in.hh.txt --chr 23 --recode vcf --out QLD_MND-chr23

#more to add


########################Important resourses#########################

##For tips using delta server:

http://igc-wiki.imb.uq.edu.au/igc/PBSPro%20usage

##All my data saved and all processing run in delta:

cd /shares/compbio/PCTG/methylation/mQTL_project/1_cohort/COHORT_genotype

plink: http://zzz.bwh.harvard.edu/plink/dataman.shtml#exclude

strand:http://www.well.ox.ac.uk/~wrayner/strand/

liftover: http://genome.sph.umich.edu/wiki/LiftOver

imputation: https://imputationserver.sph.umich.edu/start.html#!pages/help

HRC:

http://www.haplotype-reference-consortium.org/site

ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz


Book chapter:Bioinformatics Volume II: Structure, Function and Applications Second Edition

			Chapter 9 Analysis of Genome-Wide Association Data.
			
			Author: Allan F. McRea

##use an interactive node when necessary. 

qsub -l select=1:mem=8G -l walltime=12:00:00 -I -q interactive


########################Collect data from colleague#################

#better ask for the no-duplicate list too

#make a copy and name the folder as COHORT_genotype

##if get the compressed data

for tar.gz files, use tar -xvf 

for .tar files, use tar -xvf

for .gz files, use unzip


########################recode to ATGC if using 1234################

check the bim file, if using 1234, recode it.

	./plink --bfile cohort --recode --make-bed --alleleACGT --out cohort_recoded


:) Thanks to Qian.

########################Liftover####################################

check whether it's in build 37

if it's in build 36, using following method to lift it.

This step can be skiped if a minus strand flip will be carried out with update_build.sh


#web site:

http://genome.sph.umich.edu/wiki/LiftOver

Lift PLINK format

PLINK format usually referrs to .ped and .map files.


Method 1

We mainly use UCSC LiftOver binary tools to help lift over. We have a script liftMap.py, however, it is recommended to understand the job step by step:

(1) Convert .map to .bed file

By rearrange columns of .map file, we obtain a standard BED format file.

(2) LiftOver .bed file

Use method mentioned above to convert .bed file from one build to another.

(3) Convert lifted .bed file back to .map file

Rearrange column of .map file to obtain .bed file in the new build.

(4) Modify .ped file

.ped file have many column files. By convention, the first six columns are family_id, person_id, father_id, mother_id, sex, and phenotype. From the 7th column, there are two letters/digits representing a genotype at the certain marker. In step (2), as some genome positions cannot be lifted to the new version, we need to drop their corresponding columns from .ped file to keep consistency. You can use PLINK --exclude those snps, see Remove a subset of SNPs.

(5) (optionally) change the rs number in the .map file

Similar to the human reference build, dbSNP also have different versions. You may consider change rs number from the old dbSNP version to new dbSNP version depending on your needs. Such steps are described in Lift dbSNP rs numbers.

##Liftover to build 37 (hg19) use python script liftmap.py

##Download the three files:

liftMap.py

hg18ToHg19.over.chain

liftOver (a binary file that need to chmod u+x, and use it as ./liftOver)


	vi liftMap.py
	##find the following part and change the params['LIFTOVER_BIN'] and params['CHAIN']

	params = dict()
    params['LIFTOVER_BIN'] = './liftOver'
    params['OLD'] = fin
    params['CHAIN'] = 'hg18ToHg19.over.chain'
    params['NEW'] = fout
    params['UNLIFTED'] = fout + '.unlifted'

##See how to use it 
[tian]$ python liftMap.py 
usage: liftMap.py [-h] -m MAPFILE [-p PEDFILE] [-d DATFILE] -o PREFIX

liftMap.py: error: argument -m is required

##CAUTION: liftover don't recognize chr23,24. chr23 need to change to X. chr24 need to change to Y

	sed 's/^23/X/' LBCW1_recoded.map > LBCW1_recoded_Xchanged.map

	cp LBCW1_recoded.ped LBCW1_recoded_Xchanged.ped

##now we have both map and ped files ready.

##run the following command using qsub is better. it takes quite a few minutes.

	python liftMap.py -m LBCW1_recoded.map -p LBCW1_recoded.ped -o LBCW1_recoded_lifted

##the outputs:

	LBCW1_recoded_lifted.ped	
	LBCW1_recoded_lifted.map	
	LBCW1_recoded_lifted.bed	
	LBCW1_recoded_lifted.bed.unlifted	
	LBCW1_recoded.map.bed

LBCW1_recoded.map.bed is the BED format file which is reformatted from the plink bed file. Don't need it.

LBCW1_recoded_lifted.bed file is in BED format, we don't need it either. change the name into LBCW1_recoded_lifted.bed.format, so won't mass up as plink bed file.

LBCW1_recoded_lifted.map and LBCW1_recoded_lifted.ped are what we need.

check how many SNPs in the unlifted file. In my file, all of them are shown as "#Deleted in new".


:) Thanks to Yang.

########################look into the genotype data carefully####################

how many SNPs in the bim file.

how many individuals in the genotype data est. (belong to 666 families)

how many individuals in the methylation data set.

how many individuals in the information sheet

#are genotype data overlap with methylation data set?

how many of them are overlapped.

how many extra in methylation data set.

how many extra in genotyped data set.

awk '{print $1}' COHORT_subset.bim | uniq -c | less

Check how many families and whether there are multiple individuals in the same family.

########################Get the subset of genotype##################

Determine what samples to keep based on the methylation data.
#1. Use the ID and the matched family ID (which is the same of ID) in final.txt to form a table "genotype.IDs.txt".
#2. check the chromosomes included in the genotype data
#3. use plink to get the subset of genotype data.

	./plink --file LBCW1_recoded_lifted --keep genotype.IDs.txt --chr 1-23  --make-bed --out LBCW1_subset

#how many samples in raw
#how many samples to keep
#how many SNPs input
#how many SNPs kept

########################QC##########################################

	./plink --bfile cohort --mind 0.05 --geno 0.05 --hwe 0.000001 --maf 0.01  --make-bed --out Tian_cleand_cohort

#make a record of variables and samples removed and kept.

how many people input

how many SNPs input

how many people removed

how many SNPs removed as missing

how many SNPs removed due to hwe

how many SNPs removed due to low MAF

how many people kept

how many SNPs kept

Any other files than bam,bed,fam,log generated?

plink.hh file:List of heterozygous haploid genotypes (SNPs/individuals)

plink.irem	--mind	List of individuals removed for low genotyping

########################flip the minus strand########################

##Find the correct file for the chip used from the following webpage.

http://www.well.ox.ac.uk/~wrayner/strand/

http://www.well.ox.ac.uk/~wrayner/tools/

#go through the file, greb a few SNPs and double check with the bim file.

##CAUTION: look into the SNP identifiers, are there a lot of kgp? 
#go through the identifiers used in each file. and make sure.
#are they both uppercase or lower case?
#are there extra letters in the bim file or the strand file?

##1. grep the minus strand from the strand file

	grep -w '-' Human610-Quadv1_B-b37-v2.strand > Human610-Quadv1_B-b37-v2_minus.strand 

##CAUTION: Remember to use -w in grep. otherwise the SNPs with a - in its identifier will be grepped and flipped.

##2. use plink --flip function to flip the minus strand in bim file.

	./plink --bfile Tian_cleand_bsgs --flip Human610Quadv1_B-b37-strand-v2/Human610-Quadv1_B-b37-v2_minus.strand  --recode  --make-bed  --out flip_Tian_cleand

#go through the file, grep a few SNPs and double check with the bim file and strand file.

########################deal with TOPBOT strands####################

To understand the designation method TOPBOT illumina use

https://www.illumina.com/documents/products/technotes/technote_topbot.pdf

Illumina data files
For the Illumina chips the zip file contains three files, the .strand file, the .miss file and the .multiple file. The strand file contains six columns, SNP id, chromosome, position, %match to genome, strand and alleles. The SNP ids used are those from the annotation file and so are not necessarily the latest from dbSNP. The alleles listed are the Illumina TOP alleles, if you are in any doubt whether your data file can be used with these strand files a check of the non A/T G/C SNPs alleles vs the strand file should confirm this for you. If there are differences then it is likely your genotype file has been created using a different set of alleles, in this case if you can provide a list of the SNP ids and their alleles on the chip it is likely a strand file can be created for you.
The .miss file gives the ids of the SNPs that did not reach the required threshold for mapping to the genome, the position and strand of the best match are given. The .multiple file contains SNPs that had more than 1 high quality match (>90%) to the genome, in this instance the better match is taken for the .strand file. 

The strand files assume that your genotype calling algorithm has standardised the genotype allele calls to the Illumina TOP strand. 

Other Allele Sets

Illumina chips are not always called using the TOP strand and if this is the case for your data file then the most likely alternative is using the alleles derived from the Source Sequence in the annotation file, I have labelled the files, where they exist, as SourceStrand and will be adding them as time, and demand, allows.

A script developed by Neil Robertson for updating the chomosome, position and strand of binary ped files using these strand and position files can be downloaded here:

update_build.sh

Usage is 

	update_build.sh <bed-file-stem> <strand-file> <output-file-stem>

where:
 
<bed-file-stem> is the name of your binary ped set minus the .bed, .bim or .fam extension

<strand-file> is appropriate strand file for you chip and current strand orientation (TOP, SOURCE, etc)

<output-file-stem> is the name of the new output file to create again minus the .bed, .bim or .fam extension

:) Thanks to Sonia.

########################Convert the kgp SNP identifiers##################

#do this step in local PC.
#don't convert the kgp numbers to rs numbers before flip, because the strand file use the same identifiers as the original data.

scp tian.lin@delta.imb.uq.edu.au:/shares/compbio/PCTG/methylation/mQTL_project/1_COHORT/COHORT_genotype/fliped.bim  /?/

#use key_scripts/rsids.R to convert it, and output the converted bim, and a list of SNPs that were converted.

	R
	bimname <- "flip_Tian_cleaned_PD.bim"
	bim <- read.table(bimname, colClasses=c("numeric", "character", "numeric", "numeric", "character", "character"))
	bim$index <- 1:nrow(bim)

	##23-->X
	bim[(bim[,1]==23),1]="X"
	###if there are other chr than 1~23, remove them use plink,
	###Don't remove them only in bim file.

	# For each chromosome download the SNP locations and match to bim file
	a <- ddply(bim, .(V1), .progress="text", function(x)
	{
	  x <- mutate(x)
	  chr <- paste("ch", x$V1[1], sep="")
	  snps <- getSNPlocs(chr)
	  snps <- subset(snps, loc %in% x$V4, select=-c(alleles_as_ambig))
	  snps$RefSNP_id <- paste("rs", snps$RefSNP_id, sep="")
	  snps <- subset(snps, !duplicated(loc))
	  snps <- subset(snps, !duplicated(RefSNP_id))
  
	  x <- merge(x, snps, by.x="V4", by.y="loc", all.x=T)
	  x <- x[order(x$index), ]
	  index <- !is.na(x$RefSNP_id)
	  x$V2[index] <- x$RefSNP_id[index]
	  x <- subset(x, select=c(V1, V2, V3, V4, V5, V6))
	  return(x)
	})

	# If there are duplicated SNPs for any reason then label them as _2, _3 etc
	temp <- rle(a$V2)
	temp2 <- paste0(rep(temp$values, times = temp$lengths), "_", unlist(lapply(temp$lengths, seq_len)))
	temp2 <- gsub("_1", "", temp2)
	a$V2 <- temp2


	##CAUTION: the SNPs in the bim file are sorted in 1, 11, 12, .....19, 2, 21, 22, 3, 4, 5, 6, 7,..
	##reorder the SNPs

	a[a[,1]=="X",1]=23

	sorted.a=a[with(a, order(as.numeric(V1),V4)),]

	write.table(sorted.a, file="idsed.fliped.QLD_MND.bim", row=F, col=F, qu=F)

	scp converted.bim tian.lin@delta.imb.uq.edu.au:/shares/compbio/PCTG/methylation/COHORT_genotype/

##cp the fam and bed file with the same name as the converted bim file

########################separate into vcf format on chromosomes#######

##whole panel is available.

HRC.r1.GRCh37.autosomes.mac5.sites.tab 

39235158 lines in HRC panel

##softwares are available in modules in delta server. 
	module avail 					

	module load tabix/0.2.6		##for the function bgzip

	module load vcftools/0.1.15 	##for the function vcf-sort

	##files are separated based on chromosome.

	for i in $(seq 1 22); do ./plink --bfile idsed.fliped.QLD_MND  --chr $i --recode vcf --out QLD_MND-chr$i; done

	for i in $(seq 1 22); do vcf-sort QLD_MND-chr$i.vcf | bgzip -c > QLD_MND-chr$i.vcf.gz; done

	##remove the SNPs in hh file for X chromosome

	./plink --bfile idsed.fliped.QLD_MND  --exclude SNPs.in.hh.txt --chr 23 --recode vcf --out QLD_MND-chr23

	##The X chromosome is denoted as 23 in the genotype file. but X in the reference panel.

	##Remember to change 23 into X in the vcf file.

	sed 's/23/X/' QLD_MND-chr23.vcf > QLD_MND-chrX.vcf

	vcf-sort QLD_MND-chrX.vcf | bgzip -c > QLD_MND-chrX.vcf.gz

########################submit to imputation#########################

##apply for an account in the website 

https://imputationserver.sph.umich.edu/

#login and click Run

#pick Parameters:

Reference Panel: HRC r1.1 2016

select file

Phasing: Eagle v2.3 (phased output); ShapeIT for X chr

Population (for QC only): ?

Mode: Quality Control and Imputation

click yes

click yes

click Start Imputation


#CAUTION: choose ShapeIT v2.r790 (unphased) as the phasing tool, not the eagle v2.3. otherwise will fail
#Maximum of jobs can be uploaded is 3.
#Download the results with wget and save in COHORT_genotype/imputed/
	use screens to do several at the same time.
	tried with parallel and qsub, failed, always disconnect.


:) Thanks to Yang, Beben, Allan, Jian.


###############################check the imputation results#########################

##check the plot in qcreport.html

##Useful code:

/usr/bin/unzip 

zless -S chr22.dose.vcf.gz 

vcftools --gzvcf input.vcf.gz --plink --out output

##Release a group of files that share the same password at once.

	for zipfile in *.zip; do

	/usr/bin/unzip -P 'iaBx|62TibB{HX' $zipfile

	done



#################Use bcftools to check and fix the switched allels##############

This step is needed if use Sanger imputation server.

Find the plugin is at /opt/Modules/bcftools/1.4.1/libexec/bcftools/. The two files in need are fixref.c and fixref.so. 

qsub the following script. 

	#!/bin/sh
	#PBS -l walltime=48:00:00
	#PBS -l select=1:ncpus=1:mem=128gb

	cd /shares/compbio/PCTG/methylation/mQTL_project/4_PD/PD_genotype/separated_with_resorted/chrX/

	module load bcftools/1.4.1
	export BCFTOOLS_PLUGINS=/opt/Modules/bcftools/1.4.1/libexec/bcftools/

	ref='/shares/compbio/PCTG/methylation/mQTL_project/human_g1k_v37.fasta'
	input='/shares/compbio/PCTG/methylation/mQTL_project/4_PD/PD_genotype/separated_with_resorted/chrX/PD_nomismatch-chrX.vcf.gz'
	output='/shares/compbio/PCTG/methylation/mQTL_project/4_PD/PD_genotype/separated_with_resorted/chrX/PD_nomismatch-chrX_fixed.vcf.gz'

	bcftools +fixref $input -- -f $ref
	bcftools +fixref $input -Oz -o $output -- -d -f $ref -m flip
	bcftools +fixref $output -- -f $ref

:) Thanks to Irfahan.

########################sanger imputation service####################

##another option:

First step: registered a globus ID: lintian1986z@globusid.org

Second step: install globus connect on computer. 

Third step: go to sanger frontpage and submit a job, then get a run ID in the email. Click the link, so linked to the page to transfer data (https://www.globus.org/app/transfer)

Choose "my laptop at uq", which was named by myself. Find the file to upload and click on the arrow. 

Then, receive an email from globus that my file has been uploaded successfully.

Then, receive an email from Sanger that to conform the imputation. Click on the link and then on the "confirm" button.

To check the status, click on "status", and put in the Run ID.

When receive an email from Globus, click the link to download the result.

When globus send another email saying that the download is succeeded, go to the email from Sanger and click the link to remove your data from globus.

##Check the result:

:) Thanks to Yang.

###############Check the preparation with the check script.###################

cd /shares/compbio/PCTG/methylation/mQTL_project/4_PD/PD_genotype/

Requires the unzipped tab delimited HRC reference (currently v1.1 HRC.r1-1.GRCh37.wgs.mac5.sites.tab) from the Haplotype Reference Consortium Website here:

 http://www.haplotype-reference-consortium.org/site
 
Usage: perl HRC-1000G-check-bim-v4.2.pl -b <bim file> -f <Frequency file> -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

Never tried this yet.

Help page:
http://www.well.ox.ac.uk/~wrayner/tools/#Checking





#######################compare Sanger vs Michigan ###########

Michigan is easier in submit the job. Just need to apply an account and click on the webpage.

Sanger needs to install globus to submit the data.

No need to separate into 22 chromosomes when submit to Sanger.

Sanger communicate with around 6 emails per job.

Sanger results don't need to be extracted with a password.

Sanger output used as much rs IDs as possible. Only some "." and duplicate IDs need to be fixed.

Michigan output all used chr:position as SNP ID.

Both of their output are separated based on chromosomes.

Need to be merge chromosomes together if use GCTA --mlma-loco




###################### convert vcf.gz file to plink binary file ################

## CAUTION  

##Before run the following command,

##1. check the individual ID, if there's multiple _ in the FID_IID field, use --const-fid, if only the _ between FID_IID, use --double-id  --id-delim '_', so the ID in vcf file will be split into FID and IID.

##2. if imputed by Michigan, use  

	bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/R2\t%INFO/ER2\n' "chr"$i".dose.vcf.gz" > "imputed_chr"$i".info"
	
if imputed by Sanger, use
	
	bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/RefPanelAF\t%INFO/INFO\n' $vcf$i".dose.vcf.gz" > $outdir"BSGS_chr"$i".info"

Because their outputs have different names of the information in INFO column.


	#!/bin/sh
	#PBS -l walltime=48:00:00
	#PBS -l select=1:ncpus=1:mem=128gb
	cd /shares/compbio/PCTG/methylation/mQTL_project/Sanger_imputed/BSGSautosome.vcfs/

	module load bcftools/1.4.1
	outdir="plink_format/"
	for((i=1;i<=22;i++))
	do
    ./plink2 --vcf $i".vcf.gz" --double-id  --id-delim '_' --keep-allele-order --make-bed --out $outdir"BSGS_imputed_chr"$i > $outdir"BSGS_chr"$i".log" 
	bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/RefPanelAF\t%INFO/INFO\n' $i".vcf.gz" > $outdir"BSGS_chr"$i".info"
	done


################### Fix missing IDs in the bim file ########################

	R

	library(data.table)
	for (i in 1:22){
	bim=fread(paste("plink_format/BSGS_imputed_chr",i,".bim",sep=""))
	missid=which(bim$V2==".")
	bim$V2[missid]=paste("chr",bim$V1[missid],":",bim$V4[missid],sep="")
	write.table(bim,paste("fixed_ID/BSGS_imputed_chr",i,".bim",sep=""),col=F,row=F,sep="\t",quote=F)
	}

	rename.duplicate <- function (x, sep = "_dup", verbose = FALSE) 
	{
	x <- as.character(x)	
	duplix <- duplicated(x)	
	duplin <- x[duplix]	
	ix <- numeric(length = length(unique(duplin)))	
	names(ix) <- unique(duplin)	
	retval <- numeric(length = length(duplin))	
	for (i in 1:length(duplin)) {
	retval[i] <- ix[duplin[i]] <- ix[duplin[i]] + 1
	}
	
	retval <- retval + 1	
	x[duplix] <- paste(duplin, retval, sep = sep)	
	if (verbose) {
		message(sprintf("%i duplicated names", length(duplin)))
	}	
	return(list(new.x = x, duplicated.x = duplin))
	}

	for (i in 1:22){
	bim=fread(paste("fixed_ID/BSGS_imputed_chr",i,".bim",sep=""))
	dups <- unique(bim$V2[duplicated(bim$V2)])
	ndup <- length(dups)
	if(ndup == 0)
	{
	cat("No duplicate SNPs\n")	
	q()
	} else {
	cat(paste(ndup, "duplicate SNPs\n"))
	}
	bim$V2 <- rename.duplicate(bim$V2, sep="_dup")[[1]]
	write.table(bim,paste("fixed_ID/BSGS_imputed_chr",i,".bim",sep=""),col=F,row=F,sep="\t",quote=F)
	}


:) Thanks to Yang!


#################### merge imputed chromosomes ##################

vi the following file.

##CAUTION: no chr1 in this file.

##CAUTION: put the .bed filename first, then the .bim filename, then the .fam filename

##otherwise: Error: Line 1 of BSGS_chr2.bed has fewer tokens than expected.

allfiles.txt

	BSGS_chr2.bed BSGS_chr2.bim BSGS_chr2.fam 
	BSGS_chr3.bed BSGS_chr3.bim BSGS_chr3.fam 
	BSGS_chr4.bed BSGS_chr4.bim BSGS_chr4.fam 
	BSGS_chr5.bed BSGS_chr5.bim BSGS_chr5.fam 
	BSGS_chr6.bed BSGS_chr6.bim BSGS_chr6.fam 
	BSGS_chr7.bed BSGS_chr7.bim BSGS_chr7.fam 
	BSGS_chr8.bed BSGS_chr8.bim BSGS_chr8.fam 
	BSGS_chr9.bed BSGS_chr9.bim BSGS_chr9.fam 
	BSGS_chr10.bed BSGS_chr10.bim BSGS_chr10.fam 
	BSGS_chr11.bed BSGS_chr11.bim BSGS_chr11.fam 
	BSGS_chr12.bed BSGS_chr12.bim BSGS_chr12.fam 
	BSGS_chr13.bed BSGS_chr13.bim BSGS_chr13.fam 
	BSGS_chr14.bed BSGS_chr14.bim BSGS_chr14.fam 
	BSGS_chr15.bed BSGS_chr15.bim BSGS_chr15.fam 
	BSGS_chr16.bed BSGS_chr16.bim BSGS_chr16.fam 
	BSGS_chr17.bed BSGS_chr17.bim BSGS_chr17.fam 
	BSGS_chr18.bed BSGS_chr18.bim BSGS_chr18.fam 
	BSGS_chr19.bed BSGS_chr19.bim BSGS_chr19.fam 
	BSGS_chr20.bed BSGS_chr20.bim BSGS_chr20.fam 
	BSGS_chr21.bed BSGS_chr21.bim BSGS_chr21.fam 
	BSGS_chr22.bed BSGS_chr22.bim BSGS_chr22.fam 

combine BSGS_chr1 with the ones in list with following script:

	#!/bin/sh
	#PBS -l walltime=48:00:00
	#PBS -l select=1:ncpus=1:mem=128gb
	cd /home/tian.lin/mQTL_project/Imputed_data/1_BSGS/GWAS/plink_format
	./plink2 --bfile BSGS_chr1 --merge-list allfiles.txt --make-bed --out BSGS_imputed_autosomes


##################### Use GCTA to do GWAS ########################

#some of my data have family structures, so picked gcta --mlma-loco to do the GWAS.

#Generate the phenotype file first. It should have three columns: FID, IID, phenotype.

Introduction from the GCTA webpage:

"This option will implement an MLM based association analysis with the chromosome, on which the candidate SNP is located, excluded from calculating the GRM. We call it MLM leaving-one-chromosome-out (LOCO) analysis. The model is

y = a + bx + g- + e

where g- is the accumulated effect of all SNPs except those on the chromosome where the candidate SNP is located. The var(g-) will be re-estimated each time when a chromosome is excluded from calculating the GRM. The MLM-LOCO analysis is computationally less efficient but more powerful as compared with the MLM analysis including the candidate (--mlma).

The results will be saved in the *.loco.mlma file."

	#!/bin/sh
	#PBS -l walltime=48:00:00
	#PBS -l select=1:ncpus=16:mem=128gb
	cd /home/tian.lin/mQTL_project/Sanger_imputed/BSGSautosome.vcfs/GWAS
	./gcta64 --bfile BSGS_imputed_autosomes --mlma-loco  --maf 0.01 --thread-num 16 --pheno mean.beta.for.GWAS.txt  --out mean.beta.mlma

:) Thanks to Allan. 












