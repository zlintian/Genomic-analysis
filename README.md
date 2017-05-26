# Genomic-analysis

Currently there are five cohorts I have finished imputation:
BSGS, collected from Allan.
QLD_MND, collected from Beben.
CNH_MND, collected from Beben.
PD, collected from Costanza.
LBCW1, collected from Qian.

The genotype data can be found in each directory under /shares/compbio/PCTG/methylation/mQTL_project/#_cohort in Delta cluster.

There's a readme file in each cohort. You can read it and decide whether you would like to start with my imputed data, or QC and impute by yourself.
Always check with me before you use the imputed data to make sure I have updated it to the latest I have.

My email: tian.lin@imb.uq.edu.au.


########################Genotype preparation workflow###############
The following steps should fit most SNP data, but there can be special situations come up.

0. define sample list and remove duplicates 
1. recode to ATGC if original data was coded with 1234
2. liftover to build 37 if in build 36 (this step can be omited because it can be done with the flip step, step6.
3. get the subset of data of selected individuals
4. pick the subset of chr1~23 if in original has 1~26. Only 1~23 can be imputed.
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

ssh ????@delta.imb.uq.edu.au

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

## now we have both map and ped file ready.

##run the following command using qsub is better. it takes quite a few minutes.

python liftMap.py -m LBCW1_recoded.map -p LBCW1_recoded.ped -o LBCW1_recoded_lifted

##the outputs:
SUCC:  map->bed succ
Reading liftover chains
Mapping coordinates
SUCC:  liftBed succ
SUCC:  bed->map succ
SUCC:  liftPed succ

LBCW1_recoded_lifted.ped
LBCW1_recoded_lifted.map
LBCW1_recoded_lifted.bed
LBCW1_recoded_lifted.bed.unlifted
LBCW1_recoded.map.bed

LBCW1_recoded.map.bed is the BED format file which is reformatted from the plink bed file. Don't need it.

LBCW1_recoded_lifted.bed file is in BED format, we don't need it either. change the name into LBCW1_recoded_lifted.bed.format, so won't mass up as plink bed file.

LBCW1_recoded_lifted.map and LBCW1_recoded_lifted.ped are what we need.

? SNPs in the unlifted file. All of them are shown as #Deleted in new?


:) Thanks to Yang.

########################look into the QLD_MND genotype data carefully##
? SNPs in the bim file.
? individuals in the genotype data est. (belong to 666 families)
? individuals in the methylation data set.
? individuals in the information sheet

#are genotype data overlap with methylation data set?
? of them are overlapped.
? extra in methylation data set
? extra in genotyped data set.

awk '{print $1}' COHORT_subset.bim | uniq -c 
## how many families? 
## are there multiple samples in any family?

########################Get the subset of genotype##################

Determine what samples to keep based on the methylation data.
#1. Use the ID and the matched family ID (which is the same of ID) in final.txt to form a table "genotype.IDs.txt".
#2. check the chromosomes included in the genotype data
#3. use plink to get the subset of genotype data.

./plink --file LBCW1_recoded_lifted --keep genotype.IDs.txt --chr 1-23  --make-bed --out LBCW1_subset

#how many samples in raw?
#how many samples to keep?
#how many SNPs input?
#how many SNPs kept?

########################QC##########################################

./plink --bfile cohort --mind 0.05 --geno 0.05 --hwe 0.000001 --maf 0.01  --make-bed --out Tian_cleand_cohort

#make a record of variables and samples removed and kept.
? people input
? SNPs input

? people removed
? SNPs removed as missing
? SNPs removed due to hwe
? SNPs removed due to low MAF

? people kept
? SNPs kept?

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

Usage is update_build.sh <bed-file-stem> <strand-file> <output-file-stem>
where:

 
<bed-file-stem> is the name of your binary ped set minus the .bed, .bim or .fam extension

<strand-file> is appropriate strand file for you chip and current strand orientation (TOP, SOURCE, etc)

<output-file-stem> is the name of the new output file to create again minus the .bed, .bim or .fam extension

:) Thanks to Sonia.

########################Convert the kgp SNP identifiers##############
#do this step in local PC.
#don't convert the kgp numbers to rs numbers before flip, because the strand file use the same identifiers as the original data.

scp tian.lin@delta.imb.uq.edu.au:/shares/compbio/PCTG/methylation/mQTL_project/1_COHORT/COHORT_genotype/fliped.bim  /?/

#use key_scripts/rsids.R to convert it, and output the converted bim, and a list of SNPs that were converted.
#how many converted?

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

module load tabix/0.2.6			##for the function bgzip

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

Mode: Quality Control & Imputation

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




####################### Compare between Sanger and Michigan imputation server.

Sanger results don't need password.
Sanger results don't need to be extracted.
No need to separate into 22 chromosomes when submit to Sanger.

bother need to merge together if use GCTA 


######################

convert vcf.gz file to plink binary file


####################

merge
















