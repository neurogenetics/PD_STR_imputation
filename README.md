# PD_STR_imputation
Imputing STR calls from PD SNP based data using the Gymrek STR reference 1kg reference panel
# **IPDGC STR PIPELINE**

Date: March 2020
Authors: Kimberley Billingsley

## General description and purpose

This is a general description of the LNG STR imputation and GWAS pipeline. This  pipeline can be roughly divided into five steps

 1. STR Imputation
 2. STR Filtering & QC
 3. GWAS QC 
 4. GWAS
 5. META-analysis


## 1) STR Imputation 
 
 This pipeline takes directions from [here](http://gymreklab.com/2018/03/05/snpstr_imputation.html)
 -  Copy over target GWAS data plink files and remove monomorphic SNPs
 -  Download STR reference panel for each chromosome (1-22) from  [here](http://gymreklab.com/2018/03/05/snpstr_imputation.html)
-   Harmonize your target GWAS data with the reference panel
-   Perform STR imputation with Beagle
 
 *Copy over target GWAS data plink files and remove monomorphic SNPs*
  
    for f in `cat IPDGC.datasets.txt` 
    	do
    		for i in {1..22} ; do plink2 --vcf pre_impute${f}.${i}.vcf.gz --maf 0.0000001  --	double-id --recode vcf bgz  --out ${f}.${i};done
	    done
*Harmonize your target GWAS data with the reference panel*

    for f in `cat IPDGC.datasets.txt` 
    	do
		    for i in {1..22} ;do java -jar conform-gt.24May16.cee.jar gt=${f}.${i}.vcf.gz ref=1kg.snp.str.chr${i}.vcf.gz chrom=${i}  match=POS out=${f}.${i}.consistent; done
    	done
*Perform STR imputation with Beagle*

    for f in `cat IPDGC.datasets.txt` 
    	do
    		for i in {1..22} ;do java -Xmx200g -jar beagle.12Jul19.0df.jar gt=${f}.${i}.consistent.vcf.gz ref=1kg.snp.str.chr${f}.${i}.vcf.gz nthreads=20 window=10 overlap=1.0 out=${f}.${i}.IMPUTED.vcf.gz; done
    	done
    	
## 2) STR Filtering & QC 
 
 This pipeline takes directions from  [here](https://github.com/bibb/GWAS-QC-and-STR-imputation-analysis-pilenines/blob/master/STR_imputation_pipeline.md)

- Split SNPs and STRs 
- Manually change the VCF headers
- Perform multiallelic splitting on STRs with the VT software
- Create Unique IDs for each variant 
- Filter STRs by DR2 > 0.3

*Split SNPs and STRs*

    module load VCFTOOLS
    for f in `cat IPDGC.datasets.txt` 
	 	do
			for i in {1..22} do ;vcftools --gzvcf ${f}.${i}.IMPUTED.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip -c > onlySNPs.${f}.${i}.IMPUTED.vcf.gz |vcftools --gzvcf ${f}.${i}.IMPUTED.vcf.gz --keep-only-indels  --recode --recode-INFO-all --stdout | gzip -c > onlySTRs.${f}.${i}.IMPUTED.vcf.gz ;done
		done
 *VCF headers on the STRs requires manually changing to reflect this:*

OLD:
```
##INFO=<ID=DR2,Number=1
```

NEW:
```
##INFO=<ID=DR2,Number=A
```
*Manually change the VCF headers*

    for f in `cat IPDGC.datasets.txt` 
    	do	
	    	for i in {1..22} ;zcat onlySTRs.${f}.${i}.IMPUTED.vcf.gz | sed 's/##INFO=<ID=DR2,Number=1/##INFO=<ID=DR2,Number=A/g' | gzip -c > DR2A.onlySTRs.${f}.${i}.IMPUTED.vcf.gz ;done
    	done
	
*Perform multiallelic splitting on STRs with the VT software*

    module load VT
    for f in `cat IPDGC.datasets.txt`
    	do	
	    	for i in {1..22} do; vt decompose -s DR2A.onlySTRs.${f}.${i}.IMPUTED.vcf.gz -o DR2A.onlySTRs.${f}.${i}.split.IMPUTED.vcf.gz; done
    	done
*Create Unique IDs for each variant* 

    for f in `cat IPDGC.datasets.txt`
    	do	
    	for i in {1..22} 
    		do 
    			zcat DR2A.onlySTRs.${f}.split.${i}.vcf.gz | grep -F "#" > DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.vcf
    	    	zcat DR2A.onlySTRs.${f}.split.${i}.vcf.gz | grep -F -v "#" | awk 'BEGIN {OFS="\t"} {$3=$3"_"NR} {print}' >> DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.vcf 
    	    done
    	done
    	
    	
*Filter STRs by DR2 > 0.3*

  

     module load VCFTOOLS
        for f in `cat IPDGC.datasets.txt` 
        	do
        	for i in {1..22} 
        		do 
        		grep -v "#" DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.vcf  | cut -f3,8 | sed 's/DR2=/\t/g' | sed 's/;/\t/g' | awk '$2 > 0.3{print $1}' > chr${i}.${f}_newsplit_Dr2_03.temp 
        		vcftools --vcf DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.vcf  --snps chr${i}.${f}_newsplit_Dr2_03.temp --recode --recode-INFO-all --stdout | gzip -c > DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.gz
        		done
        	done\
        	
        	
        	
## 3) GWAS QC 

This pipeline takes directions from [here](https://github.com/neurogenetics/GWAS-pipeline)

*Make sure that covariate file has PCs per dataset

- Calculate PCs per cohort 
 - Reconstruct new covar file 
 - Create peopleIncludeFile
 - Create rangeFile
 -  Re-format VCF for RV-test 
 
   
        module load plink
        for f in `cat IPDGC.datasets.txt`
            do
        	    plink --vcf pre_impute${f}.${i}.vcf.gz --double-id --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --exclude exclusion_regions_hg19.txt --make-bed --out ${f}.clean_2  --memory 119500 --threads 19
        	    plink --bfile ${f}.clean_2 --indep-pairwise 1000 10 0.02 --autosome --out ${f}.pruned_data
        	    plink --bfile ${f}.clean_2 --extract ${f}.pruned_data.prune.in --make-bed --out ${f}.clean_3
        	    plink --bfile ${f}.clean_3 --pca --out ${f}.PCA
            done
*Combine the different datasets* 

    cat MCGILL.PCA.eigenvec MF.PCA.eigenvec NEUROX.PCA.eigenvec NIA.PCA.eigenvec OSLO.PCA.eigenvec PDBP.PCA.eigenvec PPMI.PCA.eigenvec PROBAND.PCA.eigenvec SHULMAN.PCA.eigenvec SPAIN3.PCA.eigenvec TUBI.PCA.eigenvec UK_GWAS.PCA.eigenvec VANCE.PCA.eigenvec COURAGEUK.PCA.eigenvec DUTCH.PCA.eigenvec FINLAND.PCA.eigenvec GERMANY.PCA.eigenvec HBS.PCA.eigenvec > All.PCA.all.IPDGC.cohorts.txt

*Reformat  to match vcf IDs* 

    awk '$1=$1"_"$2' All.PCA.all.IPDGC.cohorts.txt > Double.ID.All.PCA.all.IPDGC.cohorts.txt

*Manually add Header to Double.ID.All.PCA.all.ad.cohorts.txt with vim* 

*Format to match vcf* 

    awk '$1=$1"_"$2' IPDGC_all_samples_covariates.txt > Double.ID.IPDGC_all_samples_covariates.txt


*Merge together the PCs per dataset with covariate file in R*

    R
    old <-read.table(file="Double.ID.phenotype_file_plink_for_filtered_genotypes.txt", fill=T header=T)
    ID <-read.table(file="Double.ID.All.PCA.all.ad.cohorts.txt", header=T)
    merged< -merge(old, new, by="FID")
    write.table(merged, file="New.covariate.IPDGC.pc.per.cohort.txt", col.names=T, row.names=F, quote=F, sep="\t")
*Format new covariate file*

    awk '{print $1 ,$1 ,0 ,0 ,$3 ,$4 ,$5 ,$20 ,$21 ,$22 ,$23 ,$24 ,$25 ,$26 ,$27 ,$28 ,$29 ,$30}' New.covariate.IPDGC.pc.per.cohort.txt >  Final.covariate.AD.pc.per.cohort.txt

*Manually add Header to Final.covariate.IPDGC.pc.per.cohort.txt with vim*  

------

*Create peopleIncludeFile
 (the samples were already QCd per https://github.com/neurogenetics/GWAS-pipeline therefore we only need to generate a list of all the individuals in the vcf)*

   
    module load bcftools  
    for f in `cat IPDGC.datasets.txt`
        	do		
    	bcftools query -l  DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.gz  > ${f}.keep.txt
    	done

*Create Regionsfile, no output from imputation apart from vcf, manually removed low imputation quality so now filtering for MAF*

  
    module load plink  
    for f in `cat IPDGC.datasets.txt`
        		do
    			for i in {1..22}
    					do		
    					plink --vcf  DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.gz --vcf-half-call missing --maf 0.001 --make-bed --out ${f}.${i}.maf.001
    					awk '{print $1":"$4"-"$4}'  ${f}.${i}.maf.001.bim > keep.variant.${f}.${i}.txt
    					done
    		done

*Reformat the vcf* 


    module load bcftools 
    for f in `cat IPDGC.datasets.txt' 
    	do
    		for i in {1..22}
    		do
   
    		gunzip DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.gz 
    		bgzip -c DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf > r.vs.${f}.chr${i}.vcf.gz
    		tabix -p vcf r.vs.${f}.chr${i}.vcf.gz
    		done
    	done
    	
## 4) GWAS  

This pipeline takes directions from [here](https://github.com/neurogenetics/GWAS-pipeline)
RVtest needs to be run separately for the datasets that have AAO and for the ones that do not. 

AAO:

    module load rvtests/2.1.0
    
    for f in `cat All.for.age.rvtest.meta5.txt` 
    	do
    		for i in {1..22}
    		do
    		rvtest --noweb --hide-covar --rangeFile keep.variants.${f}.chr${i}.txt \
    		--out RV.${f}.chr${i} --single wald \
    		--inVcf r.vs.${f}.chr${i}.vcf.gz --dosage DS --pheno  New.covar.final.with.pc.per.dataset.IPDGC.txt \
    		--pheno-name pheno --covar New.covar.final.with.pc.per.dataset.IPDGC.txt \
    		--covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5 \
    		--peopleIncludeFile keep.people.${f}.chr${i}.txt
    	done
    
    done

NO AAO:

    module load rvtests/2.1.0
    
        for f in `cat No.age.age.rvtest.meta5.txt` 
        	do
        		for i in {1..22}
        		do
        		rvtest --noweb --hide-covar --rangeFile keep.variants.${f}.chr${i}.txt \
        		--out RV.${f}.chr${i} --single wald \
        		--inVcf r.vs.${f}.chr${i}.vcf.gz --dosage DS --pheno  New.covar.final.with.pc.per.dataset.IPDGC.txt \
        		--pheno-name pheno --covar New.covar.final.with.pc.per.dataset.IPDGC.txt \
        		--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5 \
        		--peopleIncludeFile keep.people.${f}.chr${i}.txt
        	done
        
        done

## 5) META-analysis 

 - Make maf001rsq03minimums_chr*.info file 
 - Merge all GWAS files per chromosome and merge all .info files
 - Reformat in R for METAL 
 - Run METAL for meta-analysis 

**Make maf001rsq03minimums_chr.info file*

    for f in `cat IPDGC.datasets.txt` 
    	do
    		for i in {1..22}
    		do
        		zcat DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.gz| grep -v "##" | awk '{print $1"_"$2"_"$5,$8}'  > DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.1
        		awk -F';IMP' '{print $1}'  DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.1 | awk '{print $1, $2}' | sed 's/DR2=//g' | sed 's/;AF=/\t/g' >  DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.2
        		sed 's/\\/&&/g;x;s/.*/1c\\/;G;q' Header.info.txt | sed -i -f - DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.2  
        		cp DR2A.onlySTRs.${f}.chr${i}.split.uniqueID.DR2_03.vcf.2 ${f}.chr${i}.info
    		done
    	done
*Merge all GWAS files per chromosome and merge all .info files*

.INFO:

    for f in `cat IPDGC.datasets.txt`;do cat  ${f}.chr1.info ${f}.chr2.info ${f}.chr3.info ${f}.chr4.info ${f}.chr5.info ${f}.chr6.info ${f}.chr7.info ${f}.chr8.info ${f}.chr9.info ${f}.chr10.info ${f}.chr11.info ${f}.chr12.info ${f}.chr13.info ${f}.chr14.info ${f}.chr15.info ${f}.chr16.info ${f}.chr17.info ${f}.chr18.info ${f}.chr19.info ${f}.chr20.info ${f}.chr21.info ${f}.chr22.info  | grep -v 'SNP' >  ${f}.all.chr.txt; done 

GWAS:

    for f in `cat IPDGC.datasets`;do cat RV.${f}.chr1.SingleWald.assoc RV.${f}.chr2.SingleWald.assoc  RV.${f}.chr3.SingleWald.assoc RV.${f}.chr4.SingleWald.assoc RV.${f}.chr5.SingleWald.assoc RV.${f}.chr6.SingleWald.assoc RV.${f}.chr7.SingleWald.assoc RV.${f}.chr8.SingleWald.assoc RV.${f}.chr9.SingleWald.assoc RV.${f}.chr10.SingleWald.assoc RV.${f}.chr11.SingleWald.assoc RV.${f}.chr12.SingleWald.assoc RV.${f}.chr13.SingleWald.assoc RV.${f}.chr14.SingleWald.assoc RV.${f}.chr15.SingleWald.assoc RV.${f}.chr16.SingleWald.assoc RV.${f}.chr17.SingleWald.assoc RV.${f}.chr18.SingleWald.assoc RV.${f}.chr19.SingleWald.assoc RV.${f}.chr20.SingleWald.assoc RV.${f}.chr21.SingleWald.assoc RV.${f}.chr22.SingleWald.assoc| grep -v 'SNP' > ${f}.allChrs_FILE.assoc; done

  *Reformat in R for METAL* 
  
*Pay attention when naming .infos colnames as it is different compared to the github LNG page, right way is SNP		Rsq		ALT_Frq

*strip header for both files before using this script 

    for f in `cat IPDGC.datasets.txt`;do tail -n +2 "maf001rsq03minimums.all.chr.$(f).info" > head.maf001rsq03minimums.all.chr.$(f).info

   Run this in R for each dataset

    infos <- read.table(paste("maf001rsq03minimums.all.chr.$(f).info"))
    names(infos) <-c("Test", "Rsq","ALT_Frq")
    infos$ALT_Frq <- as.numeric(as.character(infos$ALT_Frq))
    assoc <- read.table(file="$(f).allChrs_FILE.assoc")
    colnames(assoc) <- c("CHROM","POS","REF","ALT","N_INFORMATIVE","Old_marker","Beta","SE","Pvalue")
    assoc$Beta <- as.numeric(as.character(assoc$Beta))
    assoc$Test <- paste(assoc$CHROM, assoc$POS,assoc$ALT,sep="_")
    data <- merge(infos, assoc, by= "Test")
    dat <- subset(data, Beta < 5 & Beta > -5 & !is.na(data$Pvalue))
    dat$chr <- paste("chr",dat$CHROM, sep = "")
    dat$minorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$ALT), as.character(dat$REF))
    dat$majorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$REF), as.character(dat$ALT))
    dat$beta <- ifelse(dat$ALT_Frq <= 0.5, dat$Beta, dat$Beta*-1)
    dat$se <- dat$SE
    dat$maf <- ifelse(dat$ALT_Frq <= 0.5, dat$ALT_Frq, 1 - dat$ALT_Frq)
    dat$P <- dat$Pvalue
    dat0 <- dat[,c("SNP","minorAllele","majorAllele","beta","se","maf","P")]
    write.table(dat0, file=paste("$(f).toMeta.FILE.tab"), quote = F, sep = "\t", row.names = F)



*Run METAL for meta-analysis* 

Create a metal.txt file that looks like this:

```
#../generic-metal/metal metalAll.txt
#THIS SCRIPT EXECUTES AN ANALYSIS OF EIGHT STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Inputfile1.txt THROUGH Inputfile8.txt
SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
LABEL TotalSampleSize as N # If input files have a column for the sample size labeled as 'N'
# LOAD THE FIRST SEVEN INPUT FILES

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ   maf
EFFECT beta
STDERR se
PVALUE P
WEIGHT N 
PROCESS toMeta.GWAS1.tab

# === DESCRIBE AND PROCESS THE SECOND INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ   maf
EFFECT beta
STDERR se
PVALUE P
WEIGHT N
PROCESS toMeta.GWAS2.tab

```


    module load metal/2018-08-28

```
metal metal.txt
```


### References:

METAL:  [https://www.ncbi.nlm.nih.gov/pubmed/20616382](https://www.ncbi.nlm.nih.gov/pubmed/20616382)




