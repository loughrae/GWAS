

#### Downloads ####


# wget https://api.gdc.cancer.gov/data/9df36037-3f9a-47aa-bcab-1c5554bb0443 -O map_tcga_to_birdseed_id.txt
# wget https://api.gdc.cancer.gov/data/fdfa536a-c3c8-405d-99d9-bc9375b5084c -O ucsf_ancestry.csv
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz 

#### Shared processing of Patients ####
Rscript compare_calls.R  #compare copy number calls from ASCAT and ABSOLUTE and get consensus WGD call
Rscript phenotypes.R #put together patient clinical and WGD information
grep "EUR Sub-Cluster 1" patient_info.txt | cut -f64 | awk '{print "0 \t" $1}' > eur1.txt
cut -f1,2,3 phen.txt > sex.txt
cut -f1,2,4 phen.txt > phens.txt


#### Process variants from GWAS catalog ####
Rscript gwascat.R #filter and format GWAS catalog SNPs
../liftOver -positions catalog_snps.txt ../hg38ToHg19.over.chain.gz catalog_snps_lifted.txt unmapped #lift GWAS catalog SNPs from hg38 to hg19


for chrom in {1..22}
do
#plink --vcf chr${chrom}.dose.vcf.gz --make-bed --out chr${chrom}_store --const-fid 0
plink2 --bfile chr${chrom}_store --set-all-var-ids @:# --make-bed --out chr${chrom}_recode  #remove ref/alt from SNP names
Rscript --vanilla dedup_qual.R ${chrom} #remove multi-allelic SNPs and filter for imputation quality > 0.5
plink --bfile chr${chrom}_recode --extract chr${chrom}_dedup_qual.snp --update-sex sex.txt --pheno phens.txt --prune --make-bed --out chr${chrom}_qc #add sex and phenotype info; remove samples without phenotype
rm chr${chrom}_recode* #save storage space
plink --bfile chr${chrom}_qc --keep eur1.txt --hardy midp -out chr${chrom}_har_eur1 #compute Hardy-Weinberg on European Sub-Cluster 1 (from Sayaman's clustering)
done

#### Merge chromosomes, remove HWE deviations, LD-prune and do PCA ####
R -e "library(tidyverse) ; data.frame(paste0('chr', 1:22, '_qc.bed'), paste0('chr', 1:22, '_qc.bim'), paste0('chr', 1:22, '_qc.fam')) %>% write.table(file = 'tomerge.txt', col.names = F, row.names = F, quote = F, sep = '\t')"
grep "UNAFF" chr*_har_eur1.hwe | awk '$10 < 1e-6' | awk '{print $3}' | awk -F ":" '{print $0 "\t" "chr" $1 ":" $2 "-" $2}' | grep -v -f catalog_snps_lifted.txt | awk '{print $1}' > hwe_out.txt #find HWE deviants (except those listed as cancer GWAS hits)
plink --merge-list tomerge.txt --exclude hwe_out.txt --maf 0.01 --make-bed --out merged #merge chromosomes and exclude HWE deviants and MAF < 1% 


#### PCA ####

#merge pruned 1000G chromosomes into 1
find ../1000G -name "*.bim" | grep -e "Pruned" > ForMerge.list ;
sed -i 's/.bim//g' ForMerge.list ;
plink --merge-list ForMerge.list --out Merge ;
cut -f2 Merge.bim > mergedvars.txt #make list of variants from 1000 GP

#get intersection of 1000G and Sayaman variants
plink2 --bfile merged --set-all-var-ids @:#:\$r:\$a --extract mergedvars.txt --maf 0.1 --exclude range high-LD-regions-hg19-GRCh37.txt --make-bed --out filtered_sayaman
plink --bfile filtered_sayaman --indep-pairwise 50 5 0.2 --out filtered_sayaman_ld
plink --bfile filtered_sayaman --extract filtered_sayaman_ld.prune.in --make-bed --out pruno #extract the Sayaman filtered variants  from Sayaman
plink --bfile Merge --extract filtered_sayaman_ld.prune.in --make-bed --out intersec #extract the Sayaman filtered variants from 1000G
plink --bfile intersec --bmerge pruno --out combined #merge Sayaman and 1000G
plink2 --bfile combined --pca approx --out combined_pca #perform PCA (approx for speed)

#### Solo TCGA PCA ####
plink2 --bfile merged --maf 0.1 --exclude range high-LD-regions-hg19-GRCh37.txt --make-bed --out say_alone
plink --bfile say_alone --indep-pairwise 50 5 0.2 --out say_alone_ld
plink2 --bfile say_alone --extract say_alone_ld.prune.in --pca approx --out sayaman_alone

#### Solo 1000G PCA ####
plink2 --bfile Merge --pca approx --out 1000G_alone #note 1000G was already filtered for MAF >= 10%  


