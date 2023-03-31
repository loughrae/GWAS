# GWAS

- download TCGA and 1000 Genome Project VCFs
- impute variants using the Michigan Imputation Server
- measure ploidy phenotype: map copy number segment data to genes, calculate genome-wide ploidy, run MEDICC to infer Whole Genome Doubling events
- QC: quality-control SNPs and individuals; filter for relatedness, Hardy-Weinberg equilibrium, multi-allelism, linkage disequilibrium, and do PCA
- GWAS: Run GLM on quantitative phenotype and Firth logistic model on qualitative phenotype to identify significant SNPs
- MAGMA: Combine SNP signals to identify significant hits at the gene- and pathway levels
- Heritability analysis: calculate heritability of each phenotype using LDSC (linkage disequilibrium) and GCTA-GREML (relatedness matrix)
