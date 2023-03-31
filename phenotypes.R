library(tidyverse)
library(data.table)
library(TCGAbiolinks)

ids <- fread('../sayaman/map_tcga_to_birdseed_id.txt', header = F) #map TCGA IDs to Birdseed (SNP chip) IDs - from Sayaman
ancestry <- fread('../sayaman/ucsf_ancestry.csv', sep = ',', header = T) #get ancestry assignments and PCs for each patient
types <- fread('../ASCAT/sample_average_ploidy_ASCAT.txt', header = T) %>% #get cancer types and GDC_aliquot IDs
	separate(cases, into = letters[1:7], sep = '-', remove = F) %>% 
	mutate(pt = paste(a, b, c, sep = '-'))
med <- fread('abs_asc.txt') 


##Get clinical information (age, sex)
allproj <- getGDCprojects()
projs <- allproj[startsWith(allproj$id, 'TCGA'),]$id
#clin <- lapply(projs, FUN = GDCquery_clinic, 'clinical') %>% bind_rows()
#write.table(clin, 'clin_query.txt', col.names = T, row.names = F, quote = F, sep= '\t')
clin <- fread('clin_query.txt', header = T)

##Integrate phenotypes into one df: take patient barcode and  integrate sex, age, birdseed ID, cancer type, ancestry, and ASCAT WGD calls (phenotype)
phen <- clin %>% 
	select(bcr_patient_barcode, gender, age_at_index, age_at_diagnosis, race) %>% 
	mutate(sex = recode(gender, male = 1, female = 2, .default = 0)) %>% #plink codes sex as 1 = male, 2 = female
	left_join(ancestry, by = c('bcr_patient_barcode' = 'Patient_ID')) %>%  
	filter(!is.na(Aliquot_ID)) %>% #why?
	left_join(ids, by = c('bcr_patient_barcode' = 'V3')) %>% 
	rename(birdseed = V1.y) %>% 
	left_join(types, by = c('bcr_patient_barcode' = 'pt')) %>%
	filter(!is.na(gender), !is.na(pam.ancestry.cluster), !is.na(age_at_index), !is.na(proj)) %>%
	left_join(med, by = c('bcr_patient_barcode' = 'Patient')) %>%
	mutate(phenotype = recode(comparison, "Neither WGD" = 1, "Both WGD" = 2, "ABSOLUTE only" = -9, "ASCAT only" = -9, .default = 0)) %>% #plink codes phenotype as 1 = controls, 2 = cases, -9 or 0 = missing
	mutate(phenotype = replace_na(phenotype, -9)) %>% 
	mutate(fam = '0') %>% #add constant FID (family ID) column, set to 0
	mutate(birdseed2 = paste(birdseed, birdseed, sep = '_')) #the .fam files seem to have the ID pasted together twice, for some reason - do this to match

write.table(phen, file = 'patient_info.txt', sep = '\t', col.names = T, row.names = F, quote = F)

phen %>%
	select(fam, birdseed2, sex, phenotype) %>% 
	write.table(file = 'phen.txt', sep = '\t', quote = F, col.names = F, row.names = F)


