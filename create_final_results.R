#!/software/R-3.2.2/bin/Rscript
OUT_DIR <- Sys.getenv("OUT_DIR")
PHENO_NAMES <- Sys.getenv("PHENO_NAMES")
CORES <- Sys.getenv('CORES')
traits <- as.vector(read.table(PHENO_NAMES, header=F, stringsAsFactors=F)$V1)
#traits = read.table(paste0(OUT_DIR, "/pheno_names_bcx.tsv"), header=F, stringsAsFactors=F)$V1
# load user friendly file name function
# and user friendly trait name function
source("general_functions.R")

condsig.df <- NULL
for (trait in traits) {
  trait_condsig <- read.table(paste(OUT_DIR, "/condout/igv/", user_friendly_file_name(trait), "_condsig.gwas", sep=""), header=T, stringsAsFactors=F, sep="\t")
  #trait_condsig$trait <- rep(user_friendly_trait_name(trait), nrow(trait_condsig))
  trait_condsig$trait <- rep(trait, nrow(trait_condsig))
  trait_condsig$class <- rep(get_trait_cell_type(trait), nrow(trait_condsig))
  condsig.df <- rbind(condsig.df, trait_condsig)
}


# load vep results
vep.table <- read.table(paste(OUT_DIR, "/condout/vep/tabular_all_condsig.txt", sep=""), header=T, stringsAsFactors=F)
names(vep.table)[names(vep.table) == "ID"] <- "VARIANT"
vep.table <- vep.table[, names(vep.table) %in% c("VARIANT", "Consequence", "IMPACT", "SYMBOL", "Gene")]
print(nrow(condsig.df))
condsig.df <- merge(condsig.df, vep.table, by="VARIANT", all.x=T)
print(nrow(condsig.df))

# create Minor Allele column
condsig.df$minorallele <- condsig.df$ALT
condsig.df$minorallele[condsig.df$ALT_MINOR==FALSE] <- condsig.df$REF[condsig.df$ALT_MINOR==FALSE]

# cond effect direction
condsig.df$conddirection <- c("-", "+")[as.numeric(condsig.df$JOINT_EFFECT>0)+1]

# add infoscore information
condsig.df$INFO <- get_info_score(condsig.df$VARIANT)

# create the hg19 coordid for each snp and use this to merge with the ld clumping table
condsig.df$hg19 <- paste(condsig.df$CHR, ":", condsig.df$BP, sep="")
#hg19_clump_table <- read.table(paste(OUT_DIR, "/condout/ldclump_dosage/hg19_clump_table.tsv", sep=""), header=T, stringsAsFactors=F)
hg19_clump_table = read.table(paste0(OUT_DIR, "/astle_compare/hg19_clump_table_novel.tsv"), header=T)
condsig.df <- merge(condsig.df, hg19_clump_table, by="hg19", all.x=TRUE)
condsig.df$clump_id = gsub("chr:", "", condsig.df$clump_id)

# sort the dataframe along the genome
condsig.df$CHR_num = condsig.df$CHR
condsig.df$CHR_num[condsig.df$CHR_num %in% c("X", "XY")] = 23
condsig.df$sort=as.numeric(condsig.df$CHR_num)*10^9+as.numeric(condsig.df$BP)
condsig.df <- condsig.df[order(condsig.df$sort),]
condsig.df <- condsig.df[, !(names(condsig.df) == "sort")]
print(unique(condsig.df$CHR))
print(nrow(condsig.df))


# reorder the columns and exluce the un-needed ones
condsig.df <- condsig.df[c("trait", "novel", "clump_id", "class", "VARIANT", "ID", "INFO", "CHR", "BP", "REF", "ALT", "minorallele", "ALT_FREQ", "MA_FREQ", "EFFECT", "SE", "MLOG10P", "R2", "JOINT_EFFECT", "JOINT_SE", "JOINT_MLOG10P", "VEP_CONSEQUENCE", "VEP_IMPACT", "VEP_GENE_SYMBOL", "VEP_GENE_ENSEMBL_ID")]
#condsig.df <- condsig.df[c("trait", "class", "clump_id", "novel", "VARIANT", "ID", "INFO", "CHR", "BP", "REFALT", "minorallele", "ALT_FREQ", "MA_FREQ", "EFFECT", "SE", "MLOG10P", "R2", "JOINT_EFFECT", "JOINT_SE", "JOINT_MLOG10P", "VEP_CONSEQUENCE", "VEP_IMPACT", "VEP_GENE_SYMBOL", "VEP_GENE_ENSEMBL_ID")]
write.csv(condsig.df, "data/BCX_final_output_raw.csv", row.names=F)
condsig.df$class_reduced <- NULL
condsig.df$class <- get_trait_cell_type(condsig.df$trait, user_friendly=TRUE)
condsig.df$trait <- user_friendly_trait_name(condsig.df$trait)

names(condsig.df)[names(condsig.df) == "trait"] = "Associated Blood Index"
names(condsig.df)[names(condsig.df) == "class"] = "Associated Blood Index Class"
names(condsig.df)[names(condsig.df) == "clump_id"] = "Locus ID"
names(condsig.df)[names(condsig.df) == "novel"] = "Novel vs Astle 2016"
names(condsig.df)[names(condsig.df) == "VARIANT"] = "Unique Variant ID"
names(condsig.df)[names(condsig.df) == "ID"] = "rsID (where available)"
names(condsig.df)[names(condsig.df) == "INFO"] = "INFO SCORE"
names(condsig.df)[names(condsig.df) == "CHR"] = "Chr (GRCh37)"
names(condsig.df)[names(condsig.df) == "BP"] = "BP (GRCh37)"
names(condsig.df)[names(condsig.df) == "REF"] = "REF (GRC37)"
names(condsig.df)[names(condsig.df) == "ALT"] = "ALT (GRC37)"
names(condsig.df)[names(condsig.df) == "minorallele"] = "Minor Allele"
names(condsig.df)[names(condsig.df) == "ALT_FREQ"] = "Alternative Allele Frequency" 
names(condsig.df)[names(condsig.df) == "MA_FREQ"] = "Minor Allele Frequency"

names(condsig.df)[names(condsig.df) == "DIRECTION"] = "(UNIVAR) Estimate of Effect Direction (REF=Baseline, ALT=Effect)"
names(condsig.df)[names(condsig.df) == "EFFECT"] = "(UNIVAR) Estimate of Additive Allelic Effect (REF=Baseline, ALT=Effect)"
names(condsig.df)[names(condsig.df) == "SE"] = "(UNIVAR) Standard Error of Estimator"
names(condsig.df)[names(condsig.df) == "MLOG10P"] = "(UNIVAR) -log10 P"
names(condsig.df)[names(condsig.df) == "R2"] = "(UNIVAR) Unadjusted R2"

names(condsig.df)[names(condsig.df) == "conddirection"] = "(MULTI) Estimate of Effect Direction (REF=Baseline, ALT=Effect)"
names(condsig.df)[names(condsig.df) == "JOINT_EFFECT"] = "(MULTI) Estimate of Additive Allelic Effect (REF=Baseline, ALT=Effect)"
names(condsig.df)[names(condsig.df) == "JOINT_SE"] = "(MULTI) Standard Error of Estimate"
names(condsig.df)[names(condsig.df) == "JOINT_MLOG10P"] = "(MULTI) -log10 P"

names(condsig.df)[names(condsig.df) == "VEP_CONSEQUENCE"] = "Most Serious VEP Consequence of Variant"
names(condsig.df)[names(condsig.df) == "VEP_IMPACT"] = "VEP IMPACT of Most Serious Consequence"
names(condsig.df)[names(condsig.df) == "VEP_GENE_SYMBOL"] = "Gene Symbol(s) for Most Serious Consequence"
names(condsig.df)[names(condsig.df) == "VEP_GENE_ENSEMBL_ID"] = "Ensembl Gene ID(s) for Most Serious Consequence"
names(condsig.df)[names(condsig.df) == "PHENOSCANNER_DISEASE"] = "Phenoscanner: Associated Disease"
names(condsig.df)[names(condsig.df) == "CLINVAR_PHENOTYPE_LIST"] = "ClinVar: Phenotype List"
names(condsig.df)[names(condsig.df) == "CLINVAR_GENE_SYMBOL"] = "ClinVar: Gene Symbol"
names(condsig.df)[names(condsig.df) == "MEMBRANE_GO_CC"] = "GO Cellular Compartment: 'membrane'"

write.table(condsig.df, "data/BCX_final_output.tsv", sep="\t", quote=F, row.names=F, na="")
