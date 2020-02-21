#!/software/R-3.3.0/bin/Rscript
OUT_DIR <- Sys.getenv('OUT_DIR')
get_info_score <- function(variants) {
  library(data.table)
  all_min <- fread(paste0(OUT_DIR, "/all_min.tsv"), select=c("VARIANT", "INFO"));
  resp <- merge(data.frame(VARIANT=variants), all_min, by="VARIANT", all.x=TRUE)
  resp <- resp[match(variants, resp$VARIANT),]
  return(resp$INFO)
}

user_friendly_file_name <- function(traits) {
  resp <- c()
  for (trait in traits) {
    resp <- c(resp,
              switch(trait,
                     mcv="mcv", hgb="hgb", hct="hct", mch="mch", baso="baso", wbc="wbc", eo="eo", mrv="mrv",
                     lymph="lymph", mpv="mpv", pdw="pdw", neut="neut", "lymph_p"="lymph_p", rbc="rbc",
                     ret="ret", eo_p="eo_p", mono="mono", plt="plt", rdw_cv="rdw", pct="pct", hlr="hlsr",
                     hlr_p="hlsr_p", irf="irf", mscv="mscv", neut_p="neut_p", ret_p="ret_p", mono_p="mono_P",
                     mchc="mchc", baso_p="baso_p", stop("error")))
  }
  return(resp)
}

user_friendly_trait_name <- function(traits) {
  resp <- c()
  for (trait in traits) {
    resp <- c(resp,
              switch(trait,
                     mcv="MCV", hgb="HGB", hct="HCT", mch="MCH", baso="BASO#", wbc="WBC#", eo="EO#", mrv="MRV",
                     lymph="LYMPH#", mpv="MPV", pdw="PDW", neut="NEUT#", "lymph_p"="LYMPH%", rbc="RBC#",
                     ret="RET#", eo_p="EO%", mono="MONO#", plt="PLT#", rdw_cv="RDW", pct="PCT", hlr="HLSR#",
                     hlr_p="HLSR%", irf="IRF", mscv="MSCV", neut_p="NEUT%", ret_p="RET%", mono_p="MONO%",
                     mchc="MCHC", baso_p="BASO%", stop("error")))
  }
  return(resp)
}

get_trait_cell_type <- function(traits, user_friendly=FALSE) {
  resp <- c()
  for (trait in traits) {
    if (trait %in% c("plt", "mpv", "pdw", "pct")) { resp <- c(resp, "platelet") }
    else if (trait %in% c("rbc", "mcv", "hct", "mch", "mchc", "hgb", "rdw_cv", "mscv")) { resp <- c(resp, "red_cell") }
    else if (trait %in% c("ret", "ret_p", "irf", "hlr", "hlr_p", "mrv")) { resp <- c(resp, "red_cell") }
    else if (trait %in% c("mono", "neut", "eo", "baso")) { resp <- c(resp, "myeloid_white_cell") }
    else if (trait %in% c("lymph")) { resp <- c(resp, "lymphocyte") }
    else if (trait %in% c("wbc", "mono_p", "neut_p", "eo_p", "baso_p", "lymph_p")) { resp <- c(resp, "white_cell") }
  }
  if (user_friendly==TRUE) {
    resp <- tools::toTitleCase(gsub("_", " ", resp))
  }
  return(resp)
}
h <- function(w) if( any( grepl("running command \\'grep \\-r \\'.* had status 1", w) ) ) invokeRestart( "muffleWarning" )
h1 <- function(w) if( any( grepl("invalid factor level, NA generated", w) ) ) invokeRestart( "muffleWarning" )
# split output of grep search, required some postprocessing especially if
# multiple SNPs exist at that basepair position
process_grep_output <- function(header, output, variant) {
  output <- strsplit(output, "\t")
  snp_output <- ""
  if (length(output) == 1) {
    snp_output <- output[[1]]
    names(snp_output) <- header
    snp_output <- as.list(snp_output)
    return(snp_output)
  } else if (length(output) > 1) {
    # loop through all the outputs and find the one which is relevant
    for (item in 1:length(output)) {
      this_snp_output <- output[[item]]
      names(this_snp_output) <- header
      this_snp_output <- as.list(this_snp_output)
      bp <- as.integer(gsub(".+:(.+)_.+_.+","\\1", variant))
      ref <- gsub(".+:.+_(.+)_.+","\\1", variant)
      alt <- gsub(".+:.+_.+_(.+)","\\1", variant)
      chr <- gsub("(.+):.+_.+_.+","\\1", variant)
      if (this_snp_output$CHR == chr && this_snp_output$BP == bp
          && this_snp_output$ALLELE1 == ref && this_snp_output$ALLELE0 == alt) {
        snp_output <- this_snp_output
        return(snp_output)
      }
    }
  } else {
    stop(paste("for variant", variant, "grep output is empty", output))
  }
}

# retrieve the BETA, SE, PVAL and LOGP for this trait
get_variant_trait_boltlmm <- function(variant, trait) {
   chr <- gsub("(.+):.+_.+_.+","\\1", variant)
   bp <- as.integer(gsub(".+:(.+)_.+_.+","\\1", variant))
   filename <- paste(Sys.getenv('BOLT_OUT_DIR'), "/", trait, "_gwas_normalised_imputed.out", sep="")

   # get headers first (this is inefficient but I want my code to adapt if BOLT-LMM changes order of cols)
   cmd <- paste("head -n1", filename)
   output <- system(cmd, intern=TRUE)
   header <- strsplit(output, "\t")[[1]]
   cmd <- paste("grep -r '", bp, "' ", filename, sep="")
   withCallingHandlers(output <- system(cmd, intern=TRUE), warning=h)
   split_output <- process_grep_output(header, output, variant)

   # check output to ensure it's what we expect
   ref <- gsub(".+:.+_(.+)_.+","\\1", variant)
   alt <- gsub(".+:.+_.+_(.+)","\\1", variant)
   if (split_output$CHR != chr || split_output$BP != bp
       || split_output$ALLELE1 != ref || split_output$ALLELE0 != alt)
   {
     stop(paste("Error in get_variant_trait_boltlmm(), with", variant, "for trait,", trait, "ref is", split_output$ALLELE0, "alt is", split_output$ALLELE1, "output length", length(split_output)))
   }
   
   # flip the beta, because BOLT-LMM was given the ref and alt in reverse
   split_output$BETA <- as.numeric(split_output$BETA) * -1
   
   # calculate MLOG10P
   split_output$MLOG10P=-pchisq((split_output$BETA/as.numeric(split_output$SE))^2,df=1,lower.tail=FALSE,log.p=TRUE)/log(10)
   # calculate MAF
   split_output$ALT_FREQ=1-as.numeric(split_output$A1FREQ)
   split_output$ALT_MINOR <- as.numeric(split_output$ALT_FREQ)<0.5
   split_output$MA_FREQ <- as.numeric(split_output$ALT_FREQ)
   split_output$MA_FREQ[!split_output$ALT_MINOR] <- 1 - split_output$MA_FREQ[!split_output$ALT_MINOR]
   out_list <- list(variant, as.numeric(split_output$BETA), as.numeric(split_output$SE), as.numeric(split_output$P_BOLT_LMM_INF), as.numeric(split_output$MLOG10P), as.numeric(split_output$MA_FREQ))
   names(out_list) <- c("VARIANT", "BETA", "SE", "P", "MLOG10P", "MAF") 
   # get a df of the variant, beta, se and p value for this variant
   return(as.data.frame(out_list, stringsAsFactors=F))
}

get_rsid_from_variantid <- function(variantids, mapping_file=NULL) {
    variantids = gsub("XY:", "X:", variantids)
  if (is.null(mapping_file)) {
    mapping_file <- fread(Sys.getenv('MAPPING_FILE_RSID'))
    if (sum(is.element(variantids, mapping_file$COORDID)) != length(variantids)) {
      mapping_file <- fread(mapping_file)
    }
  }
  return(mapping_file$dbSNP[match(variantids, mapping_file$COORDID)])
}

get_variantid_from_rsid <- function(rsids, mapping_file=NULL) {
  library(data.table)
  if (is.null(mapping_file)) {
    mapping_file <- fread(Sys.getenv('MAPPING_FILE_RSID'))
    if (sum(is.element(rsids, mapping_file$dbSNP_49)) != length(rsids)) {
      mapping_file <- fread(paste(Sys.getenv('SUPPORT_FILES'), "/id_mapping_table_rsids.tsv", sep=""), drop=c("dbSNP_47"))
    }
  }
  return(mapping_file$COORDID[match(rsids, mapping_file$dbSNP_49)])
}

