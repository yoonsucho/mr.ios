#' Calculate index of suspicion
#'
#' If a SNP influences multiple other traits then it could be 'suspicious', and more likely to be pleiotropic. This function implements two basic approaches to estimate IOS
#'
#' - ios1: A summary of the SNP r2 with the other traits (r2_gu)
#' - ios2: A summary of the ratio of r2_gu / r2_gx, where r2_gx is the variance explained by the SNP on the exposure. Estimates the index of suspicion, whereupon SNPs which have a larger effect on a set of traits given their effect on the exposure are deemed more suspicious
#'
#'
#' \code{background_ids} Returns list of IDs of available GWAS summary statistics
#' @param id_exp ID for the exposure in IEU GWAS database. 
#' @param id_out ID for the outome in IEU GWAS database.
#' 
#' 
#' \code{ios} Summarising across multiple traits can be dune using mean, sd, iqr, median, 95% value, maximum value
#' @param exp_dat Instruments for the exposure, obtained using \code{extract_instruments}
#' @param bg_dat Effects for the instruments on a set of variables, used to calculate index of suspicion
#'
#' @export
#' @return Data.frame

# choose exposure
id_exp <- "ieu-a-2"
exp_dat <- TwoSampleMR::extract_instruments(id_exp)

# choose outcome
id_out <- "ieu-a-7"
out_dat <- TwoSampleMR::extract_outcome_data(exp_dat$SNP, id_out)

# choose background dataset list
# if info about category / n samples .... are missing -> exclude 

background_ids <- function(id_exp=id_exp, id_out=id_out, type=c("default", "advanced")[1]){
  ao <- suppressMessages(TwoSampleMR::available_outcomes())
  if(type[1] == "default") {
    ids <- subset(ao) %>%
      arrange(desc(sample_size)) %>%
      filter(!duplicated(trait), mr == 1) %>%
      filter(grepl("ukb-b", id)) %>%
      filter(! id %in% c(id_exp, id_out))
    
    ids <- subset(ids, category == "NA" & unit == "NA")
    id_list <- ids$id 
    
    message("Using default list of ", nrow(ids), " traits")
  }
  
  if(type[1] == "advanced"){
    ids <- subset(ao) %>%
      arrange(desc(sample_size)) %>%
      filter(!duplicated(trait), mr == 1) %>%
      filter(!grepl("ukb-a", id)) %>%
      filter(! id %in% c(id_exp, id_out))
    
    ids <- subset(ids, category == "NA" & unit == "NA")
    id_list <- ids$id 
    
    message("Using default list of ", nrow(ids), " traits")
  }

  return(id_list)

  }

# make background dataset
#bg_dat <- make_background(snplist = exp_dat$SNP, id_bg = id_bg)
make_background <- function(exp = exp_dat, id_bg = id_bg) {
  snplist <- exp$SNP
  
  bdat <- TwoSampleMR::extract_outcome_data(snps = snplist, outcomes = id_bg)
  
  #call the IEU GWAS database to look up the information about "sd", "category" and "unit"
  ao <- TwoSampleMR::available_outcomes()
  info <- ieugwasr::gwasinfo(unique(c(bdat$id.outcome, exp_dat$id.exposure[1])))
  
  #Add columns of "category", "unit" and "SD" to the background data frame
  ##category
  bdat$category[bdat$id.outcome %in% info$id] <- info$category[info$id %in% bdat$id.outcome]
  ##unit
  bdat$unit[bdat$id.outcome %in% info$id] <- info$unit[info$id %in% bdat$id.outcome]
  ##SD
  bdat$sd.outcome <- NA
  
  #for continuous outcome
  #if(! (bdat$category %in% c('Disease','Binary'))){
  #    for(i in 1:nrow(bdat)){
  #      bdat$sd.outcome[bdat$id.outcome[i] %in% ao$id[i]] <- ao$sd[i]
  #    }
  #}
  
  bdat <-
    bdat %>% 
    mutate(sd.outcome = ifelse(category %in% c("Risk factor", "Continuous", "Categorical Ordered", "Metabolites") & 
                                unit == "SD", 
                              1, 
                             ifelse(category %in% c("Disease", "Binary") &
                                      unit == "SD",
                                    1,
                                    ifelse(category %in% c("Risk factor", "Continuous") &
                                             !unit == "SD",
                                           ao$sd[ao$id %in% bdat$id.outcome],
                                           ifelse(category %in% c("Disease", "Binary") &
                                                    !unit == "SD",
                                                  ao$sd[ao$id %in% bdat$id.outcome], NA
                                           )))))
  ##ncase
  bdat$ncase[bdat$id.outcome %in% info$id] <- info$ncase[info$id %in% bdat$id.outcome]
  ##ncontrol
  bdat$ncontrol[bdat$id.outcome %in% info$id] <- info$ncontrol[info$id %in% bdat$id.outcome]
  
  
  #Add columns of "category", "unit" and "SD" to the exposure data frame
  ##category
  exp$category[exp$id.exposure%in% info$id] <- info$category[info$id %in% exp$id.exposure]
  ##unit
  exp$unit[exp$id.exposure %in% info$id] <- info$unit[info$id %in% exp$id.exposure]
  ##SD
  exp$sd.exposure <- NA
  exp <-
    exp %>% 
    mutate(sd.exposure = ifelse(category %in% c("Risk factor", "Continuous", "Categorical Ordered", "Metabolites") & 
                                 unit == "SD", 
                               1, 
                               ifelse(category %in% c("Disease", "Binary") &
                                        unit == "SD",
                                      1,
                                      ifelse(category %in% c("Risk factor", "Continuous") &
                                               !unit == "SD",
                                             ao$sd[ao$id %in% exp$id.exposure],
                                             ifelse(category %in% c("Disease", "Binary") &
                                                      !unit == "SD",
                                                    ao$sd[ao$id %in% exp$id.exposure], NA
                                             )))))
  ##ncase
  exp$ncase[exp$id.exposure %in% info$id] <- info$ncase[info$id %in% exp$id.exposure]
  ##ncontrol
  exp$ncontrol[exp$id.exposure %in% info$id] <- info$ncontrol[info$id %in% exp$id.exposure]
  
    
    
  #Calculate variance explained by the candidiate traits 
  bdat <-
    bdat %>% 
    mutate(vgu1 = ifelse(category %in% c("Risk factor", "Continuous", "Categorical Ordered", "Metabolites"),
                         TwoSampleMR::get_r_from_pn(bdat$pval.outcome, bdat$samplesize.outcome), 
                         ifelse(category %in% c("Disease", "Binary"),
                                TwoSampleMR::get_r_from_lor(bdat$beta.outcome, bdat$eaf.outcome, bdat$ncase, bdat$ncontrol, 0.1, model = "logit", correction = F), 
                                NA))) %>%
    mutate(vgu2 = ifelse(category %in% c("Risk factor", "Continuous", "Categorical Ordered", "Metabolites"),
                         bdat$beta.outcome^2 * 2 * bdat$eaf.outcome * (1 - bdat$eaf.outcome) / bdat$sd.outcome, 
                         ifelse(category %in% c("Disease", "Binary"),
                                TwoSampleMR::get_r_from_lor(bdat$beta.outcome, bdat$eaf.outcome, bdat$ncase, bdat$ncontrol, 0.1, model = "logit", correction = F), 
                                NA)))  
    
  #Calculate variance explained by the original exposure
  exp <-
    exp %>% 
    mutate(vgx1 = ifelse(category %in% c("Risk factor", "Continuous", "Categorical Ordered", "Metabolites"),
                         TwoSampleMR::get_r_from_pn(exp$pval.exposure, exp$samplesize.exposure), 
                         ifelse(category %in% c("Disease", "Binary") &
                                  is.na(exp$ncase) &
                                  is.na(exp$ncontrol),
                                TwoSampleMR::get_r_from_lor(exp$beta.exposure, exp$eaf.exposure, exp$ncase, exp$ncontrol, 0.1, model = "logit", correction = F), 
                                NA))) %>%
    mutate(vgx2 = ifelse(category %in% c("Risk factor", "Continuous", "Categorical Ordered", "Metabolites"),
                         exp$beta.exposure^2 * 2 * exp$eaf.exposure * (1 - exp$eaf.exposure) / exp$sd.exposure, 
                         ifelse(category %in% c("Disease", "Binary")&
                                  is.na(exp$ncase) &
                                  is.na(exp$ncontrol),
                                TwoSampleMR::get_r_from_lor(exp$beta.exposure, exp$eaf.exposure, exp$ncase, exp$ncontrol, 0.1, model = "logit", correction = F), 
                                NA)))  
  
  bg_dat <- merge(bg_dat, subset(exp_dat, select=c(SNP, vgx1, vgx2)), by="SNP")
  bg_dat$r2_ratio <- bg_dat$vgu / bg_dat$vgx
  
  return(bdat)

}




# generate ios
#ios_dat <- ios(exp_dat=exp_dat, bg_dat=bg_dat)
ios <- function(exp_dat=exp_dat, bg_dat=bg_dat){
  require(dplyr)
  require(reshape2)
  ios <- dplyr::group_by(bg_dat, SNP) %>%
    dplyr::summarise(
      ios1_mean = sum(vgu, na.rm=TRUE), #sum or mean
      ios1_sd = sd(vgu, na.rm=TRUE),
      ios1_iqr = quantile(vgu, 0.75, na.rm=TRUE) - quantile(vgu, 0.25, na.rm=TRUE),
      ios1_median = median(vgu, na.rm=TRUE),
      ios1_95 = quantile(vgu, 0.95, na.rm=TRUE),
      ios1_max = max(vgu, na.rm=TRUE),
      ios2_mean = sum(r2_ratio, na.rm=TRUE),
      ios2_sd = sd(r2_ratio, na.rm=TRUE),
      ios2_iqr = quantile(r2_ratio, 0.75, na.rm=TRUE) - quantile(r2_ratio, 0.25, na.rm=TRUE),
      ios2_median = median(r2_ratio, na.rm=TRUE),
      ios2_95 = quantile(r2_ratio, 0.95, na.rm=TRUE),
      ios2_max = max(r2_ratio, na.rm=TRUE)
    ) %>% melt

  # Reshape IOS
  temp <- reshape(ios, timevar="variable", idvar="SNP", direction="wide")
  names(temp)[-1] <- as.character(unique(ios$variable))
  
  return(temp)
}



# perform mr
#dat <- harmonise_data(exp_dat, out_dat)
#mr_res <- mr.ios(dat=dat, ios=ios_dat)

mr.ios <-function(ios = ios_dat, dat=dat, r_input,alpha,weights,tol){
 
  radat <- format_radial(ios = ios_dat, dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
  rares <- RadialMR::ivw_radial()
  
  return(rares)
  
}
