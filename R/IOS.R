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
    id_list <- ids$id 
    message("Using default list of ", nrow(ids), " traits")
  }
  
  if(type[1] == "advanced"){
    ids <- subset(ao) %>%
      arrange(desc(sample_size)) %>%
      filter(!duplicated(trait), mr == 1) %>%
      filter(!grepl("ukb-a", id)) %>%
      filter(! id %in% c(id_exp, id_out))
    id_list <- ids$id 
    message("Using default list of ", nrow(ids), " traits")
  }
  return(id_list)
}

# make background dataset
#bg_dat <- make_background(snplist = exp_dat$SNP, id_bg = id_bg)
make_background <- function(snplist= exp_dat$SNP, id_bg = id_bg) {
  bdat <- TwoSampleMR::extract_outcome_data(snps = snplist, outcomes = id_bg)
  
  #caculate R2 reference- ieugwasr::gwasinfo(idlist)
  
  return(bdat)

}

# generate ios
#ios_dat <- ios(exp_dat=exp_dat, bg_dat=bg_dat)
ios <- function(exp_dat=exp_dat, bg_dat=bg_dat){
  require(dplyr)
  require(reshape2)
  bg_dat$vgu <- bg_dat$beta.outcome^2 * 2 * bg_dat$eaf.outcome * (1 - bg_dat$eaf.outcome)
  exp_dat$vgx <- exp_dat$beta.exposure^2 * 2 * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure)
  bg_dat <- merge(bg_dat, subset(exp_dat, select=c(SNP, vgx)), by="SNP")
  bg_dat$r2_ratio <- bg_dat$vgu / bg_dat$vgx
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
