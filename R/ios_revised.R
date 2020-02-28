ios <- function(exp_dat=exp_dat, bg_dat=bg_dat)
{
  require(dplyr)
  require(reshape2)
  
  #if(bg_dat$Category != "log odds"){
    # Calculate SD from SE and the number of SNPs if unit is not SD
    #### if (bg_dat$unit != "SD"){} ? 
    bg_dat$sd.outcome <-bg_dat$se.outcome * sqrt(length("number of SNPs tested in the original dataset")) #needs to be confirmed
    exp_dat$sd.exposure <-exp_dat$se.exposure * sqrt(length("number of SNPs tested in the original dataset"))
    
    #Calculate variance explained by the candidiate traits 
    bg_dat$vgu1 <- TwoSampleMR::get_r_from_pn(bg_dat$pval.outcome, bg_dat$samplesize.outcome)
    
    bg_dat$vgu2 <- bg_dat$beta.outcome^2 * 2 * bg_dat$eaf.outcome * (1 - bg_dat$eaf.outcome) / bg_dat$sd.outcome
    
    #Calculate variance explained by the original exposure
    exp_dat$vgx <- exp_dat$beta.exposure^2 * 2 * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure) / exp_dat$sd.exposure
    
}

TwoSampleMR::get_r_from_lor

  
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
