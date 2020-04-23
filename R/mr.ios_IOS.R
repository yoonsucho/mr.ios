#' Calculate index of suspicion
#'
#' Calculate IOS across based on R squared multiple traits and summarising estimates using mean, sd, iqr, median, 95% value, maximum value
#' @param exp_dat Instruments for the exposure, obtained using \code{extract_instruments}
#' @param bg_dat Effects for the instruments on a set of variables, used to calculate index of suspicion
#' @examples ios(exp=exp_dat, bg=bg_dat)
#'
#' @export
#' @return Data.frame


# calculate ios
ios <- function(exp=exp_dat, bg=bg_dat){
  require(dplyr)
  require(reshape2)
  ios <- dplyr::group_by(bg, SNP) %>%
    dplyr::summarise(
      ios1_mean = sum(rsq.outcome, na.rm=TRUE), #sum or mean
      ios1_sd = sd(rsq.outcome, na.rm=TRUE),
      ios1_iqr = quantile(rsq.outcome, 0.75, na.rm=TRUE) - quantile(rsq.outcome, 0.25, na.rm=TRUE),
      ios1_median = median(rsq.outcome, na.rm=TRUE),
      ios1_95 = quantile(rsq.outcome, 0.95, na.rm=TRUE),
      ios1_max = max(rsq.outcome, na.rm=TRUE),
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
