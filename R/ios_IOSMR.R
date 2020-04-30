#' Calculate index of suspicion
#'
#' Calculate IOS across based on R squared multiple traits and summarising estimates using mean, sd, iqr, median, 95% value, maximum value
#' @param exp Instruments for the exposure, obtained using \code{TwoSampleMR::extract_instruments}.
#' @param bg Effects for the instruments on a set of variables, used to calculate index of suspicion, obtained using \code{make_background()}.
#' 
#'
#' @export
#' @return Data.frame


# calculate ios
ios <- function(exp=exp_dat, bg=bg_dat){
  requireNamespace("reshape2", quietly = TRUE)
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
  temp <- stats::reshape(ios, timevar="variable", idvar="SNP", direction="wide")
  names(temp)[-1] <- as.character(unique(ios$variable))
  
  return(temp)
}




#' Perform MR accounting for pleiotropy using IOS
#'
#' Performs Radial inverse variance weighted (IVW) model using IOS to detect the most pleiotropic SNPs based on IOS score. The SNPs with higher score for IOS are downweighted by multiplying IOS by inverse variance weight implemented in Radial MR package.
#'  
#' @param dat Harmonised dataset of the exposure and the outcome, obtained using \code{TwoSampleMR::harmonise_data()}
#' @param ios IOS score obtained using \code{ios()}
#' @param alpha Statistial significance threshold for identifying outliers. Default value is \code{0.05}.
#' @param weights Inverse variance weights used to calculate IVW estimate and Cochran's Q statistic, considering IOS in the model. Detailed information of the definition of each weight is provided here: \code{https://github.com/WSpiller/RadialMR}. Select first order (\code{1}), second order (\code{2}) or modified second order weights (\code{3}).
#' @param tol Tolerance threshold for performing the iterative IVW approach. Default value is \code{0.0001}.
#' 
#'
#' @export
#' @return List of estimates

# perform mr
#dat <- harmonise_data(exp_dat, out_dat)
#mr_res <- mr.ios(dat=dat, ios=ios_dat)

mr.ios <-function(dat=dat, ios = ios, alpha = 0.05, weights, tol = 0.0001){
  
  weights <- weights + 3
  dat_rmr <- RadialMR::format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP, ios$ios1_mean, ios$SNP)
  rares <- RadialMR::ivw_radial(dat_rmr, alpha, weights, tol, external_weight = TRUE)
  
  return(rares)
  
}
