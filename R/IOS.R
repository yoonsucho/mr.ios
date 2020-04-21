#' MR.IOS:: Calculate index of suspicion
#' 
#' A package for sensitivity analysis in two-samle Mendelian randomization setting when horizontal pleiotropy exists. 
#' It uses [IEU GWAS database](https://gwas.mrcieu.ac.uk/) to obtain data automatically, 
#' and is implemented to RadialMR to run the analysis. 
#' 
#' If a SNP influences multiple other traits then it could be 'suspicious', and more likely to be pleiotropic. This function implements two basic approaches to estimate IOS
#'
#' - ios1: A summary of the SNP r2 with the other traits (r2_gu).
#' - ios2: A summary of the ratio of r2_gu / r2_gx, where r2_gx is the variance explained by the SNP on the exposure. Estimates the index of suspicion, whereupon SNPs which have a larger effect on a set of traits given their effect on the exposure are deemed more suspicious.
#' 
#'
#' @name mr.ios-package
#' @aliases TwoSampleMR twosamplemr
#' @docType package
NULL
