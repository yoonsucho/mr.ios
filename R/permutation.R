#' Permutation test. 
#' 
#' @export
#' @return Simulation result in data frame

perm_sim <- function(dat = dat, bg = bg_dat, num_test = 1000){
sim <- list()
invisible(capture.output(sim <- lapply(1:num_test, function(x) {
  tryCatch(
    {
  dat <- dat
  bg_dat_temp <-transform(bg_dat, rsq.outcome = sample(rsq.outcome))
  clust <- hclust_instruments(bg_dat_temp, value_column = "rsq.outcome", kmax=min(50, length(unique(bg_dat_temp$SNP))))
  dat_clust <- merge(dat, clust, by = "SNP")
  hclust <- list()
  hclust$dat <- dat_clust
  hclust$estimate <- mr_cluster_estimate(dat = hclust$dat)
  l <- list()
  l$dat <- clust
  l$Q <- mr_cluster_heterogeneity(cluster = hclust)
  return(l)
  },
  error = function(error){
  return(NULL)
  }
  )
}
))
)
}
