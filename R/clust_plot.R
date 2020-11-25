#' Generate a plot for the simulation result
#'
#' Generate a plot for the result of permutation test
#'
#' @param simulation Result from the simulation (permutation test) 
#' @param Q_radial Q statistics obtained from original radial MR
#' @param Q_mrios Q statistics obtained from \code{mr_cluster_heterogeneity}
#'
#' @export
#' @return Plot

permutation_plot <- function(simulation = NULL, Q_radial = NULL, Q_mrios = NULL)
{
  sim <- simulation %>% 
    purrr::map_df(dplyr::as_tibble)
  
  p <- plot(density(sim$value), main="Distribution of Q statistics")
  
  abline(v = Q_radial, col = "red", lwd = 1)
  abline(v = Q_mrios, col = "blue", lwd = 1)

  return(p)
  
  }



#' Forest plot
#' 
#' Generate forest plot of the MR estimates from each cluster
#' 
#' @param result Results obtained using \code{mr_cluster}
#'
#' @export
#' @return Plot


mr_cluster_forest <- function(result = NULL)
{
  test <-  result$estimate %>% 
    bind_rows() %>%
    dplyr::mutate(ID = row_number())
  
  temp <- meta::metagen(b,
                        se,
                        data=test,
                        studlab=paste(ID),
                        comb.fixed = TRUE,
                        comb.random = FALSE,
                        prediction=TRUE,
                        sm="SMD")
  
  p <- meta::forest(temp,
                    layout = "JAMA",
                    text.predict = "95% PI",
                    col.predict = "black",
                    colgap.forest.left = unit(15,"mm"))
  return(p)
}


