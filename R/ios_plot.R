#' Plot of IOS MR result
#'
#' A function for generating scatter plot of IOS-Radial IVW MR result, obtained from \code{mr.ios()}. The size of each circle represents the relative weight that SNP carried in the IVW MR analysis. The x-axis indicates Wj, and the y axis indicates betaWj.
#' 
#' @param iosmr Result of Radial IVW MR analysis obtained using \code{mr.ios()}. 
#' @param ios IOS score obtained using \code{ios()}
#' @param ios_type Types of IOS estimators. Defalt value is \code{"ios1_mean"}. 
#'
#' @export
#' @return Plot

#plot - IOS
ios_plot <- function(iosmr = ios_mr, ios = ios_dat, ios_type="ios1_mean"){
  dat1 <- iosmr$data
  temp <- merge(dat1, ios, by = "SNP")
  
  p <- temp %>%
        dplyr::arrange(Qj) %>%
        ggplot2::ggplot(ggplot2::aes(y=BetaWj, x=Wj)) +
        ggplot2::geom_smooth(method='lm', formula = y ~ x, se=FALSE, size=0.5) +
        ggplot2::geom_point(ggplot2::aes(size = ios[[ios_type]]), colour = "deepskyblue", alpha = 1) +
        #geom_line(aes(colour=as.factor(method))) +
        #geom_hline(yintercept=648.0386, linetype="dashed", color = "grey") +
        ggplot2::labs(x="Wj", y="BetaWj", size="ios score")  
  }

