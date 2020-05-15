#plot - IOS
ios_plot <- function(iosmr = ios_mr, ios = ios_dat, ios_type="ios1_mean"){
  dat1 <- iosmr$data
  temp <- merge(dat1, ios, by = "SNP")
  
  p <- temp %>%
        arrange(Qj) %>%
        ggplot(aes(y=BetaWj, x=Wj)) +
        geom_smooth(method='lm', formula = y ~ x, se=FALSE, size=0.5) +
        geom_point(aes(size = ios[[ios_type]]), colour = "deepskyblue", alpha = 1) +
        #geom_line(aes(colour=as.factor(method))) +
        #geom_hline(yintercept=648.0386, linetype="dashed", color = "grey") +
        labs(x="Wj", y="BetaWj", size="ios score")  
  
}

