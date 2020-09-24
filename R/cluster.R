
#' Convert bg_dat to wide format
#'
#' Imputes missing values using missMDA::imputePCA function
#'
#' @param bg_dat Output from make_background
#' @param value_column Which column to use to fill the values of the matrix
#'
#' @export
#' @return Matrix
bg_to_wide <- function(bg_dat, value_column = "rsq.outcome")
{
	message("Restructuring")
	bg_dat$val <- bg_dat[[value_column]]
	wide <- bg_dat %>%
		dplyr::select(SNP, id.outcome, val) %>%
		tidyr::spread(key = id.outcome, value=val)
	rownames(wide) <- wide$SNP
	wide <- wide[,-1] %>% as.matrix()
	message("Imputing missing values")
	wide <- missMDA::imputePCA(wide)$completeObs
	message("Scaling data")
	wide <- as.matrix(scale(wide))
	return(wide)
}


#' Cluster instruments using k-means
#'
#'
#' @param bg_dat <what param does>
#' @param value_column <what param does>
#' @param kmax <what param does>
#' @param nstart <what param does>
#' @param iter.max <what param does>
#'
#' @export
#' @return
kmeans_instruments <- function(bg_dat, value_column, kmax=min(15, length(unique(bg_dat$SNP))), nstart=50, iter.max=15)
{
  
  #kmax=min(15, length(unique(bg_dat$SNP)))
	wide <- bg_to_wide(bg_dat, value_column)

	message("Performing k-means clustering for up to ", kmax, " clusters")
	l <- lapply(1:kmax, function(k)
	{
		message(k)
		kmeans(wide, k, nstart=nstart, iter.max=iter.max)
	})
	
	mse <- sapply(l, function(x) x$tot.withinss)
	cluster <- sapply(l, function(x) x$cluster)
	diff_mse <- diff(mse) * -1
	clust <- which.max(diff_mse) + 2
	dat <- dplyr::tibble(SNP = rownames(cluster), cluster = cluster[,clust])
	return(dat)
}


#' Cluster instruments using hclust
#'
#' @param bg_dat <what param does>
#' @param value_column <what param does>
#' @param kmax <what param does>
#'
#' @export
#' @return
hclust_instruments <- function(bg_dat, value_column, kmax=min(50, length(unique(bg_dat$SNP))))
{
	get_mse <- function(d, cuts)
	{
		dm <- as.matrix(d)
		mse <- sapply(1:ncol(cuts), function(j)
		{
			message(j)
			sapply(unique(cuts[,j]), function(i)
			{
				sn <- rownames(cuts)[cuts[,j] == i]
				sum(dm[rownames(dm) %in% sn, colnames(dm) %in% sn])
			}) %>% sum
		})
		return(mse)
	}
	wide <- bg_to_wide(bg_dat, value_column)
	pc <- prcomp(wide)
	d <- dist(pc$x)
	hc <- hclust(d)
	if (length(unique(bg_dat$SNP)) < 50)
	  {
	  cuts <- cutree(hc, k=1:kmax)
	  }
	if (length(unique(bg_dat$SNP)) >= 50)
	  {
	  cuts <- cutree(hc, k=1:50)
	  }
	mse <- get_mse(d, cuts)
	diff_mse <- diff(mse) * -1
	clust <- which.max(diff_mse) + 2
	dat <- dplyr::tibble(SNP = rownames(cuts), cluster = cuts[,clust])
	return(dat)
}

#' Cluster instruments using the pvclust method
#'
#' Hierarchical clustering model
#'
#' @param bg_dat <what param does>
#' @param value_column <what param does>
#' @param alpha <what param does>
#' @param method.hclust <what param does>
#' @param method.dist <what param does>
#' @param nboot <what param does>
#' 
#'
#' @export
#' @return

pvclust_instruments <- function(bg_dat, value_column, alpha=0.95, method.hclust="ward.D", method.dist="euclidean", nboot=200)
{
  kmax=min(15, length(unique(bg_dat$SNP)))
  wide <- bg_to_wide(bg_dat, value_column)
	pc <- prcomp(wide)
	fit <- pvclust::pvclust(t(pc$x), method.hclust=method.hclust, method.dist=method.dist, nboot=nboot)
	plot(fit)
	pvclust::pvrect(fit, alpha=alpha)
	#cluster <- pvclust::pvpick(fit, alpha = alpha)
	get_mse <- function(d, cuts)
	{
	  dm <- as.matrix(d)
	  mse <- sapply(1:ncol(cuts), function(j)
	  {
	    message(j)
	    sapply(unique(cuts[,j]), function(i)
	    {
	      sn <- rownames(cuts)[cuts[,j] == i]
	      sum(dm[rownames(dm) %in% sn, colnames(dm) %in% sn])
	    }) %>% sum
	  })
	  return(mse)
	}
	if (length(unique(bg_dat$SNP)) < 50)
	{
	  cuts <- cutree(fit$hclust, k=1:kmax)
	}
	if (length(unique(bg_dat$SNP)) >= 50)
	{
	  cuts <- cutree(fit$hclust, k=1:50)
	}
	mse <- get_mse(dist(pc$x), cuts)
	diff_mse <- diff(mse) * -1
	clust <- which.max(diff_mse) + 2
	dat <- dplyr::tibble(SNP = rownames(cuts), cluster = cuts[,clust])
	return(dat)
}


#' Cluster using partitioning around medoids
#'
#'
#' @param bg_dat <what param does>
#' @param value_column <what param does>
#' @param kmax=50 <what param does>
#' @param criterion="asw" <what param does>
#'
#' @export
#' @return
pam_instruments <- function(bg_dat, value_column, kmax=min(50, length(unique(bg_dat$SNP))), criterion="asw")
{
	wide <- bg_to_wide(bg_dat, value_column)
	o <- fpc::pamk(wide, krange=1:(kmax-1), criterion=criterion)$pamobject
	dat <- dplyr::tibble(SNP = names(o$clustering), cluster=o$clustering)
	return(dat)
}


#' MR in each cluster group
#' 
#' 
#' @param bg Effects for the instruments on a set of variables, used to calculate index of suspicion, obtained using \code{make_background()}.
#' @param dat Harmonised dataset of the exposure and the outcome, obtained using \code{TwoSampleMR::harmonise_data()}.
#' @param method Choose the method to generate cluster. 
#' \itemize{
#'   \item kmean See ?\code{kmeans_instruments}
#'   \item hclust See ?\code{hclust_instruments}
#'   \item pv See ?\code{pvclust_instruments}
#'   \item pam See ?\code{pam_instruments}
#' }
#' 
#' @export
#' @return 

mr_cluster <- function(bg = bg_dat, dat = dat, method = c("kmean", "hclust", "pv", "pam")){
  output <- list()
  #generate cluster using selected method
  if(method == "kmean"){
    clust <- kmeans_instruments(bg, value_column = "rsq.outcome", kmax=min(15, length(unique(bg$SNP))), nstart=50, iter.max=15)
  }
  
  if(method == "hclust"){
    clust <- hclust_instruments(bg, value_column = "rsq.outcome", kmax=min(50, length(unique(bg$SNP))))
  }
  
  if(method == "pv"){
    clust <- pvclust_instruments(bg, value_column = "rsq.outcome", alpha=0.95, method.hclust="ward.D", method.dist="euclidean", nboot=200)
  }
  
  if(method == "pam"){
    clust <- pam_instruments(bg, value_column = "rsq.outcome", kmax=min(50, length(unique(bg$SNP))), criterion="asw")
  }
  
  # Add cluster info to the harmonised data
  dat_clust <- merge(dat, clust, by = "SNP")
  
  output$dat <- dat_clust
  
  
  # Perform mr within cluster
  mr_cluster_estimate <- function(dat = dat){
    res_clust <- list()
    for (i in 1:max(dat[ , ncol(dat)]))
    {
      temp <- subset(dat, dat[ , ncol(dat)] == i) 
      if(nrow(temp) > 3)
      {
        Ratios <- temp$beta.outcome / temp$beta.exposure
        W <- ((temp$beta.exposure^2) / (temp$se.outcome^2))
        Wj <- sqrt(W)
        BetaWj <- Ratios * Wj
        IVW.Model <- lm(BetaWj~-1+Wj)
        EstimatesIVW <- summary(lm(IVW.Model))
        IVW.Slope <- EstimatesIVW$coefficients[1]
        IVW.SE <- EstimatesIVW$coefficients[2]
        IVW_CI<-confint(IVW.Model)
        DF <- length(temp$SNP)-1
        Qj <- W*(Ratios-IVW.Slope)^2
        Total_Q <- sum(Qj)
        Total_Q_chi <- pchisq(Total_Q,length(temp$beta.exposure)-1,lower.tail = FALSE)
        W<- ((temp$se.outcome^2+(IVW.Slope^2*temp$se.exposure^2))/temp$beta.exposure^2)^-1
        Wj<-sqrt(W)
        BetaWj<-Ratios*Wj
        IVW.Model<-lm(BetaWj~-1+Wj)
        EstimatesIVW<-summary(lm(BetaWj~-1+Wj))
        IVW.Slope<-EstimatesIVW$coefficients[1]
        IVW.SE<-EstimatesIVW$coefficients[2]
        IVW_CI<-confint(IVW.Model)
        Qj<-W*(Ratios-IVW.Slope)^2
        Total_Q<-sum(Qj)
        Total_Q_chi<-pchisq(Total_Q,length(temp$beta.exposure)-1,lower.tail = FALSE)
        b <- EstimatesIVW$coefficients[1,1]
        se <- EstimatesIVW$coefficients[1,2]
        pval <- 2 * pnorm(abs(b/se), lower.tail = FALSE)
        Q_df <- DF
        Q <- Total_Q
        Q_pval <- pchisq(Q, Q_df, lower.tail=F)
        nsnp <- nrow(temp)
        res <- data.frame(id.exposure = temp$id.exposure[1], id.outcome = temp$id.outcome[1], outcome = temp$outcome[1],  
                          method = c("Radial MR"), nsnp = nsnp, 
                          b = b, se = se, pval = Q_pval)   
      }
      
      if(nrow(temp) > 1 & nrow(temp) <3)
      {
        res <- TwoSampleMR::mr(temp, method_list=c("mr_ivw"))     
      }
      
      if(nrow(temp) == 1)
      {
        res <- TwoSampleMR::mr(temp, method_list=c("mr_wald_ratio"))     
      }
      
      res_clust[[i]] <- res
    }
    
    return(res_clust)
  }
  
  output$estimate <- mr_cluster_estimate(dat = output$dat)
  return(output)
}


# Heterogeneity

mr_cluster_heterogeneity <- function(cluster = NULL, weights = 3){

  dat <- cluster$dat
  estimate <- cluster$estimate
  
  Q <- 0

  for (i in 1:max(dat[ , ncol(dat)]))
  {
    temp <- subset(dat, dat[ , ncol(dat)] == i)  
    temp$beta <- temp$beta.outcome / temp$beta.exposure
    
    if(nrow(temp) < 1){
      Qsum = 0
    }
    
    if(nrow(temp) > 1){
    #temp$w <- weight[match(temp$SNP, ios_dat$SNP)]
    
    if(weights == 1){
      W <- ((temp$beta.exposure^2) / (temp$se.outcome^2))
    }
    
    if(weights == 2){
      W<-((temp$se.outcome^2/temp$beta.exposure^2)+((temp$beta.outcome^2*temp$se.exposure^2)/temp$beta.exposure^4))^-1
    }
    
    if(weights==3){
      W <- ((temp$beta.exposure^2) / (temp$se.outcome^2))
      Wj<-sqrt(W)
      BetaWj<- (temp$beta.outcome/temp$beta.exposure)*Wj
      IVW.Model<-lm(BetaWj~-1+Wj)
      EstimatesIVW<-summary(lm(IVW.Model))
      IVW.Slope<-EstimatesIVW$coefficients[1]
      W <-  ((temp$se.outcome^2+(IVW.Slope^2*temp$se.exposure^2))/temp$beta.exposure^2)^-1
    bi <- estimate[[i]]$b 
    
    Qj <- (W * ( bi - temp$beta)^2)
    Qsum <- sum(Qj)
    Q <- Q + Qsum
    }
    #Q_pval <- pchisq(Q, nrow(dat)-1, lower.tail = FALSE)
    }
    
  }
  
  return(Q)
  
  }
  
