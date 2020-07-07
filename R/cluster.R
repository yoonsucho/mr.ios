
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
kmeans_instruments <- function(bg_dat, value_column, kmax=15, nstart=50, iter.max=15)
{
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
	which.max(diff_mse) + 2
	dat <- dplyr::tibble(instrument=rownames(cluster), cluster=cluster[,which.max(diff_mse) + 2])
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
hclust_instruments <- function(bg_dat, value_column, kmax=50)
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
	cuts <- cutree(hc, k=1:50)
	mse <- get_mse(d, cuts)
	diff_mse <- diff(mse) * -1
	which.max(diff_mse) + 2
	dat <- dplyr::tibble(instrument=rownames(cuts), cluster=hc[[which.max(diff_mse) + 2]]$cluster)
	return(dat)
}

#' Cluster
#'
#' <full description>
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
	wide <- bg_to_wide(bg_dat, value_column)
	pc <- prcomp(wide)
	fit <- pvclust::pvclust(t(pc$x), method.hclust=method.hclust, method.dist=method.dist, nboot=nboot)
	plot(fit)
	pvclust::pvrect(fit, alpha=alpha)
}


