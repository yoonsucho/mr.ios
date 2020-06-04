#'Extract information from IEU GWAS database
#'
#'Returns traits which will be used in making background dataset. The data is obtained from [IEU GWAS database](https://gwas.mrcieu.ac.uk/).
#' 
#' @param snplist A list of instrumental genetic variants for the exposure.
#' @param id_bg A list of available GWAS summary statistics.
#' @param nsnp_per_chunk Number of SNPs to be used per phewas scanning.
#' 
#' 
#' @export
#' @return Data frame

extract_phewas <- function(snplist = NULL, id_bg = id_bg, nsnp_per_chunk = 5){
  
  batchlist <- sapply(strsplit(id_bg, "-"), function(x) paste(x[1], x[2], sep="-"))
  #basplit <- rep(1:ceiling((length(batchlist)/1500)), each = 1500)[1:length(batchlist)]
  #balist_split <- split(batchlist, basplit)
  #basplits <- data.frame(id = batchlist, chunk_id=basplit)
  
  n = length(snplist)
  nsplit <- rep(1:(ceiling(n/nsnp_per_chunk)), each = nsnp_per_chunk)[1:n]
  snplist_split <- split(snplist, nsplit)
  splits <- data.frame(snps = snplist, chunk_id=nsplit)
  
  
  #if(max(basplits$chunk_id) > 1)
  #  {
    l <- list()
    for(i in 1:length(nsnp_per_chunk))
    {
      message(max(splits$chunk_id), " chunks were generated out of ", n, " SNPs")
      
      l[[i]] <- plyr::ddply(splits, c("chunk_id"), function(x)
      {
        x <- plyr::mutate(x)
        message(" [>] ", x$chunk_id[1], " of ", max(splits$chunk_id), " chunks")
        
        #for(j in 1:max(basplits$chunk_id))
        #  {
           d <-  ieugwasr::phewas(x$snps, pval=1, batch=unique(batchlist)) %>% subset(id %in% id_bg)
        #  }
        
        if(!is.data.frame(d)) d <- data.frame()
        return(d)
      })
    }
    temp <- bind_rows(l)
    
  #}
    
  
  #remove duplicates
  #remove traits without eaf info
  #remove traits with se < 0

  temp <- subset(temp) %>%
    filter(!duplicated(temp)) %>%
    filter(!is.na(eaf)) %>%
    filter(se > 0)

  
  #Rearrange dataset
  col_order <- c("rsid", "chr", "position", "beta", "se", "n", "p", "eaf", "ea", "nea", "trait", "id")
  out <- temp[, col_order]
  
  out <- plyr::rename(out, c("rsid" = "SNP",
                             "position" = "pos",
                             "beta" = "beta.outcome",
                             "se" = "se.outcome",
                             "n" = "samplesize.outcome",
                             "p" = "pval.outcome",
                             "eaf" = "eaf.outcome",
                             "ea" = "effect_allele.outcome",
                             "nea" = "other_allele.outcome",
                             "trait" = "outcome",
                             "id" = "id.outcome"))
  
  
  
  if(!is.data.frame(out)) out <- data.frame()
  return(out)
  
}
