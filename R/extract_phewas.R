extract_phewas <- function(snplist = NULL, id_bg = id_bg, nsnp_per_chunk = 30){
  
  batchlist <- sapply(strsplit(id_bg, "-"), function(x) paste(x[1], x[2], sep="-"))
  
  n <- length(snplist)
  nsplit <- round(length(snplist) / nsnp_per_chunk)
  snplist_split <- split(snplist, 1:nsplit)
  
  splits <- data.frame(snplist, chunk_id=rep(1:(ceiling(nsplit)), each=nsnp_per_chunk)[1:n])

  
  l <- list()
  for(i in 1:length(nsnp_per_chunk))
  {
    message(nsplit, " chunks were generated out of ", length(snplist), " SNPs")
    
    l[[i]] <- plyr::ddply(splits, c("chunk_id"), function(x)
      {
      x <- plyr::mutate(x)
      message(" [>] ", x$chunk_id[1], " of ", max(splits$chunk_id), " chunks")
      
      ieugwasr::phewas(snplist_split[[i]], pval=1, batch=unique(batchlist)) %>% subset(id %in% id_bg)
    })
  }
  temp <- bind_rows(l)
  
  
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
