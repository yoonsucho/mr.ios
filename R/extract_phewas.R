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

extract_phewas <- function(snplist = NULL, id_bg = id_bg, nsnp_per_chunk = 10){
  
  #Define the list of batches for Phewas
  batchlist <- sapply(strsplit(id_bg, "-"), function(x) paste(x[1], x[2], sep="-"))
  batches <- unique(batchlist)

  #Split the SNPs into each chunk
  n = length(snplist)
  nsplit <- rep(1:(ceiling(n/nsnp_per_chunk)), each = nsnp_per_chunk)[1:n]
  snplist_split <- split(snplist, nsplit)
  splits <- data.frame(snps = snplist, chunk_id=nsplit)
  
  
  #Phewas scanning of SNPs in each chunk and each batch
    l <- list()
   
    k <- 1
    for(i in 1:max(nsplit))
      {
      message(max(nsplit), " chunks were generated out of ", n, " SNPs")
      for(j in batches)
        {
        message(" [>] ", i, " of ", max(splits$chunk_id), " chunks; search for ", j)
        l[[k]] <- ieugwasr::phewas(variants=snplist_split[[i]], batch=j) %>% subset(id %in% id_bg)
        k <- k + 1
         }
      }
      
    temp <- dplyr::bind_rows(l)


  #remove duplicates
  #remove traits without eaf info
  #remove traits with se < 0

  temp <- subset(temp) %>%
    dplyr::filter(!duplicated(temp)) %>%
    dplyr::filter(!is.na(eaf)) %>%
    dplyr::filter(se > 0)

  
  #Rearrange dataset
  col_order <- c("rsid", "chr", "position", "beta", "se", "n", "p", "eaf", "ea", "nea", "trait", "id")
  out <- temp[, col_order]
  
  out <- dplyr::rename(out, c("SNP" = "rsid",
                              "pos" = "position",
                              "beta.outcome" = "beta",
                              "se.outcome" = "se",
                              "samplesize.outcome" = "n",
                              "pval.outcome" = "p",
                              "eaf.outcome" = "eaf",
                              "effect_allele.outcome" = "ea",
                              "other_allele.outcome" = "nea",
                              "outcome" = "trait",
                              "id.outcome" = "id"))
  
  
  
  if(!is.data.frame(out)) out <- data.frame()
  return(out)
  
}
