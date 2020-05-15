extract_phewas <- function(id_exp = id_exp, id_out = id_out, exp = exp_dat, id_bg = id_bg){
  #split SNPs
  n <- length(exp$SNP)
  snplist <- exp$SNP
  splits <- data.frame(exp$SNP, chunk_id=rep(1:(ceiling(n/50)), each=50)[1:n])
  d <- list()
    for(i in 1:n)
    {
      message(i, "of", n, "SNPs")
      
      d[[i]] <- plyr::ddply(splits, c("chunk_id"), function(x)
      {
        x <- plyr::mutate(x)
        message(" [>] ", x$chunk_id[1], " of ", max(splits$chunk_id), " chunks")
        
        if(any(id_bg %like% c("ebi")))
        {
          #id <- grep("ebi", id_bg, value = T)
          out1 <- phewas(snplist, pval=1, batch="ebi-a", 
                         access_token = check_access_token())
        }
        
        if(any(id_bg %like% c("ieu-a")))
        {
          out1 <- phewas(snplist, pval=1, batch="ieu-a", 
                        access_token = check_access_token())
        }
        
        if(any(id_bg %like% c("ukb-b")))
        {
          out2 <- phewas(snplist, pval=1, batch="ukb-b", 
                       access_token = check_access_token())
        }
        
        if(any(id_bg %like% c("eqtl-a")))
        {
          out3 <- phewas(snplist, pval=1, batch="eqtl-a", 
                       access_token = check_access_token())
        }
        
        if(any(id_bg %like% c("met")))
        {
          meta <- phewas(snplist, pval=1, batch="met-a", 
                       access_token = check_access_token())
          metb <- phewas(snplist, pval=1, batch="met-b", 
                         access_token = check_access_token())
          metc <- phewas(snplist, pval=1, batch="met-c", 
                        access_token = check_access_token())
          out4 <- rbind(meta, metb, metc)
        }
        if(any(id_bg %like% c("prot")))
        {
          prota <- phewas(snplist, pval=1, batch="prot-a", 
                         access_token = check_access_token())
          protb <- phewas(snplist, pval=1, batch="prot-b", 
                         access_token = check_access_token())
          out5 <- rbind(prota, protb)
        }
        if(any(id_bg %like% c("ubm")))
        {
          out6 <- phewas(snplist, pval=1, batch="ubm-a", 
                       access_token = check_access_token())
        }
        
        df_list <- mget(ls(pattern = "^out\\d+"))
        #out <- plyr::rbind(df_list)
        temp <- bind_rows(df_list)
        
        #Remove traits 
        temp <- subset(temp) %>%
          filter(! temp$id %in% c(id_exp, id_out)) %>%
          filter(test$id.outcome %in% c(id_bg)) %>%
          arrange(desc(sample_size)) %>%

        
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
      })
    }
  d <- bind_rows(d)
  return(d)
}
   
  
