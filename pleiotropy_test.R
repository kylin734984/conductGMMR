if(!'MRPRESSO' %in% installed.packages()[,1]){
  if(!'devtools' %in% installed.packages()[,1]){
    install.packages('devtools', repos='https://cloud.r-project.org')
  }
  devtools::install_github("rondolab/MR-PRESSO")
}
if (!'data.table' %in% installed.packages()[,1]){
  install.packages('data.table', repos='https://cloud.r-project.org')
}
args <- commandArgs(T)

pleiotropy <- function(exp_fp, out_fp, expname, outname, method='mr_presso'){
  exp_df <- data.table::fread(exp_fp, header=T, stringsAsFactors = F, col.names= c('SNP', 'effect_allele','other_allele', 'eaf', 'beta','se', 'pval', 'samplesize'))
  exp_dat <- TwoSampleMR::format_data(exp_df, type="exposure")
  out_df <- data.table::fread(out_fp, header=T, stringsAsFactors = F, col.names= c('SNP', 'effect_allele','other_allele', 'eaf', 'beta','se', 'pval', 'samplesize'))
  out_dat <-  TwoSampleMR::format_data(out_df, type="outcome")
  dat <-  TwoSampleMR::harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = out_dat,
    action=3
  )
  dat$exposure=expname
  dat$outcome=outname
  if (method == 'mr_presso' ){
    if (nrow(dat)<=3){
      cat('Not enough SNPs to perform MR-PRESSO global test.\n')
    }else{
      for (i in 2:10){
        error <- tryCatch(
          presso_res <<- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = c("se.exposure"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = i*1000,  SignifThreshold = 0.05)
          ,warning=function(e) e
        )
        if (!inherits(error, 'warning')){
          print(paste(i*1000, ' loops for presso'))
          break
        }
      }
      pleio_indices = presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      #nonpleio <- setdiff(row.names(dat), pleio_indices)
      #pleio <- intersect(row.names(dat), pleio_indices) 
      #print(warning(pleio_indices))
	  #print(pleio_indices)
	  #print(warning(presso_res))
	  #print(presso_res)
      if (is.null(pleio_indices)){
        cat(paste('No pleiotropic SNPs for ', expname, ' and ', outname, '!\n'))
        pleiodf <- data.frame(exposure=character(), outcome=character(), method=character(), RS=character())  
      }else if (pleio_indices == 'No significant outliers'){
        cat(paste('No pleiotropic SNPs for ', expname, ' and ', outname, '!\n'))
        pleiodf <- data.frame(exposure=character(), outcome=character(), method=character(), RS=character())  
      }else{
        pleiosnps <- dat[intersect(row.names(dat), pleio_indices),]$SNP
        pleiodf <- data.frame(exposure=rep(expname, length(pleiosnps)), outcome=rep(outname, length(pleiosnps)), method=rep('mr_presso',length(pleiosnps)), RS=pleiosnps)
        cat(paste('Pleiotropic SNPs for ', expname, ' and ', outname, 'using MR-PRESSO global test!\n'))
        print(pleiodf)
      }
      write.table(pleiodf, file='mr.pleiosnp', sep = '\t', quote = F,row.names = F, col.names = T, fileEncoding = 'utf-8')
      nonpleiosnps <- dat[setdiff(row.names(dat), pleio_indices),]$SNP
      write.table(nonpleiosnps, file = 'include', sep = '\t', quote=F, row.names = F, col.names = F, fileEncoding = 'utf-8')
      #重新写exposure.txt和outcome.txt
      new_exposure = exp_df[exp_df$SNP %in% nonpleiosnps,]
      colnames(new_exposure) <- c('RS', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P', 'N')
      new_outcome  = out_df[out_df$SNP %in% nonpleiosnps,]
      colnames(new_outcome) <- c('RS', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P', 'N')
      write.table(new_exposure, file= 'exposure.txt', sep = '\t', quote=F, row.names = F, col.names = T, fileEncoding = 'utf-8')
      write.table(new_outcome, file= 'outcome.txt', sep = '\t', quote=F, row.names = F, col.names = T, fileEncoding = 'utf-8')
      
    }
    
  }
}

if (F){
  args[1] <- 'exposure.txt'
  args[2] <- 'outcome.txt'
  args[3] <- 'exposure'
  args[4] <- 'outcome'
  args[5] <- 'mr_presso'
}
pleiotropy(args[1], args[2], args[3], args[4], args[5])


