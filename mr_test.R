if (!'TwoSampleMR' %in% installed.packages()[,1]){
  if(!'devtools' %in% installed.packages()[,1]){
    install.packages('devtools', repos='https://cloud.r-project.org')
  }
  devtools::install_github("MRCIEU/TwoSampleMR")
}
if(!'MRPRESSO' %in% installed.packages()[,1]){
  if(!'devtools' %in% installed.packages()[,1]){
    install.packages('devtools', repos='https://cloud.r-project.org')
  }
  devtools::install_github("rondolab/MR-PRESSO")
}
if (!'data.table' %in% installed.packages()[,1]){
  install.packages("data.table",, repos='https://cloud.r-project.org')
}

args <- commandArgs(T)
mr_all <- function(exp_fp, out_fp, expname, outname, method=c('mr_wald_ratio', 'mr_two_sample_ml','mr_egger_regression', 'mr_weighted_median', 'mr_ivw', 'mr_presso'), pleio_test='no'){
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
  twosamplemr_method = method[method %in% c('mr_presso', 'mr_gsmr')==F]
  if (length(twosamplemr_method)==0){
    res <- data.frame(outcome=character(), exposure=character(), method=character(), nsnp=character(), b=character(), se=character(), pval=character())
  }else{
    res <-  TwoSampleMR::mr(dat, method_list =method[method %in% c('mr_presso', 'mr_gsmr')==F])
    res$method <- chartr(' ', '_', res$method)
    res <- subset(res, select=outcome:pval)
  }

  
  if ('mr_presso' %in% method){
    if (nrow(dat) <= length("beta.exposure") +2){
      cat("Not enough intrumental variables for MR PRESSO analysis")
    }else{
      #presso_res <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = c("se.exposure"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 5000,  SignifThreshold = 0.05)
      if (pleio_test == 'yes'){
        test_flag = TRUE
      }else if (pleio_test == 'no'){
        test_flag = FALSE
      }  
      for (i in 1:10){
        error <- tryCatch(
          presso_res <<- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = c("se.exposure"), OUTLIERtest = test_flag, DISTORTIONtest = test_flag, data = dat, NbDistribution = i*1000,  SignifThreshold = 0.05)
          ,warning=function(e) e
        )
        if (!inherits(error, 'warning')){
          cat(paste(i*1000, ' loops for presso\n'))
          break
        }
      }
      main_presso_res <- presso_res$`Main MR results`
      n_snps <- c(nrow(dat), nrow(dat)-length(presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`))
      colnames(main_presso_res) <- c('exposure', 'method', 'b', 'se', 'T','pval')
      main_presso_res$outcome <- outname
      main_presso_res$nsnp <- n_snps
      main_presso_res$method <- c('presso','presso_outlier_correction')
      main_presso_res$exposure <- expname
      main_presso_res <- main_presso_res[,c(7,1,2,8,3,4,6)]
      res <- rbind(res, main_presso_res)
    }
  }
  write.table(res, file='mr.result', sep = '\t', quote = F,row.names = F, col.names = T, fileEncoding = 'utf-8')
  cat('\n\n--Main MR results--\n')
  print(res)
  if (nrow(dat)==1){
    cat('No enough SNPs for heterogeneity test or pleiotropy test\n')
  }else{
    heter_test<- TwoSampleMR::mr_heterogeneity(dat,method_list = intersect(subset(TwoSampleMR::mr_method_list(), heterogeneity_test)$obj ,method))
    heter_test$method <- chartr(' ', '_', heter_test$method)
    heter_test <- subset(heter_test, select=outcome:Q_pval)
    write.table(heter_test, file='mr.heter', sep = '\t', quote = F,row.names = F, col.names = T, fileEncoding = 'utf-8')
    cat('\n\n--Main heterogeneity test for MR results--\n')
    print(heter_test)
    if ("mr_egger_regression" %in% method ){
      pleio_test<- TwoSampleMR::mr_pleiotropy_test(dat)
      pleio_test<-subset(pleio_test, select=outcome:pval)
      write.table(pleio_test, file='mr.pleio', sep = '\t', quote = F, row.names = F, col.names = T, fileEncoding = 'utf-8')
      cat('\n\n--Egger pleiotropy: Intercetpt--\n')
      print(pleio_test)
    }
  }
  write.table(dat, file='mr.data', sep='\t', quote=F, row.names=F, col.names=T, fileEncoding='utf-8')
  return(dat)
}

method = unlist(strsplit(args[5], split=','))
dat <- mr_all(args[1], args[2], args[3] , args[4], method=method, args[6])
quit(save='no')



