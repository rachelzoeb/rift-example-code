### Functions for perform Rare Variant Influential Filtering Tool

### Library loads
library(SKAT)
library(dplyr)

## Function which currently takes in the raw_data file output from PLINK with both genotype and phenotype information
rift_skato <- function(raw_dat, covariates=NULL){
  
  ## Load functions
  require(SKAT)
  require(dplyr)
  
  ## Save phenotype information
  PHENOTYPE <- raw_dat$PHENOTYPE - 1
  #COV <- raw_dat$SEX
  GENO <- as.matrix(raw_dat[,!(names(raw_dat) %in% c("FID", "IID", "PAT", "MAT", "PHENOTYPE", covariates))])
  SNPS <- names(raw_dat)[!(names(raw_dat) %in% c("FID", "IID", "PAT", "MAT", "PHENOTYPE", covariates))]
  
  if (is.null(covariates) == TRUE){
    
    obj <- SKAT_Null_Model(PHENOTYPE ~ 1, out_type = "D", Adjustment = FALSE) 
    
  } else {
    
    COV <- raw_dat[,covariates]
    obj <- SKAT_Null_Model(PHENOTYPE ~ as.matrix(COV), out_type = "D", Adjustment = FALSE) 
    
  }
  
  ## Create dataframe to hold results
  res <- NULL
  
  ## Loop through to remove one SNPS and re-run SKAT, saving results in iterative_results dataframe
  for (idx in c(0:dim(GENO)[2])){
    
    ## Set up genotype data to remove 0 or 1 snps
    if (idx == 0){
      
      GENO_sub <- GENO
      SNP_excluded <- "NONE"
      
      ## Run SKATO to find optimal parameter
      out_skato <- SKAT(GENO_sub, obj, method = "optimal.adj")
      rho_est <- out_skato$param$rho_est
      pval_rho <- out_skato$param$minp
      
    } else {
      
      GENO_sub <- GENO[,-idx]
      SNP_excluded <- SNPS[idx]
      
      ## Run SKATO with rho from entire data 
      out_skato <- SKAT(GENO_sub, obj, r.corr = rho_est)
      pval_rho <- out_skato$p.value
      
    }
    
    
    ## Save results to dataframe
    res <- rbind(data.frame(SNP_idx = idx, SNP_excluded = SNP_excluded, p.value = pval_rho, rho_est = rho_est), res)
    
  }
  
  ## Calculate Chi-Square statistics 
  res$chisq <- sapply(res$p.value, FUN = function(x){return(qchisq(x, df = 1, lower.tail = FALSE))})
  
  ## Calculate delta Chi-Square statistics
  res$delta <- sapply(res$chisq, FUN = function(x){return(x - res$chisq[res$SNP_excluded == "NONE"])})
  
  ## Resort results
  res <- arrange(res, SNP_idx)
  
  ## Return results 
  return(list(ind.stats = res))
  
}

# Function to calculate Tukey fences 
calculateTukeyFences <- function(results){
  
  # Extract full parameter
  full_delta <- results$ind.stats$delta[results$ind.stats$SNP_excluded == "NONE"]
  
  # Extract jackknife etsimates
  ind_delta <- results$ind.stats$delta[results$ind.stats$SNP_excluded != "NONE"]
  n <- length(ind_delta)
  
  # Calculate IQR and quantiles
  lowerq = quantile(ind_delta)[2]
  upperq = quantile(ind_delta)[4]
  iqr = upperq - lowerq
  
  # Compute bounds for inner fences
  mild.threshold.upper = (iqr * 1.5) + upperq
  mild.threshold.lower = lowerq - (iqr * 1.5)
  
  # Compute bounds for outer fences
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  
  # Save thresholds
  results$iqr <- list(iqr = iqr, lowerq = lowerq, upperq = upperq, 
                      mild.threshold.lower = mild.threshold.lower, mild.threshold.upper = mild.threshold.upper, 
                      extreme.threshold.lower = extreme.threshold.lower, extreme.threshold.upper = extreme.threshold.upper)
  
  # Identify outlier
  results$ind.stats$fence <- sapply(results$ind.stats$delta, function(x){
    if ((x < extreme.threshold.lower) | (x > extreme.threshold.upper)) {
      return("extreme")
    } else if ((x < mild.threshold.lower) | (x > mild.threshold.upper)) {
      return("mild") 
    } else {
      return("none")
    }
  })
  
  # Make the fence variable a factor
  results$ind.stats$fence <- factor(results$ind.stats$fence, levels = c("extreme", "mild", "none"))
  
  # Make test vectors 
  # results$ind.stats$tukey.mild <- ifelse(results$ind.stats$fence %in% c("mild", "extreme"), TRUE, FALSE)
  results$ind.stats$tukey.mild <- ifelse((results$ind.stats$fence %in% c("mild", "extreme")) & (results$ind.stats$delta < 0), TRUE, FALSE)
  results$ind.stats$tukey.extreme <- ifelse((results$ind.stats$fence %in% c("extreme")) & (results$ind.stats$delta < 0), TRUE, FALSE)
  
  # Return thresholds
  return(results)
  
}