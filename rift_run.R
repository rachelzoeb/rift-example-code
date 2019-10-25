### Code to run SNP Influential Filtering Tool
### Last updated September 17, 2019

## set directory 
setwd("~/rift_example_code/")

## Load functions 
source("rift_functions.R")

### Read in example dataset (in PLINK format using --recodeA)
raw_dat <- read.table(file = "rift_example_data.txt")

## Run function to calculate delta chi-square scores for each variant
## This function assumes data formatted as sift_example_data.txt (output from PLINK --recodeA)
rift.results <- rift_skato(raw_dat)

## Sort results by delta chi-square score
rift.results$ind.stats[order(rift.results$ind.stats$delta, decreasing = FALSE),]

## Calculate Tukey Fences to identify IV
rift.results <- calculateTukeyFences(rift.results)
rift.results$ind.stats$IV <- ifelse(rift.results$ind.stats$tukey.mild, "Yes", "No")

## plot delta value by position
# extract position from SNP name
rift.results$ind.stats$pos <- as.numeric(sapply(as.character(rift.results$ind.stats$SNP_excluded), FUN = function(x){strsplit(x, split = "[.]")[[1]][2]}))

# plot delta by position
plot(rift.results$ind.stats$delta ~ rift.results$ind.stats$pos, col = ifelse(rift.results$ind.stats$IV=="Yes","red","black"), pch = 19, xlab = "Position", ylab = "Delta chi-square score")
abline(h = 0, col = "red")
