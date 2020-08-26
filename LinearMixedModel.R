# Linear Mixed Model in R
# Last Edited: May 10, 2019 by Alyssa Holman

# Set directory and read in data with fly lines
setwd("/local/workdir/arh223/GWAS")
dat = read.delim(file="25BM.pheno.txt", as.is=T)

# Load Libraries
library(parallel)
library(lme4)

# Read in .traw file after processing in Plink
snpdat = read.delim(file="dgrp2.traw",as.is=T)
snp.m <- as.matrix(snpdat[, -c(1:6)])
col.line.vec <- as.integer(sub("^line_(\\d+)_.*$", "\\1", colnames(snp.m)))

# Run linear mixed model
ptm <- proc.time()
test.out <- mclapply(1:nrow(snp.m), function(snp.idx)
{
  c.dat <- dat
  c.dat = c.dat[c.dat$line %in% col.line.vec, ]
  c.dat$gt <- snp.m[snp.idx, match(c.dat$line, col.line.vec)]
  
  c.dat <- c.dat[which(!is.na(c.dat$gt)), ]
  
  lmer.out <- lmer(area ~ gt + (1 | line), data=c.dat, REML=F)
  null.out <- lmer(area ~ (1 | line), data=c.dat, REML=F)
  a1 <- anova(null.out, lmer.out)
  p1 <- a1[["Pr(>Chisq)"]][2]
  
  return(p1)
}, mc.cores=10)

proc.time() - ptm

# Results from model in "text.out"
test.out.vec <- unlist(test.out)
save(test.out.vec, file="test_out_vec.rdata")
 
# Load back into R to plot and extract SNP data
load("test_out_vec.rdata")

# Make a Manhattan Plot
# Using "manhattan_plot.r" file
manhat.df <- data.frame(chr=snpdat[,1], pos=snpdat[,4], pval=test.out.vec, stringsAsFactors=F)
source("manhattan_plot.r")
manhattan.plot(manhat.df)
abline(h=5, col="red", lty=2)

# Get SNP hits
manhat.df$adjusted.pval <- p.adjust(manhat.df$pval, method="BH")
snp.hits <- manhat.df[which(manhat.df$adjusted.pval < 0.05), ]
snp.hits

