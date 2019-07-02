# Alyssa Holman - Last Updated: 6/27/2019
# Script used for a GWAS
# # Here we compare phenotypes to SNPs to find SNPs of interest

# dgrp2 files from http://dgrp2.gnets.ncsu.edu/
# downloaded Plink formatted genotype data

# ran files through Plink with filtering to get SNPs
# gives dgrp2.traw file
/programs/plink-1.9-x86_64-beta3.46/plink --bfile dgrp2 --maf 0.05 --geno 0.2 --indep-pairwise 50 5
0.5 --out gwas_plink

/programs/plink-1.9-x86_64-beta3.46/plink --bfile dgrp2 --maf 0.05 --geno 0.2 --extract
gwas_plink.prune.in --assoc --make-bed --recode A-transpose --out dgrp2

# take phenotype data and Plink data into R

# uses lme4 for analysis
library(lme4)

# loads snp data generated from Plink
snpdat = read.delim(file="dgrp2.traw",as.is=T)
snp.m <- as.matrix(snpdat[, -c(1:6)])
col.line.vec <- as.integer(sub("^line_(\\d+)_.*$", "\\1", colnames(snp.m)))

# load phenotype data from txt file
dat = read.delim(file="25BM.pheno.txt", as.is=T)

# association test using cores
test.out <- mclapply(1:nrow(snp.m))
function(snp.idx)
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
  
  # specifies cores used for this
}, mc.cores=10)

# result is in test.out, as a list from lapply()
test.out.vec <- unlist(test.out)
save(test.out.vec, file="test_out_vec.rdata")

# making a manhattan plot
manhat.df <- data.frame(chr=snpdat[,1], pos=snpdat[,4], pval=test.out.vec, stringsAsFactors=F)
source("manhattan_plot.r")
manhattan.plot(manhat.df)

# chromosomes to numbers
# 1 = 2L, 2 = 2R, 3 = 3L, 4 = 3R, 5 = X, 6 = 4

# adding abline at e-05
abline(h=5, col="red", lty=2)

# snp hits
manhat.df$adjusted.pval <- p.adjust(manhat.df$pval, method="BH")

# hits with adjusted p-value < 0.05
snp.hits <- manhat.df[which(manhat.df$adjusted.pval < 0.05), ]
