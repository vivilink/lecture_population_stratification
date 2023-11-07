#ensure that cmake is installed on your computer
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("snpStats")

library("snpStats")
setwd("/data/class_stratification/simulations2")


simulate_phenotypes_genetic_confounding <- function(genotypes_polymorphic, summary_stats, num_causal_variants, heritability){
  variants_causal <- sample(x=colnames(genotypes_polymorphic), size=num_causal_variants)
  allele_freq_variants_causal <- summary_stats[variants_causal, colnames(summary_stats) == "MAF"]

  # simulate causal effect sizes according to allele frequency 
  heritability <- 0.7
  sd = (heritability / num_causal_variants) / sqrt(2 * allele_freq_variants_causal * (1 - allele_freq_variants_causal))
  # sd <- 10
  effect_sizes_causal = rnorm(mean=0, sd=sd, n=num_causal_variants)
  
  # assign effect sizes 
  effect_sizes <- numeric(num_polymorphic_snps)
  names(effect_sizes) <- colnames(genotypes_polymorphic)
  effect_sizes[variants_causal] <- effect_sizes_causal 
  
  
  # multiply effect sizes with genotypes
  genotypes_polymorphic_raw <- as.matrix(genotypes_polymorphic@.Data)
  genotypes_polymorphic_raw <- matrix(as.numeric(as.vector(genotypes_polymorphic_raw)), nrow = nrow(genotypes_polymorphic_raw))
  colnames(genotypes_polymorphic_raw) <- colnames(genotypes_polymorphic)
  genotypes_polymorphic_raw <- genotypes_polymorphic_raw - 1
  phenotypes <- genotypes_polymorphic_raw[,variants_causal] %*% effect_sizes[variants_causal]
  
  # add noise to get requested heritability
  var_genetic <- var(phenotypes)
  var_env <- var_genetic * (1 - heritability) / heritability
  phenotypes <- phenotypes + rnorm(n=num_ind, mean=0, sd=sqrt(var_env))
  
  return(phenotypes)
}


simulate_phenotypes_environmental_confounding <- function(name_fam, offset){
  fam <- read.table(name_fam)
  colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
  num_inds <- length(fam$FID)
  phenotypes <- rnorm(n=num_inds, mean=300, sd=100)
  names(phenotypes) <- fam$IID
  TSI <- fam$IID[fam$FID == "TSI"]
  phenotypes[TSI] <- phenotypes[TSI] + offset
  
  return(phenotypes)
}

# extract CDX and TSI from total data
# system(paste("plink2 --bfile allChroms_ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes --maf 0.05 --geno 0.1 --hwe 1e-6 --mind 0.1 --keep sample_names_populations.txt --make-bed --out subset", sep=''))


# plink_prefix <- "/data/class_stratification/Data/CDX_10kSnps_ldPruned/CDX_unrels"
# plink_prefix <- "/data/class_stratification/Data/FullDataLDPruned/allChroms_ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes"
plink_prefix <- "subset"

# separate plink files by chromosomes
system(paste("plink2 --bfile ",plink_prefix, " --chr 1-5 --make-bed --out ", plink_prefix, "_firstFiveChr", sep=''))
system(paste("plink2 --bfile ",plink_prefix, " --chr 6-10 --make-bed --out ", plink_prefix, "_allButFirstFiveChr", sep=''))


# read data
# data <- read.plink(paste(plink_prefix, "_firstFiveChr", sep=''))
data <- read.plink(plink_prefix)

show(data$genotypes)

# access genotypes
num_snp <- ncol(data$genotypes)
num_ind <- nrow(data$genotypes)

# look at data
summary_stats <- col.summary(data$genotypes)
summary_stats

#? what are the meanings of the different columns? Hint: Check the snpStat manual section row.summary

#? Why are the z.HWE values for some SNPs NA?
##because one of the alleles has fixed, these SNPs are actually not polymorphic in this sample

#? What is the meaning of RAF?
##Rist allele frequency, the frequency of allele B

# only keep polymorphic sites
polymorphic_snps <- row.names(summary_stats[summary_stats$MAF > 0.01,])
num_polymorphic_snps <- length(polymorphic_snps)
genotypes_polymorphic <- data$genotypes[,polymorphic_snps]

#add population label to original fam
# pop_info <- read.csv("/data/class_stratification/Data/FullDataLDPruned/provided/sample_names.csv", header=TRUE, sep=',')
# colnames(pop_info) <- c("IID", "FID")
# fam <- read.table(paste(plink_prefix, ".fam", sep=''), header=FALSE)
# colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
# fam <- merge(pop_info, fam, by="IID")
# fam <- fam[,-which(colnames(fam) == "FID.y")]
# colnames(fam)[colnames(fam) == "FID.x"] <- "FID"
# fam <- cbind(fam$FID, fam$IID, fam$PID, fam$MID, fam$Sex, fam$Phenotype)
# write.table(x=fam, file="subset.fam", row.names = FALSE, sep=" ", quote=FALSE, col.names = FALSE)

# phenotypes <- simulate_phenotypes_genetic_confounding(genotypes_polymorphic, summary_stats, num_causal_variants=200, heritability=0.7)
phenotypes <- simulate_phenotypes_environmental_confounding(name_fam=paste(plink_prefix, ".fam", sep=''), offset=1000)

# add phenotype to plink .fam file
fam <- read.table(paste(plink_prefix, ".fam", sep=''))
fam[,6] <- phenotypes
write.table(fam, file=paste(plink_prefix, ".fam", sep=''), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

