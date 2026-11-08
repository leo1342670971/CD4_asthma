rm(list=ls())
library(tidyverse)
library(data.table)

library(coloc)
library(pacman)
library(TwoSampleMR)
library(biomaRt)
p_load(data.table,dplyr,tidyr)

WINDOW_BP <- 1000000
RUN_ID <- paste0(gsub("[^0-9_]", "", format(Sys.time(), "%Y%m%d_%H%M%OS6")), "_", Sys.getpid())
OUT_DIR <- file.path("", "", RUN_ID)
OUT_FILE <- file.path(OUT_DIR, "GBMI_asthma.txt")
SKIPPED_FILE <- file.path(OUT_DIR, "GBMI_asthma.txt")

reserve_output_directory <- function(directory) {
  if (file.exists(directory) || dir.exists(directory))
    stop("Output directory already exists; refusing to overwrite: ", directory)
  parent <- dirname(directory)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE) && !dir.exists(parent))
    stop("Cannot create output parent directory: ", parent)
  if (!dir.create(directory, recursive = FALSE, showWarnings = FALSE))
    stop("Cannot reserve a new output directory: ", directory)
  invisible(directory)
}
reserve_output_directory(OUT_DIR)
message("Output directory: ", OUT_DIR)

window_chr <- function(x) {
  x <- sub("^CHR", "", toupper(trimws(as.character(x))))
  x[x == "M" & !is.na(x)] <- "MT"
  x[!x %in% c(as.character(1:22), "X", "Y", "MT")] <- NA_character_
  x
}
window_position <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)), fixed = TRUE)))
}

gene_window <- function(a, flank_bp = WINDOW_BP) {
  stopifnot(length(flank_bp) == 1L, is.finite(flank_bp), flank_bp >= 0)
  chr <- window_chr(a$GeneChr)
  pos <- window_position(a$GenePos)
  valid <- !is.na(chr) & is.finite(pos) & pos >= 1 & pos == floor(pos)
  annotation <- unique(data.frame(chr = chr[valid], pos = pos[valid]))

  list(ok = TRUE, chr = annotation$chr[1L], center = annotation$pos[1L],
       lower = max(1, annotation$pos[1L] - flank_bp),
       upper = annotation$pos[1L] + flank_bp, flank_bp = flank_bp)
}

in_gene_window <- function(chr, pos, window) {
  stopifnot(isTRUE(window$ok), length(chr) == length(pos))
  chr <- window_chr(chr)
  pos <- window_position(pos)
  valid <- !is.na(chr) & is.finite(pos) & pos >= 1 & pos == floor(pos)
  keep <- rep(FALSE, length(pos))
  keep[valid] <- chr[valid] == window$chr &
    pos[valid] >= window$lower & pos[valid] <= window$upper
  keep
}

result <- data.frame(
  source = character(), eqtl_gene = character(), gene = character(), id = character(),
  p12 = numeric(), nsnps = numeric(), PP.H0.abf = numeric(), PP.H1.abf = numeric(),
  PP.H2.abf = numeric(), PP.H3.abf = numeric(), PP.H4.abf = numeric(), topsnp = character(),
  gene_chr = character(), gene_pos = numeric(), window_lower = numeric(),
  window_upper = numeric(), window_bp = numeric()
)
skipped <- data.frame(id = character(), reason = character())

# p1=p2=1e-4，change p12；

 p12_values <- c(1e-5, 5e-5, 5e-6)
gt <-fread("", header = T)

###CELL_TYPE	SNP	SNPID	type	id	chr	pos	other allele	effect allele	eaf	A2_FREQ_HRC	beta	S_STATISTICS	pval	Q_VALUE	FDR	RSQUARE	GENOTYPED	ROUND	se	samplesize
#gt=separate(gt,colnames(gt),into = c('typename','snp','snpID','type','id','chr','pos','ref','alt','eaf','hrc','beta',
# 'STATISTICS','pvalue','qvalue','FDR','RSQUARE','GENOTYPED','ROUND','se','samplesize'),sep = ',')
#effect allele对应alt 
colnames(gt)<-c("pval","snp","chr","pos","alt","ref","Zscore","geneid","id","GeneChr","GenePos","NrCohorts","samplesize","FDR","BonferroniP","eaf","beta","se")
gt$beta <- as.numeric(trimws(gt$beta))
gt$beta=as.numeric(gt$beta)
gt$pval=as.numeric(gt$pval)
gt$eaf=as.numeric(gt$eaf)
gt$pos=as.numeric(gt$pos)
gt$maf = ifelse(gt$eaf < 0.5, 
                gt$eaf,
                1 - gt$eaf)
gt$se=as.numeric(gt$se)
#gt$samplesize=91
gt$samplesize=as.numeric(gt$samplesize)


asthma_gwas_3<-fread("/data/Asthma_GBMI.txt", header = T)
head(asthma_gwas_3)
asthma_gwas_3=as.data.frame(asthma_gwas_3)

setnames(asthma_gwas_3,
         old = c("other allele", "effect allele"),
         new = c("ref", "alt"))
df2=asthma_gwas_3
df2 <- data.frame(
  SNP = df2$SNP,
  chrom = df2$chr,
  pos = df2$pos,
  A1 = df2$alt,      
  A2 = df2$ref,
  beta = as.numeric(df2$beta),
  se = as.numeric(df2$se),
  pval=df2$pval
)
epositions<-fread("", header = T)
head(epositions)
setDT(epositions)

head(epositions)
epositions[, gene_id := sub("_.*$", "", id.exposure)]

epositions[, gene_symbol := sub("^[^_]+_", "", id.exposure)]

head(epositions[, .(id.exposure, gene_id, gene_symbol)])
epositions_filtered <- epositions[
  (!is.na(FDR) & FDR <= 0.05) &
    id.outcome == "Asthma_GBMI"
]

epositions1=epositions_filtered
#type_list <- unique(epositions_filtered$type)


#for (tp in type_list) {
#for(e in 1:length(d)){
#  gt1=gt[gt$type==tp,]
#  epositions1=epositions[epositions$type==tp,]

genelist1 <- unique(epositions1$gene_id[
  !is.na(epositions1$gene_id) & nzchar(epositions1$gene_id)
])
gt1=gt

if (length(genelist1) == 0) {
  
} else {
  for(i in seq_along(genelist1)){
    a=gt1[gt1$geneid==genelist1[i],]
    if (length(a$snp) == 0) {
      
      reason_now <- "No eQTL data are available for this gene"
      warning(paste(genelist1[i], reason_now), call. = FALSE)
      skipped <- rbind(skipped, data.frame(id = as.character(genelist1[i]), reason = reason_now))
    } else {
   
      window <- gene_window(a)
      if (!isTRUE(window$ok)) {
        reason_now <- window$reason
        warning(paste(genelist1[i], reason_now), call. = FALSE)
        skipped <- rbind(skipped, data.frame(id = as.character(genelist1[i]), reason = reason_now))
        next
      }
      a <- a[in_gene_window(a$chr, a$pos, window), ]
      df2_window <- df2[in_gene_window(df2$chrom, df2$pos, window), , drop = FALSE]
      if (nrow(a) == 0L || nrow(df2_window) == 0L) {
        reason_now <- paste0("no available SNP：eQTL=",
                             nrow(a), "; GWAS=", nrow(df2_window))
        warning(paste(genelist1[i], reason_now), call. = FALSE)
        skipped <- rbind(skipped, data.frame(id = as.character(genelist1[i]), reason = reason_now))
        next
      }
      df1=data.frame(SNP=a$snp,chrom=a$chr,
                     pos=a$pos,A1=a$alt,A2=a$ref,beta=a$beta,se=a$se,MAF=a$maf,N=a$samplesize)
      dfall=merge(df1,df2_window,by="SNP")

      required_cols <- c("SNP", "A1.x", "A2.x", "A1.y", "A2.y",
                         "beta.x", "beta.y", "se.x", "se.y", "MAF")
      dfall <- dfall[complete.cases(dfall[, required_cols, drop = FALSE]), , drop = FALSE]
      dfall <- dfall[
        nzchar(as.character(dfall$SNP)) &
          nzchar(as.character(dfall$A1.x)) & nzchar(as.character(dfall$A2.x)) &
          nzchar(as.character(dfall$A1.y)) & nzchar(as.character(dfall$A2.y)) &
          dfall$A1.x != dfall$A2.x & dfall$A1.y != dfall$A2.y &
          is.finite(dfall$beta.x) & is.finite(dfall$beta.y) &
          is.finite(dfall$se.x) & dfall$se.x > 0 &
          is.finite(dfall$se.y) & dfall$se.y > 0 &
          is.finite(dfall$MAF) & dfall$MAF > 0 & dfall$MAF <= 0.5,
        , drop = FALSE
      ]
      dfall = dfall[!duplicated(dfall$SNP),]
      dfall = dfall %>% filter((A1.x==A1.y&A2.x==A2.y)|(A1.x==A2.y&A2.x==A1.y)) 
      dfall = dfall %>% mutate(beta.y = ifelse(A1.x==A1.y,beta.y,-beta.y))
      
      
      dfall$VAR1 = dfall$se.x^2
      dfall$VAR2 =dfall$se.y^2
      dfall = dfall[is.finite(dfall$VAR1) & dfall$VAR1 > 0 &
                      is.finite(dfall$VAR2) & dfall$VAR2 > 0, , drop = FALSE]
      
      if (nrow(dfall) == 0L) {
        reason_now <- "NA"
        warning(paste(genelist1[i], reason_now), call. = FALSE)
        skipped <- rbind(skipped, data.frame(id = as.character(genelist1[i]), reason = reason_now))
        next
      }
      
      cdf1 <- list(
        beta    = dfall$beta.x,
        varbeta = dfall$VAR1,
        snp     = dfall$SNP,
        MAF     =dfall$MAF,  
        N       = 30596,
        type    = "quant"
      )
      cdf2 <- list(
        beta    = dfall$beta.y,
        varbeta = dfall$VAR2,
        snp     = dfall$SNP,
        s       = 121940/1376071, 
        N       = 1376071,
        type    = "cc"
      )
      
      cdf1$type = "quant"#quant
      cdf2$type = "cc"
      genename1 <- as.character(genelist1[i])
      eqtl_symbols <- unique(as.character(a$id))
      genename2 <- paste(eqtl_symbols[!is.na(eqtl_symbols) & nzchar(eqtl_symbols)], collapse = ";")
      mr_symbols <- as.character(epositions1$gene_symbol[
        !is.na(epositions1$gene_id) & epositions1$gene_id == genename1
      ])
      genename3 <- paste(unique(mr_symbols[!is.na(mr_symbols) & nzchar(mr_symbols)]), collapse = ";")
      
      for (p12_now in p12_values) {
        res = coloc.abf(cdf1,cdf2,p1=1e-4,p2=1e-4,p12=p12_now)
        need_result <- res$results %>%
          dplyr::arrange(desc(SNP.PP.H4))
        color_result1=res$summary
        result=rbind(result,data.frame(
          source="eQTLGen", eqtl_gene=genename2, gene=genename3, id=genename1,
          p12=p12_now, nsnps=color_result1[["nsnps"]],
          PP.H0.abf=color_result1[["PP.H0.abf"]],
          PP.H1.abf=color_result1[["PP.H1.abf"]],
          PP.H2.abf=color_result1[["PP.H2.abf"]],
          PP.H3.abf=color_result1[["PP.H3.abf"]],
          PP.H4.abf=color_result1[["PP.H4.abf"]],
          topsnp=need_result$snp[1],
          gene_chr=window$chr, gene_pos=window$center,
          window_lower=window$lower, window_upper=window$upper,
          window_bp=WINDOW_BP
        ))
      }
    }
  }
}
#}

write.table(result, OUT_FILE, sep = "\t", row.names = F, quote = F)
write.table(skipped, SKIPPED_FILE, sep = "\t", row.names = F, quote = F)
message("Results: ", OUT_FILE, "\nSkipped: ", SKIPPED_FILE)

