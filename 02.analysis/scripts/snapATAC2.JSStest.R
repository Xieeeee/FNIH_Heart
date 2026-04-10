#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load h5ad")
parser$add_argument("-a", "--attr", required=TRUE, help="attribute's column name")
parser$add_argument("-l", "--lrTest", required=TRUE, help="LR test results from snapatac2")
parser$add_argument("--path_to_peak", default="/tscc/projects/ps-renlab2/y2xie/projects/77.LC/81.FNIH_DPT_IGM_240827/reference/FNIH_peaks/", help="path to peak set for each cluster")
parser$add_argument("-s", "--sig", default = 0.01, help="significant cutoff for LR test [default %(default)s]")
parser$add_argument("-e", "--eFDR", default = 0.05, help="empirical FDR cutoff [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
parser$add_argument("-n", "--norm", required=FALSE, help="norm_df_scale")
parser$add_argument("-c", "--scale", required=FALSE, help="prob_df_scale")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("anndata"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("philentropy"))
suppressPackageStartupMessages(library("fdrtool"))
suppressPackageStartupMessages(library("mixtools"))
library("tictoc")

RDataF = args$input
path_to_peak = args$path_to_peak
lrF = args$lrTest
attr = args$attr
sig = as.numeric(args$sig)
eFDR = as.numeric(args$eFDR)
outF = args$output
if(!is.null(args$norm)){
	prop.norm.df = read.csv(args$norm, row.names = 1)	
	prop.scale.df = read.csv(args$scale, row.names = 1)
}
#---------
# init
x.sp <- read_h5ad(RDataF)
lr.df <- read.csv(lrF)

#-----------------------
# Jensen-Shannon divergence
#----------------
# filter used peaks
# type.lst <- unique(x.sp$obs[,attr])
# if(length(type.lst) == 1){
#   inf = list.files(path = path_to_peak, pattern = paste("*", type.lst, "*.bed",sep=""), full.names=T )
#   peak.lst <- read.table(inf, sep="\t", header=F) %>% mutate(range = paste0(V1, ":", V2, "-", V3))
#   peak.lst <- unique(peak.lst$range)
# }else{
#   peak.lst <- list()
#   for(i in type.lst){
#     inf = list.files(path = path_to_peak, pattern = paste("*", i, "*.bed",sep=""), full.names=T)
#     inp = read.table(inf, sep="\t", header=F);
#     peak.lst[[i]] <- inp %>% mutate(range = paste0(V1, ":", V2, "-", V3))
#   }
#   peak.lst <- do.call(rbind, peak.lst)
#   peak.lst <- unique(peak.lst$range)
# }

# 1. proprotion for each cluster
attr.lst <- unique(x.sp$obs[,attr])

tic("get/norm prop")
if(is.null(args$norm)){
prop.lst <- list()
medacc.lst <- list()

for(t in attr.lst){
  idx <- which(x.sp$obs[, attr] == t)
  prop <- Matrix::colSums(x.sp$X[idx, ])
  medacc <- Matrix::rowSums(x.sp$X[idx, ])
  medacc <- median(medacc)
  prop <- prop / sum(prop) * 10^4
  prop.lst[[t]] <- prop
  medacc.lst[[t]] <- medacc
}

prop.df <- do.call(cbind, prop.lst)
medacc.df <- do.call(cbind, medacc.lst)
scalef.df <- mean(log10(medacc.df[1,])) / log10(medacc.df)

# 3. correct/norm and scale
prop.norm.df <- sweep(prop.df, 2, scalef.df[1,], "*")
prop.scale.df <- t(apply(prop.norm.df, 1, function(x){x/sum(x)}))
prop.scale.df[is.na(prop.scale.df)] <- 0
}
toc()

# 4. Jensen–Shannon divergence & Speciﬁcity
calJSD <- function(p, q){
    P <- p
    Q <- q
    x <- rbind(P,Q)
    suppressMessages(JSD(x))
}

tic("JS calculation")

jss.df <- data.frame(peaks=as.character(),
                     JSD=as.numeric(),
                     JSS=as.numeric(),
                     z=as.numeric(),
                     normCPM = as.numeric(),
                     JSscore=as.numeric(),
                     attr=as.character(),
                     group=as.character(),
                     zPval=as.numeric(),
                     zQval=as.numeric(),
                     zLocFDR=as.numeric(),
                     zCUT=as.numeric(),
                     JSSeFDRcut=as.numeric(),
                     passJSS=as.character())

for(i in 1:length(attr.lst)){
  attr.sel <- attr.lst[i]
  diffgene <- lr.df %>% filter(celltype == attr.sel & adjusted.p.value < 0.05) %>% select(feature.name) %>% unlist
  Q <- rep(0, length(attr.lst))
  Q[i] <- 1
  jsd <- apply(prop.scale.df, 1, calJSD, q=Q)
  jss <- 1 - sqrt(jsd)

  out <- data.frame(peaks = rownames(prop.scale.df)) 
  out$JSD <- jsd
  out$JSS <- jss
  out$z <- ( jss - mean(jss) ) / sd(jss)
  out$normCPM <- prop.norm.df[,i]
  out$JSscore <- jss^2 * prop.norm.df[,i]
  out$attr <- attr.sel
  out$group <- "bgd"
  idx <- which(gsub(":", "-", out$peaks) %in% gsub(":", "-", diffgene))
  out[idx, "group"] <- "diff"

  # fit
  fdrOUT <- fdrtool(out$z, plot=FALSE)
  zCUT <- fdrOUT$param[,"cutoff"]
  out$zPval <- fdrOUT$pval
  out$zQval <- fdrOUT$qval
  out$zLocFDR <- fdrOUT$lfdr
  out$zCUT<- zCUT
  JSSeFDRcut <- quantile(out[which(out$group=="bgd"), ][["JSS"]], 1-eFDR)
  out$JSSeFDRcut <- JSSeFDRcut
  out$passJSS <- "no"
  out[which(out$JSS > JSSeFDRcut), "passJSS"] <- "yes"
  jss.df <- rbind(jss.df, out)
}

toc()

# 5. output 
fwrite(jss.df, file=paste(outF, "JSStest.out.tsv", sep="."), quote =F, col.names = T, row.names = F, sep="\t")

pdfname <- paste(outF, "JSStest.fit.pdf", sep=".")
pdf(pdfname)

for(i in 1:length(attr.lst)){
  attr.sel <- attr.lst[i]
  dat <- jss.df[which(jss.df$attr == attr.sel), ]
  zCUT <- unique(dat$zCUT)
  JSSeFDRcut <- unique(dat$JSSeFDRcut)

p1 <- ggplot(dat, aes(x=JSS)) +
  geom_density(aes(fill=group), alpha=0.5) +
  theme_bw() + ggtitle(paste(attr.sel, " (JSS)", sep="")) +
  geom_vline(xintercept = JSSeFDRcut, color="black") +
  theme(axis.text.x=element_text(colour = "black", size = 8),
        axis.text.y=element_text(colour = "black", size = 8),
        axis.title.x=element_text(colour = "black", size = 8),
        axis.title.y=element_text(colour = "black", size = 8))
  
p2 <- ggplot(dat, aes(x=z)) +
  geom_density(fill="grey", alpha=0.5) +
  theme_bw() + ggtitle(paste(attr.sel, " (z-score of JSS)", sep="")) +
  geom_vline(xintercept = zCUT) +
  theme(axis.text.x=element_text(colour = "black", size = 8),
        axis.text.y=element_text(colour = "black", size = 8),
        axis.title.x=element_text(colour = "black", size = 8),
        axis.title.y=element_text(colour = "black", size = 8))
print(p1)
print(p2)
}
dev.off()




