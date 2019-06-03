#=============================================================================#
# eQTL_dataProcessing_and_Analysis.R                                          #
#                                                                             #
# Version: 1.0                                                                #
# Date: Jan 21, 2019                                                          #
# Author: Michiel Adriaens, PhD; MaCSBio, Maastricht University               #
# History:                                                                    #
#  1.0: Creation                                                              #
#                                                                             #
#=============================================================================#

#-----------------------------------------------------------------------------#
# 1. Initialization
#-----------------------------------------------------------------------------#

# Load libraries
require(readxl)
require(edgeR)
require(biomaRt)
require(pcaMethods)
require(WGCNA)
options(stringsAsFactors = F)

# Set directories
DATA.DIR1 <- "~/macsbio_research/michiel.adriaens/Data/Diogenes/RNA-seq"
DATA.DIR2 <- "~/macsbio_research/michiel.adriaens/Data/Johan/VitD/Data"
RES.DIR <- "~/macsbio_research/michiel.adriaens/Data/Johan/VitD/Results"

# Load GX data
setwd(DATA.DIR1)
.t <- read.delim("rawcnts_CID1.txt", as.is = T, row.names = 1)
gxData.original <- .t
.t2 <- data.frame(read_excel("matching_RNAseqID_matches_overview.xlsx", 
        sheet = "missing_RNAseqID_matches_overvi"))
.t2[.t2[, "hasNoMatch"] == "TRUE", 2] <- .t2[.t2[, "hasNoMatch"] == "TRUE", 5]
gxSampleTable <- .t2[, 1:2]; rm(.t2)
colnames(.t) <- gxSampleTable[match(colnames(.t),
        gxSampleTable[, 1]), 2]
gxData <- .t[, !is.na(colnames(.t))]; rm (.t)

# Normalize gene expression data using scaling factors (edgeR)
expr <- DGEList(counts = as.matrix(gxData))
expr <- calcNormFactors(expr)
gxData.norm <- cpm(expr, normalized.lib.sizes = TRUE, log = TRUE)

# Load SNP data (coordinates are in in HG19)
setwd(DATA.DIR2)
clinData <- data.frame(read_excel("Dataset VDR Diogenes N=553.xlsx", 
        sheet = "Sheet1"))
.t <- data.frame(read_excel("Dataset VDR Diogenes N=553.xlsx", 
    sheet = "Coding for SNPs", range = "B3:E7"))
snpAnn <- cbind(.t[, 1:3], 'alt.allele' = gsub(".{1}(.)", "\\1", .t[, 4]))
rm (.t)
vdrTargetGeneData <- data.frame(read_excel("VitDtargetGenes_21012019.xlsx", 
        sheet = "Sheet1"))

# Recode SNP data: ref = 0, alt = 1
.t <- clinData[, 2:5]
colnames(.t) <- paste(colnames(.t), "alleles", sep = ".")
.t2 <- cbind(clinData, .t); rm (.t)
for (i in 1:4) {
    .t2[, i + 1] <- as.numeric(gsub("0/1", 1, gsub("1/1", 2, gsub("0/0", 0, 
        gsub(snpAnn[i, 4], 1, gsub(snpAnn[i, 3], 0, clinData[, i + 1]))))))
}
colnames(.t2)[1] <- "ID"
.t2[, 1] <- gsub("-", "_", .t2[, 1])
snpData <- .t2[match(colnames(gxData.norm), .t2[, 1]), ]; rm (.t2)
# FIXME: a lot of missing values?

# Define VitD receptor ID
vitdId <- "ENSG00000111424"

# Filter out lowly expressed genes
plot(density(as.vector(as.matrix(gxData.norm)))) 
# anything with median below -5 seems noise
gxData.norm.f <- gxData.norm[rowMedians(gxData.norm) > -5, ]

#-----------------------------------------------------------------------------#
# 2. PCA analysis of transcriptomics data
#-----------------------------------------------------------------------------#
pcaRes <- pca(t(gxData.norm.f), nPcs = 10)

# Highest ~ 4-5% (nothing much then)

#-----------------------------------------------------------------------------#
# 3. VDR cis eQTL exploration
#-----------------------------------------------------------------------------#

# Make exploratory boxplots
for (s in c("rs731236", "rs1544410", "rs7975232", "rs10735810")) { #*
    setwd(RES.DIR)
    png(paste0("VDR_eQTL_boxplot_", s, ".png"), 
        height = 1200, width = 1200, res = 188)
    boxplot(split(gxData.norm.f[vitdId, ], 
        snpData[, paste0(s, ".alleles")]),
        col = colorRampPalette(c("dodgerblue2", "gold1"))(5)[c(1,3,5)],
        border = "lightgrey", frame = F, outline = F, notch = T,
        main = expression(italic("VDR")),
        xlab = s, ylab = "Gene expression (log2)",
        ylim = c(-1, 3.5))
    dev.off()
}

#* "rs10735810" was merged into rs2228570
# Nothing significant or striking

#-----------------------------------------------------------------------------#
# 4. VDR SNPs: all cis eQTLs (+/- 1 Mb)
#-----------------------------------------------------------------------------#

# Function to get all genes in cis (HG38, for completeness sake)
snp2cisgenes <- function(snp) {
    if (!exists("grch38.snp") | !exists("grch38.gene")) {
        grch38.snp <<- useMart(
            biomart = "ENSEMBL_MART_SNP", 
            dataset = "hsapiens_snp")
        grch38.gene <<- useMart(
            biomart = "ENSEMBL_MART_ENSEMBL", 
            dataset = "hsapiens_gene_ensembl")
    }
    Sys.sleep(0.5)
    pos <- getBM(attributes = c("chr_name", "chrom_start", "chrom_end"), 
        filters = "snp_filter", 
        values = snp, 
        mart = grch38.snp)
    Sys.sleep(0.5)
    gene <- getBM(attributes = c("hgnc_symbol", "external_gene_name", 
            "ensembl_gene_id"), 
        filters = c("chromosome_name", "start", "end"), 
        values = list(pos$chr_name, pos$chrom_start - 1E6, pos$chrom_end + 1E6), 
        mart = grch38.gene)
    return(gene)
}

# Function to perform eQTL test (additive model)
eqtlTest <- function(gene, snp) {
    if (is.element(gene, rownames(gxData.norm.f))) {
    gx <- as.vector(gxData.norm.f[gene, ])
    snp <- as.vector(snpData[, snp])
    sex <- as.vector(snpData[, "sex"])
    center <- as.vector(snpData[, "center"])
    age <- as.vector(snpData[, "age"])
    summary(lm(gx ~ snp + sex + center + age))$coefficients[2, 1:4] 
    } else {
        rep(NA, 4)
    }
}

# Function to perform all cis eQTL tests for a given SNP
eqtlCisTest <- function(snp) {
    if (snp == "rs10735810") {
        genes <- snp2cisgenes("rs2228570")
        # rs10735810 was merged into rs2228570
    } else {
        genes <- snp2cisgenes(snp)
    }
    genes <- genes[!is.na(genes[, 3]), ]
    resTable <- c()
    for (i in seq(along = genes[, 3])) {
        resTable <- rbind(resTable, eqtlTest(genes[i, 3], snp)) 
    }
    cbind(genes, resTable)
}

# Run eQTL analysis for all cis genes
for (s in c("rs731236", "rs1544410", "rs7975232", "rs10735810")){ 
    eqtlRes <- try(eqtlCisTest(s))
    print(eqtlRes[which(eqtlRes[, 7] < 0.05), 3])
    try({
    for (g in eqtlRes[which(eqtlRes[, 7] < 0.05), 3]) {
        # Make eQTL boxplots for significant eQTLs
        setwd(RES.DIR)
        gn <- eqtlRes[match(g, eqtlRes$ensembl_gene_id), "external_gene_name"]
        p <- eqtlRes[match(g, eqtlRes$ensembl_gene_id), 7]
        png(paste0(g, "_eQTL_boxplot_", s, ".png"), 
            height = 1200, width = 1200, res = 188)
        boxplot(split(gxData.norm.f[g, ], 
                snpData[, paste0(s, ".alleles")]),
            col = colorRampPalette(c("dodgerblue2", "gold1"))(5)[c(1,3,5)],
            border = "lightgrey", frame = F, outline = F, notch = T,
            main = substitute(expression(paste(italic(gn))),
                env = list(gn = gn)),
            xlab = paste(s, "\n(p = ", signif(p, 2), ")"), 
            ylab = "Gene expression (log2)")
        dev.off()
    }})
}

#-----------------------------------------------------------------------------#
# 5. VDR target genes: co-expression network analysis
#-----------------------------------------------------------------------------#

# First assess the soft thresholding power (to force scale free topology)
# (adapted from WGCNA tutorial)
#------------------------------------------------------------------------
wgcnaInput <- t(gxData.norm.f[is.element(rownames(gxData.norm.f),
            vdrTargetGeneData$Ensembl), ])
powers <- c(1:20)
sft <- pickSoftThreshold(wgcnaInput, powerVector = powers, verbose = 5,
		RsquaredCut = 0.8)

# Plot the results:
setwd(RES.DIR)
setwd("VDR target genes")
png("WGCNA_threshold_assessment.png", width = 1200, height = 600)
par(mfrow = c(1,2))
cex1 <- 0.9
# Left plot: Scale-free topology fit index as a function of the soft-thresholding power:
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
		xlab = "Soft Threshold (power)",
		ylab = "Scale Free Topology Model Fit,signed R^2",
		type = "n",
		main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[, 2],
		labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "red") # R-squared cut-off of 0.8

# Right plot: Mean connectivity as a function of the soft-thresholding power:
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
		xlab = "Soft Threshold (power)",
		ylab = "Mean Connectivity", 
		type = "n",
		main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
		labels = powers, cex = cex1, col = "red")
dev.off()
# Conclusion: which power to use? 	
power <- 5

# Run WGCNA
#-----------
wgcnaRes <- blockwiseModules(wgcnaInput, power = power,
		TOMType = "unsigned", minModuleSize = 10, verbose = 3)

# Summarize modules:
table(wgcnaRes$colors)

# REMARK: grey module is the "dumpster"

# Visualize WGCNA results
#------------------------
# Plot the dendrogram and the module colors underneath
png("WGCNA_res.png", height = 1200, width = 1200)
plotDendroAndColors(wgcnaRes$dendrograms[[1]], 
		wgcnaRes$colors[wgcnaRes$blockGenes[[1]]],
		"Module colors", dendroLabels = FALSE, hang = 0.03, 
		addGuide = TRUE, guideHang = 0.05)
dev.off()

# REMARK: hierarchical average linkage clustering: the distance between clusters
#         (reflected in the height) is defined as the average distance between 
#         all inter-cluster pairs;
# 		  distance here is defined as the Pearson correlation statistic to
#         the power of 5

# Plot invididual modules
# FIXME: since genes can have opposite directions, isn't the addition of
#        'colMeans' incorrect? (used to be matplot, genes separately)
#        there needs to be some sorting and ordering
plotModule <- function(moduleId) {
	exprData <- t(wgcnaInput) - rowMeans(t(wgcnaInput))
	if (is.character(moduleId)) {
		mod <- moduleId
	} else {
		mod <- names(table(wgcnaRes$colors))[moduleId]
	}
	x <- colMeans(exprData[wgcnaRes$color == mod, ])
	o <- order(x)
	plot(1:length(x), x[o], type = "b", xlab = "Individual", 
        ylab = "Scaled GX value",
			main = paste("GX of individual genes\n (module ", 
					mod, ", n = ", table(wgcnaRes$colors)[[mod]],
					")", sep = ""), lwd = 2,
			pty = 1, lty = 3, pch = 1,
			col = colorRampPalette(c("lightgrey", mod))(9)[5])
	
}
# 3 modules in total (excluding the grey module)
for (i in unique((wgcnaRes$colors))) {
	png(paste("WGCNA_module", i, ".png", sep = ""), height = 600, width = 600)
	plotModule(i)
	dev.off()
}	

# Test for association
eqtlMEtest <- function(color, snp) {
    gx <- wgcnaRes$MEs[, paste0("ME", color)]
    message(colnames(wgcnaRes$MEs)[i])
    snp <- as.vector(snpData[, snp])
    sex <- as.vector(snpData[, "sex"])
    center <- as.vector(snpData[, "center"])
    age <- as.vector(snpData[, "age"])
    summary(lm(gx ~ snp + sex + center + age))$coefficients[2, 1:4] 
}
eqtlMEtest("blue", "rs731236")
eqtlMEtest("blue", "rs1544410")
eqtlMEtest("blue", "rs7975232")
eqtlMEtest("blue", "rs10735810")

eqtlMEtest("brown", "rs731236")
eqtlMEtest("brown", "rs1544410")
eqtlMEtest("brown", "rs7975232")
eqtlMEtest("brown", "rs10735810")

eqtlMEtest("turquoise", "rs731236")
eqtlMEtest("turquoise", "rs1544410")
eqtlMEtest("turquoise", "rs7975232")
eqtlMEtest("turquoise", "rs10735810")

# Nothing even approaching significance

# Test all target genes individually
for (s in c("rs731236", "rs1544410", "rs7975232", "rs10735810")){
    for (g in colnames(wgcnaInput)) {
      eqtlRes <- eqtlTest(g, s)
      if(eqtlRes[4] < 0.05) {
        # Make eQTL boxplots for significant eQTLs
        gn <- vdrTargetGeneData[match(g, vdrTargetGeneData[, 2]), 3]
        p <- eqtlRes[4]
        png(paste0(g, "_eQTL_boxplot_", s, ".png"), 
            height = 1200, width = 1200, res = 188)
        boxplot(split(gxData.norm.f[g, ], 
                snpData[, paste0(s, ".alleles")]),
            col = colorRampPalette(c("dodgerblue2", "gold1"))(5)[c(1,3,5)],
            border = "lightgrey", frame = F, outline = F, notch = T,
            main = substitute(expression(paste(italic(gn))),
                env = list(gn = gn)),
            xlab = paste(s, "\n(p = ", signif(p, 2), ")"), 
            ylab = "Gene expression (log2)")
        dev.off()
      }
    }
}
# Only two significant, would be lost after correction for mt

#-----------------------------------------------------------------------------#
# 6. Co-expression network analysis on all (expressed) genes
#-----------------------------------------------------------------------------#

# Rationale: SNP function/effect unknown: may have downstream impact 
# First assess the soft thresholding power (to force scale free topology)
# (adapted from WGCNA tutorial)
#------------------------------------------------------------------------
wgcnaInput.full <- t(gxData.norm.f)
powers <- c(1:20)
enableWGCNAThreads(12)
sft.full <- pickSoftThreshold(wgcnaInput.full, powerVector = powers, verbose = 5,
		RsquaredCut = 0.8)

# Plot the results:
setwd(RES.DIR)
setwd("Co-expression results all genes")
png("WGCNA_threshold_assessment.png", width = 1200, height = 600)
par(mfrow = c(1,2))
cex1 <- 0.9
# Left plot: Scale-free topology fit index as a function of the soft-thresholding power:
plot(sft.full$fitIndices[, 1], -sign(sft.full$fitIndices[, 3]) * sft.full$fitIndices[, 2],
		xlab = "Soft Threshold (power)",
		ylab = "Scale Free Topology Model Fit,signed R^2",
		type = "n",
		main = paste("Scale independence"))
text(sft.full$fitIndices[,1], -sign(sft.full$fitIndices[,3]) * sft.full$fitIndices[, 2],
		labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "red") # R-squared cut-off of 0.8

# Right plot: Mean connectivity as a function of the soft-thresholding power:
plot(sft.full$fitIndices[, 1], sft.full$fitIndices[, 5],
		xlab = "Soft Threshold (power)",
		ylab = "Mean Connectivity", 
		type = "n",
		main = paste("Mean connectivity"))
text(sft.full$fitIndices[, 1], sft.full$fitIndices[, 5], 
		labels = powers, cex = cex1, col = "red")
dev.off()
# Conclusion: which power to use? 	
power.full <- 15

# Run WGCNA
#-----------
wgcnaRes.full <- blockwiseModules(wgcnaInput.full, power = power.full,
		TOMType = "unsigned", minModuleSize = 10, verbose = 3,
        nThreads = 12)

# Summarize modules:
table(wgcnaRes.full$colors)

# REMARK: grey module is the "dumpster"

# Visualize WGCNA results
#------------------------
# Plot the dendrogram and the module colors underneath
png("WGCNA_res.png", height = 1200, width = 1200)
plotDendroAndColors(wgcnaRes.full$dendrograms[[1]], 
		wgcnaRes.full$colors[wgcnaRes.full$blockGenes[[1]]],
		"Module colors", dendroLabels = FALSE, hang = 0.03, 
		addGuide = TRUE, guideHang = 0.05)
dev.off()

# REMARK: hierarchical average linkage clustering: the distance between clusters
#         (reflected in the height) is defined as the average distance between 
#         all inter-cluster pairs;
# 		  distance here is defined as the Pearson correlation statistic to
#         the power of 15

# Plot all modules
plotModule.full <- function(moduleId) {
	exprData <- t(wgcnaInput.full) - rowMeans(t(wgcnaInput.full))
	if (is.character(moduleId)) {
		mod <- moduleId
	} else {
		mod <- names(table(wgcnaRes.full$colors))[moduleId]
	}
	x <- colMeans(exprData[wgcnaRes.full$color == mod, ])
	o <- order(x)
	plot(1:length(x), x[o], type = "b", xlab = "Individual", 
        ylab = "Scaled GX value",
			main = paste("GX of individual genes\n (module ", 
					mod, ", n = ", table(wgcnaRes.full$colors)[[mod]],
					")", sep = ""), lwd = 2,
			pty = 1, lty = 3, pch = 1,
			col = colorRampPalette(c("lightgrey", mod))(9)[5])
	
}
for (i in unique(wgcnaRes.full$colors)) {
	png(paste("WGCNA_module", i, ".png", sep = ""), height = 600, width = 600)
	plotModule.full(i)
	dev.off()
}	

# Test for association
eqtlMEtest.full <- function(color, snp) {
    gx <- wgcnaRes.full$MEs[, paste0("ME", color)]
    message(colnames(wgcnaRes.full$MEs)[i])
    snp <- as.vector(snpData[, snp])
    sex <- as.vector(snpData[, "sex"])
    center <- as.vector(snpData[, "center"])
    age <- as.vector(snpData[, "age"])
    summary(lm(gx ~ snp + sex + center + age))$coefficients[2, 4] 
}
for (color in unique(wgcnaRes.full$colors)) {
    message(paste("\n------------------------------------", "color:", color))
    for (s in c("rs731236", "rs1544410", "rs7975232", "rs10735810")) {
        print(eqtlMEtest.full(color, s)[1])
    }
    invisible(readline(prompt="Press [enter] to continue"))
}

# Nothing even approaching significance. Lowest is 0.06 (on 0.9 and on 0.10), 
#  but considering number of tests, this is flukey at best.


# END
