# ================================================================
# Aysha_Class_3B_Preprocessing_Limma_Final.R
# AI and Omics Research Internship (2025)
# Module II - Microarray Preprocessing, Normalization & DE Analysis
# ================================================================

# -----------------------------
# 0. Clear memory and set directory
# -----------------------------
rm(list=ls())
gc()
setwd("C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_II")
print(getwd())

# -----------------------------
# 1. Load libraries
# -----------------------------
library(limma)
library(pheatmap)
library(ggplot2)

# -----------------------------
# 2. Read raw microarray data
# -----------------------------
raw_path <- "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_II/raw_data/"
stopifnot(dir.exists(raw_path))

raw_files <- list.files(raw_path, pattern="\\.gpr$", full.names=TRUE)
stopifnot(length(raw_files) > 0)
cat("Found", length(raw_files), "GPR files.\n")

columns_list <- list(G="F635 Median", Gb="B635 Median")
targets <- data.frame(FileName=raw_files)

data_raw <- read.maimages(files=targets$FileName, source="genepix", columns=columns_list)

# Assign unique probe names
probe_names <- make.unique(as.character(data_raw$genes$Name))
rownames(data_raw$E) <- probe_names

# -----------------------------
# 3. Background correction
# -----------------------------
data_bg <- backgroundCorrect(data_raw, method="normexp", offset=50)
exprs <- data_bg$E

# Log2 transform
log_expr <- log2(exprs)

# Quantile normalization
norm_expr <- normalizeBetweenArrays(log_expr, method="quantile")

# -----------------------------
# 4. Filter probes
# -----------------------------
valid_probes <- apply(norm_expr, 1, function(x) all(is.finite(x)) && !any(is.nan(x)) && !any(is.na(x)))
filtered_expr <- norm_expr[valid_probes, ]
cat("Probes remaining after filtering for finite and NA/NaN:", nrow(filtered_expr), "\n")

# -----------------------------
# 5. Heatmap of top variable probes
# -----------------------------
var_probes <- apply(filtered_expr, 1, var, na.rm=TRUE)
valid_var <- var_probes[!is.na(var_probes) & var_probes > 0]

if(length(valid_var) == 0) stop("No probes with positive variance remain.")

top_n <- min(50, length(valid_var))
top_probes <- names(sort(valid_var, decreasing=TRUE))[1:top_n]
heatmap_data <- filtered_expr[top_probes, , drop=FALSE]

# Save Heatmap
if(!dir.exists("Results")) dir.create("Results")
png("Results/TopVariableProbes_Heatmap.png", width=1200, height=1000)
pheatmap(as.matrix(heatmap_data),
         scale="row",
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         main = paste("Top", nrow(heatmap_data), "Variable Probes"))
dev.off()
cat("Heatmap saved successfully.\n")

# -----------------------------
# 6. Differential Expression Analysis
# -----------------------------
# Define sample groups (replace with actual annotation)
groups <- factor(c(rep("normal", 19), rep("cancer", 19)))

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(filtered_expr, design)
contrast.matrix <- makeContrasts(cancer_vs_normal = cancer - normal, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, coef="cancer_vs_normal", number=Inf, adjust.method="BH")

# Classify DEGs
deg_results$threshold <- factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "Not Significant")
))

cat("Upregulated genes:", sum(deg_results$threshold == "Upregulated"), "\n")
cat("Downregulated genes:", sum(deg_results$threshold == "Downregulated"), "\n")

# -----------------------------
# 7. Save results
# -----------------------------
write.csv(filtered_expr, file="Results/Processed_Expression_Data.csv", row.names=TRUE)
write.csv(deg_results, file="Results/DEG_Results.csv", row.names=TRUE)
cat("Processed data and DEG results saved to 'Results/' folder.\n")

# -----------------------------
# 8. MA Plot (first two samples)
# -----------------------------
plotMAcustom <- function(exprs, sample1=1, sample2=2, file_name="Results/MA_Plot_Sample1_vs_2.png"){
  A <- (exprs[, sample1] + exprs[, sample2]) / 2
  M <- exprs[, sample1] - exprs[, sample2]
  png(file_name, width=1200, height=1000)
  plot(A, M, pch=20, cex=0.3, col="grey",
       main=paste("MA Plot: Sample", sample1, "vs", sample2),
       xlab="A = Average log2 expression",
       ylab="M = log2 fold change")
  abline(h=0, col="red", lty=2)
  dev.off()
  cat("MA Plot saved as", file_name, "\n")
}
plotMAcustom(filtered_expr)

# -----------------------------
# 9. DEG Summary Bar Plot
# -----------------------------
deg_counts <- table(deg_results$threshold)
deg_df <- as.data.frame(deg_counts)
colnames(deg_df) <- c("Category", "Count")

png("Results/DEG_Summary_BarPlot.png", width=1200, height=1000)
ggplot(deg_df, aes(x=Category, y=Count, fill=Category)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  ggtitle("DEG Summary") +
  ylab("Number of Genes") +
  xlab("") +
  scale_fill_manual(values=c("Upregulated"="red", "Downregulated"="blue", "Not Significant"="grey")) +
  geom_text(aes(label=Count), vjust=-0.5)
dev.off()
cat("DEG Summary Bar Plot saved.\n")

# -----------------------------
# 10. Volcano Plot
# -----------------------------
volcano <- function(deg_df, logFC_thresh=1, adjP_thresh=0.05, file_name="Results/Volcano_Plot.png"){
  deg_df$Significance <- "Not Significant"
  deg_df$Significance[deg_df$adj.P.Val < adjP_thresh & deg_df$logFC > logFC_thresh] <- "Upregulated"
  deg_df$Significance[deg_df$adj.P.Val < adjP_thresh & deg_df$logFC < -logFC_thresh] <- "Downregulated"
  
  p <- ggplot(deg_df, aes(x=logFC, y=-log10(adj.P.Val), color=Significance)) +
    geom_point(alpha=0.6, size=1.5) +
    scale_color_manual(values=c("Upregulated"="red","Downregulated"="blue","Not Significant"="grey")) +
    theme_minimal() +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value") +
    ggtitle("Volcano Plot of DEGs") +
    geom_hline(yintercept=-log10(adjP_thresh), linetype="dashed", color="black") +
    geom_vline(xintercept=c(-logFC_thresh, logFC_thresh), linetype="dashed", color="black")
  
  ggsave(file_name, p, width=8, height=6)
  cat("Volcano Plot saved as", file_name, "\n")
}
volcano(deg_results)

# -----------------------------
# Analysis Completed
# -----------------------------
cat("All analysis finished. Check the 'Results/' folder for CSVs and PNG plots.\n")
