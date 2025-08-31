# ---------- 1. Define helper function ----------
classify_gene <- function(logFC, padj) {
  logFC <- as.numeric(logFC)
  padj  <- as.numeric(padj)
  
  if (is.na(padj)) padj <- 1  # missing padj -> treat as non-sig
  
  if (!is.na(logFC) && padj < 0.05) {
    if (logFC > 1) {
      return("Upregulated")
    } else if (logFC < -1) {
      return("Downregulated")
    } else {
      return("Not significant")
    }
  }
  return("Not significant")
}

# ---------- 2. Input files (absolute paths) ----------
input_files <- c(
  "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/raw_data/DEGs_Data_1.csv",
  "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/raw_data/DEGs_Data_2.csv"
)

# ---------- 3. Ensure results folder exists ----------
results_dir <- "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results"
if (!dir.exists(results_dir)) dir.create(results_dir)

# ---------- 4. Helper to find matching column names ----------
find_col <- function(df, candidates) {
  found <- intersect(candidates, names(df))
  if (length(found) > 0) return(found[1])
  return(NA_character_)
}

# ---------- 5. Possible column names ----------
log_choices  <- c("log2FoldChange", "logFC", "log_fold_change", "log2FC")
padj_choices <- c("padj", "p.adjust", "pvalue", "p_value", "pVal", "adj.P.Val")

# ---------- 6. Loop over each file ----------
for (infile in input_files) {
  cat("\nProcessing:", infile, "\n")
  
  if (!file.exists(infile)) {
    cat("  ❌ File not found:", infile, "\n")
    next
  }
  
  # read data
  df <- read.csv(infile, header = TRUE)
  
  # find correct columns
  log_col  <- find_col(df, log_choices)
  padj_col <- find_col(df, padj_choices)
  
  if (is.na(log_col) | is.na(padj_col)) {
    cat("  ⚠️ Could not find required columns in:", infile, "\n")
    next
  }
  
  # classify
  df$Regulation <- mapply(classify_gene, df[[log_col]], df[[padj_col]])
  
  # output file name
  outfile <- file.path(results_dir, paste0("classified_", basename(infile)))
  
  write.csv(df, outfile, row.names = FALSE)
  cat("  ✅ Saved results to:", outfile, "\n")
}

cat("\nAll done. Processed files are in the 'results/' folder.\n")
classified <- read.csv("C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/classified_DEGs_Data_1.csv")

table(classified$Regulation)
classified <- read.csv("C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/classified_DEGs_Data_2.csv")

table(classified$Regulation)
summary_counts <- table(classified$Regulation)
write.csv(as.data.frame(summary_counts), 
          "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/summary_DEGs_Data_1.csv")
summary_counts <- table(classified$Regulation)
write.csv(as.data.frame(summary_counts), 
          "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/summary_DEGs_Data_2.csv")
summary(classified$padj)
summary(classified$logFC)

# How many padj < 0.05?
sum(classified$padj < 0.05, na.rm = TRUE)

# How many abs(logFC) > 1?
sum(abs(classified$logFC) > 1, na.rm = TRUE)
# Step 1. Load classified datasets
classified1 <- read.csv("C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/classified_DEGs_Data_1.csv")
classified2 <- read.csv("C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/classified_DEGs_Data_2.csv")

# Step 2. Check classification counts
cat("\n--- Dataset 1 Regulation counts ---\n")
print(table(classified1$Regulation))

cat("\n--- Dataset 2 Regulation counts ---\n")
print(table(classified2$Regulation))

# Step 3. Explore padj and logFC
cat("\n--- Dataset 1 padj summary ---\n")
print(summary(classified1$padj))
cat("\n--- Dataset 1 logFC summary ---\n")
print(summary(classified1$logFC))

cat("\n--- Dataset 2 padj summary ---\n")
print(summary(classified2$padj))
cat("\n--- Dataset 2 logFC summary ---\n")
print(summary(classified2$logFC))

# Step 4. Count significant genes
cat("\n--- Dataset 1 Significant Counts ---\n")
cat("padj < 0.05:", sum(classified1$padj < 0.05, na.rm = TRUE), "\n")
cat("|logFC| > 1:", sum(abs(classified1$logFC) > 1, na.rm = TRUE), "\n")

cat("\n--- Dataset 2 Significant Counts ---\n")
cat("padj < 0.05:", sum(classified2$padj < 0.05, na.rm = TRUE), "\n")
cat("|logFC| > 1:", sum(abs(classified2$logFC) > 1, na.rm = TRUE), "\n")

# Step 5. (Optional) Combine both datasets
combined <- rbind(
  transform(classified1, Dataset = "Dataset1"),
  transform(classified2, Dataset = "Dataset2")
)

# Save combined file
write.csv(combined,
          "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/combined_classified_DEGs.csv",
          row.names = FALSE)

cat("\n✅ All steps done. Combined file saved.\n")
# Count regulation categories for each dataset
counts1 <- table(classified1$Regulation)
counts2 <- table(classified2$Regulation)

summary_all <- data.frame(
  Regulation = union(names(counts1), names(counts2)),
  Dataset1 = as.numeric(counts1[match(union(names(counts1), names(counts2)), names(counts1))]),
  Dataset2 = as.numeric(counts2[match(union(names(counts1), names(counts2)), names(counts2))])
)

summary_all[is.na(summary_all)] <- 0  # replace NA with 0

print(summary_all)

# Save to results folder
write.csv(summary_all,
          "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/summary_all_DEGs.csv",
          row.names = FALSE)
# Filter DEGs with padj < 0.05 and |logFC| > 1
sig_DEGs <- subset(classified, padj < 0.05 & abs(logFC) > 1)

# Save
write.csv(sig_DEGs, "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/significant_DEGs_Data_1.csv", row.names = FALSE)
# Filter DEGs with padj < 0.05 and |logFC| > 1
sig_DEGs <- subset(classified, padj < 0.05 & abs(logFC) > 1)

# Save
write.csv(sig_DEGs, "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/results/significant_DEGs_Data_2.csv", row.names = FALSE)
