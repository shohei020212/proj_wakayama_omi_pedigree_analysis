#!/usr/bin/env Rscript

# ==============================================================================
# Sequoia Pedigree Reconstruction Script
# ==============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Error: No VCF file provided. Usage: Rscript sequoia_analysis.R <vcf_file> [life_history_file]")
}

vcf_file <- args[1]
# Check if a second argument exists and is not the placeholder "NO_FILE"
lh_file  <- if (length(args) > 1 && args[2] != "NO_FILE") args[2] else NULL

# Load required libraries
library(sequoia)
library(vcfR)
library(adegenet)

message("--- Starting Sequoia Analysis ---")
message(paste("Input VCF:", vcf_file))

# ------------------------------------------------------------------------------
# 1. Convert VCF to Sequoia Input Format
# ------------------------------------------------------------------------------
message("Converting VCF to genotype matrix...")

# GenoConvert handles VCF import. It requires the 'vcfR' package for VCFs.
# Note: Ensure your VCF is filtered (LD pruned) before this step for performance.
GenoM <- GenoConvert(InFile = vcf_file, InFormat = "vcf")
GenoM.checked <- CheckGeno(GenoM, Return="GenoM")

# ------------------------------------------------------------------------------
# 2. Load Life History Data (Optional but Recommended)
# ------------------------------------------------------------------------------
LH_data <- NULL

if (!is.null(lh_file)) {
  message(paste("Loading Life History data from:", lh_file))
  # Assuming standard format: ID, Sex, BirthYear
  LH_data <- read.table(lh_file, header = TRUE, stringsAsFactors = FALSE)
} else {
  message("No Life History data provided. Running with genotype data only.")
}

# ------------------------------------------------------------------------------
# 3. Run Pedigree Reconstruction
# ------------------------------------------------------------------------------
message("Running pedigree reconstruction...")

# 'Module = "ped"' performs the full pedigree reconstruction
ParOUT <- sequoia(GenoM = GenoM.checked, 
                  LifeHistData = LH_data, 
                  Module = "ped",           
                  Plot = FALSE) # Turn off auto-plotting to save manually later

# ------------------------------------------------------------------------------
# 4. Save Results
# ------------------------------------------------------------------------------
message("Saving results...")

# Save the inferred pedigree
write.table(ParOUT$Pedigree, file = "Sequoia_Pedigree.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

# Save a summary of the inference
sink("Sequoia_Summary.txt")
SummarySeq(ParOUT)
sink()

# Save diagnostic plots
pdf("Sequoia_Diagnostic_Plots.pdf")
PlotAgePrior(ParOUT$AgePriors)
# You can add: PlotRelPairs(ParOUT)
dev.off()

message("--- Analysis Completed Successfully ---")