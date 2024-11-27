#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script, named global.R, is designed to load necessary libraries, define global variables,
# and create utility functions that are used throughout the Shiny application.
# It prepares the environment with essential packages and configurations needed for data handling,
# visualization, and user interface elements in the application.
#-----------------------------------------------------------------------------------------------

# Load required libraries
# These libraries provide various functionalities needed for data manipulation, visualization,
# and user interface components in the Shiny app.
suppressPackageStartupMessages({
  library(bslib)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
  library(DBI)
  library(DT)
  library(formattable)
  library(ggplot2)
  library(ggseqlogo)
  library(GenomicRanges)
  library(gridExtra)
  library(kableExtra)
  library(memes)
  library(optparse)
  library(patchwork)
  library(pheatmap)
  library(plotly)
  library(rtracklayer)
  library(RColorBrewer)
  library(RSQLite)
  library(Signac)
  library(shiny)
  library(shinyFiles)
  library(shinyWidgets)
  library(shinyalert)
  library(shinycssloaders)
  library(shinydashboard)
  library(shinyjs)
  library(tidyr)
  library(universalmotif)
  library(uuid)
  library(zip)
})

# Define global variables
# These variables are used throughout the application for various purposes.
# They include paths to data directories, cutoff values, and default settings.

# p-value cutoff used to define the opacity of SNPs in the SuREViz plot panel
pvalcutoff <- 0.004581764808767
# Cutoff value for SuRE signal intensity
SureSignalcutoff <- 3
# Set alpha ranges
# ifelse(suredata$wilcox.p.value < pvalcutoff,
#    -log10(suredata$wilcox.p.value) * 100,
#    -log10(suredata$wilcox.p.value) / 100
#  ) -log10(0.004581764808767) / 100
alpha_max <- -log10(1.175478e-40) * 100 # Minimum p value
alpha_min <- -log10(0.9998009) / 100 # Second max p value (as max will be 0, as max p value is 1)


# Directories for data and plots
DBQueryScriptsDIR <- file.path(appDIR, "modules" ,"query")

# Define lists of bigwig files for different patients
bigwigs <- list(
  "SuREX38" = file.path(dataDIR, "data", "SuREX38.bw"),
  "SuREX57" = file.path(dataDIR, "data", "SuREX57.bw"),
  "SuREX59" = file.path(dataDIR, "data", "SuREX59.bw"),
  "SuREX67" = file.path(dataDIR, "data", "SuREX67.bw"),
  "SuREX68" = file.path(dataDIR, "data", "SuREX68.bw"),
  "SuREX86" = file.path(dataDIR, "data", "SuREX86.bw")
)

# Define list of ATAC-seq data files
AC16ATACbw <- list(
  "AC16\nATAC" = file.path(dataDIR, "data", "AC16.ATACseq.hg38.bw")
)

# Define hg38 PhyloP 30 way conservation track
Consbw <- list(
  "Conservation\nphyloP 30 way" = file.path(dataDIR, "data", "hg38.phyloP30way.bw")
)

# Utility functions

# Function to convert values to kilo bases (kb)
to_kb <- function(x) {
  paste0(x / 1000, "kb")
}

# Function to convert 'kb' values back to numeric base pairs
to_number <- function(x) {
  num_str <- gsub("kb", "", x, ignore.case = TRUE)
  num <- as.numeric(num_str)
  num * 1000
}

# Define custom ticks and labels for plotting
custom_ticks <- c(1000, 5000, 10000, 20000, 50000, 100000)
custom_tick_labels <- to_kb(custom_ticks)

# Parameters for motif visualization
motifvh <- 1000 # Height of motif plots in pixels
motifs <- 1 # Number of motifs to be considered
