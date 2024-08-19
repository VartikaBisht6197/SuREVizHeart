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
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(DBI)
library(DT)
library(formattable)
library(ggplot2)
library(ggseqlogo)
library(GenomicRanges)
library(gridExtra)
library(here)
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


# Define global variables
# These variables are used throughout the application for various purposes.
# They include paths to data directories, cutoff values, and default settings.

# p-value cutoff used to define the opacity of SNPs in the SuREViz plot panel
pvalcutoff <- 0.004581764808767
# Cutoff value for SuRE signal intensity
SureSignalcutoff <- 3

# Directories for data and plots
DBQueryScriptsDIR <- file.path(appDIR, "modules" ,"query")
DBDIR <- file.path(dirname(appDIR), "data")

# Define lists of bigwig files for different patients
bigwigs <- list(
  "Patient1" = file.path(DBDIR, "data", "HLHS.X57.bw"),
  "Patient2" = file.path(DBDIR, "data", "HLHS.X59.bw"),
  "Patient3" = file.path(DBDIR, "data", "HLHS.X67.bw"),
  "Patient4" = file.path(DBDIR, "data", "HLHS.X68.bw"),
  "Patient5" = file.path(DBDIR, "data", "HLHS.X86.bw")
)

# Define list of ATAC-seq data files
AC16ATACbw <- list(
  "AC16\nATAC" = file.path(DBDIR, "data", "AC16.ATACseq.hg38.bw")
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
