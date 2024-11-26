#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script, named jaspar.analysis.R, is used to analyze and visualize the impact of genomic
# variants on transcription factor (TF) binding motifs using data from the JASPAR database.
# The script generates plots comparing reference and alternate motifs and creates a table
# summarizing the effects of these variants on TF binding. It uses the 'ggplot2' and 'patchwork'
# packages for visualization and assumes the presence of necessary data files and motif matrices.
#
#-----------------------------------------------------------------------------------------------

# Constants for plot layout
plots_in_one_set <- 3
plots_cols <- 2

# Function to perform JASPAR analysis
jaspar.analysis <- function(jaspar.data, JASPAR.PPM) {
  # Extract transcription factors from the motif_id column
  TFs <- unlist(strsplit(jaspar.data$motif_id, ","))

  # Retrieve the sequence around the variant position
  seq <- DNAStringSet(getSeq(
    BSgenome.Hsapiens.UCSC.hg38, jaspar.data$CHROM,
    (as.numeric(jaspar.data$POS) - 15),
    (as.numeric(jaspar.data$POS) + nchar(jaspar.data$REF) + 15 - 1)
  ))
  names(seq) <- jaspar.data$POS

  # Filter JASPAR PPMs based on TFs present in the data
  JASPAR.PPM <- JASPAR.PPM[which(unlist(lapply(JASPAR.PPM, function(x) x@name %in% TFs)))]
  names(JASPAR.PPM) <- unlist(lapply(JASPAR.PPM, function(x) x@name))

  # Prepare the FIMO data frame by splitting rows and selecting the highest impact per TF
  FIMO <- data.frame(separate_rows(jaspar.data, start, end, strand, pos, motif_id, motif_alt_id,
    refs.score, alts.score, refs.pval, alts.pval, absdiff,
    max.score, effect.JASPAR, tax_group, tf_family, tf_class,
    data_type, uniprot_ids, pubmed_ids, Gene_Name,
    sep = ","
  ))
  FIMO$refs.score <- round(as.numeric(FIMO$refs.score), 0)
  FIMO$alts.score <- round(as.numeric(FIMO$alts.score), 0)
  FIMO$effect.JASPAR <- round(as.numeric(FIMO$effect.JASPAR), 3)
  # Only keep the highest impact alignment per TF
  FIMO <- rbindlist(lapply(split(FIMO, FIMO$motif_alt_id), function(x) x[which(abs(x$effect.JASPAR) == max(abs(x$effect.JASPAR))), ]))

  # Create SNP identifier
  snp <- paste(jaspar.data$CHROM, jaspar.data$POS, jaspar.data$REF, jaspar.data$ALT, sep = "-")

  # Convert sequences to Position Probability Matrices (PPMs)
  ref_str <- as.character(seq[[1]])
  alt_str <- paste0(substr(seq[[1]], 1, 15), jaspar.data$ALT, substr(seq[[1]], 17 + nchar(jaspar.data$REF) - 1, 31 + nchar(jaspar.data$REF) - 1))

  seqlen <- nchar(seq) - abs(nchar(alt_str) - nchar(ref_str))

  if (nchar(ref_str) > nchar(alt_str)) {
    ref_str <- substr(ref_str, 1, nchar(alt_str))
  }

  if (nchar(alt_str) > nchar(ref_str)) {
    alt_str <- substr(alt_str, 1, nchar(ref_str))
  }

  ref_ppm <- create_motif(ref_str)@motif
  alt_ppm <- create_motif(alt_str)@motif

  # Generate plots for reference and alternate motifs
  ref_plot <- ggplot() +
    geom_rect(aes(
      xmin = 16 - 0.5, xmax = 16 + 0.5 + nchar(jaspar.data$REF) - 1,
      ymin = 3, ymax = -1
    ), fill = "black", size = 0.5, alpha = 0.09) +
    geom_logo(as.matrix(ref_ppm), method = "bits") +
    theme_logo() +
    labs(y = "REF") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )

  alt_plot <- ggplot() +
    geom_rect(aes(
      xmin = 16 - 0.5, xmax = 16 + 0.5 + nchar(jaspar.data$ALT) - 1,
      ymin = 3, ymax = -1
    ), fill = "black", size = 0.5, alpha = 0.09) +
    geom_logo(as.matrix(alt_ppm), method = "bits") +
    theme_logo() +
    labs(y = "ALT") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )

  # Create individual plots for each TF motif
  wrplt <- list()

  for (i in 1:nrow(FIMO)) {
    highlight <- NULL
    p_motif <- NULL

    motif_id <- FIMO$motif_id[i]
    motif_alt_id <- FIMO$motif_alt_id[i]
    start_pos <- as.numeric(FIMO$start[i])
    end_pos <- as.numeric(FIMO$end[i])
    matched_seq <- FIMO$matched_sequence[i]
    motif_strand <- FIMO$strand[i]
    effect_JASPAR <- FIMO$effect.JASPAR[i]

    # Get the motif and adjust for strand orientation
    motif <- JASPAR.PPM[[motif_id]]@motif
    if (motif_strand == "-") {
      motif <- motif_rc(JASPAR.PPM[[motif_id]])@motif
    } else {
      motif <- JASPAR.PPM[[motif_id]]@motif
    }

    # Create an empty PPM matrix of the same length as the sequence
    seq_length <- 31
    empty_matrix <- matrix(0, nrow = 4, ncol = seq_length)
    rownames(empty_matrix) <- c("A", "C", "G", "T")

    # Insert the motif into the correct position in the empty matrix
    motif_length <- ncol(motif)
    empty_matrix[, start_pos:(start_pos + motif_length - 1)] <- motif

    # Create plot title based on the effect of the variant
    if (effect_JASPAR > 0) {
      motif_title <- paste(snp, "disrupts TF", motif_alt_id, "(", motif_id, ")", "binding to REF", sep = " ")
    } else {
      motif_title <- paste(snp, "enhances TF", motif_alt_id, "(", motif_id, ")", "binding to ALT", sep = " ")
    }

    # Create the motif plot
    p_motif <- ggplot() +
      geom_logo(empty_matrix, method = "bits") +
      theme_logo() +
      labs(y = motif_alt_id) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
      )

    plts <- list(ref_plot + labs(title = motif_title), alt_plot, p_motif)

    # Add the plot to the list
    wrplt[[i]] <- wrap_plots(plts, ncol = 1)
  }

  # Combine all individual plots into one plot
  combined_plot <- wrap_plots(wrplt, ncol = 1)

  # Create and format the summary table
  FIMO <- FIMO[, c("motif_id", "motif_alt_id", "tf_family", "tf_class", "effect.JASPAR")]
  colnames(FIMO) <- c("JASPAR 2022 Motif ID", "JASPAR 2022 Motif Name", "TF Family", "TF Class", "REF/ALT TF Binding")

  # Return the combined plot and the summary table
  return(list(plot = combined_plot, table = FIMO))
}

# Initialize output list
jaspar_output <- list(
  "plot" = NULL,
  "table" = data.frame()
)

# Initialize jaspar.data as NULL
jaspar.data <- NULL

# Process data if position is specified and suredata is not empty
if (!is.na(pos) && nrow(suredata) != 0) {
  # Run grep command to check for rows in the JASPAR data file
  grep_output <- system(paste0("grep -w ", pos, " ", file.path(dataDIR, "JASPAR2022.HLHS.raQTL.TF.filter.txt")), intern = TRUE)

  # Check if the grep command returns any rows
  if (length(grep_output) > 0) {
    # Read the data into jaspar.data
    jaspar.data <- fread(cmd = paste0("grep -w ", pos, " ", file.path(dataDIR, "JASPAR2022.HLHS.raQTL.TF.filter.txt")), sep = "\t")
    colnames(jaspar.data) <- colnames(fread(cmd = paste("head -n1", file.path(dataDIR, "JASPAR2022.HLHS.raQTL.TF.filter.txt")), header = TRUE))

    if (nrow(jaspar.data) == 1 & (jaspar.data$CHROM == chr & jaspar.data$POS == pos & jaspar.data$REF == query_snps$REF & jaspar.data$ALT == query_snps$ALT)) {
      # Load JASPAR PPM data
      JASPAR.PPM <- read_meme(file.path(dataDIR, "data", "JASPAR2022_CORE_vertebrates_non-redundant_v2.meme"))
      names(JASPAR.PPM) <- unlist(lapply(JASPAR.PPM, function(x) x@name))

      # Run the JASPAR analysis
      jaspar_output <- jaspar.analysis(jaspar.data, JASPAR.PPM)
      motifs <- (nrow(jaspar_output$table) + 2) / plots_in_one_set
      message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "JASPAR Plot : Plotted JASPAR TFs affected by query variant. ✅"))
    }
  } else {
    # Render an empty plot if no data is found
    if (query_snps$raQTL == "raQTL") {
      jaspar_output$plot <- ggplot() +
        theme_minimal() +
        annotate("text", x = 0.5, y = 0.5, label = "The variant does not affect any TFBS\n(as defined by JASPAR 2022 CORE) significantly.", size = 6, hjust = 0.5, vjust = 0.5) +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()
        )
      message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "JASPAR Plot : Plotted empty plot. The variant does not affect any TFBS (as defined by JASPAR 2022 CORE) significantly. ⚠️"))
    } else {
      jaspar_output$plot <- ggplot() +
        theme_minimal() +
        annotate("text", x = 0.5, y = 0.5, label = "The variant selected is not an raQTL\nThis view is only possible for raQTLs.", size = 6, hjust = 0.5, vjust = 0.5) +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()
        )
      message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "JASPAR Plot : Plotted empty plot. The variant selected is not an raQTL. ⚠️"))
    }
  }
} else {
  # Render an empty plot if no position is provided
  jaspar_output$plot <- ggplot() +
    theme_minimal() +
    annotate("text", x = 0.5, y = 0.5, label = "Display is only possible for the variant view.", size = 6, hjust = 0.5, vjust = 0.5) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "JASPAR Plot : Plotted empty plot. No variant selected. ⚠️"))
}
