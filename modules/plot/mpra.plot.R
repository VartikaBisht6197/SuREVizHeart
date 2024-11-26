#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script is designed to generate plots for MPRA (Massively Parallel Reporter Assay) data
# in the form of a signal plot, based on a user-provided TSV file. The file is expected to
# contain data related to genomic variants with the following columns:
# - chr: Chromosome (in the format "chrX" where X is a digit or Y)
# - pos: Position on the chromosome (numeric)
# - ref: Reference allele (character)
# - alt: Alternative allele (character)
# - ref.mean: Mean signal for the reference allele (numeric)
# - alt.mean: Mean signal for the alternative allele (numeric)
# - pvalue: p-value for the association (numeric)
#
# The script defines two functions: `check_format` and `plot_tsv_file`. The `check_format` function
# ensures that the input data meets the expected format, while the `plot_tsv_file` function generates
# a plot for the provided genomic region using the ggplot2 library.
#
#-----------------------------------------------------------------------------------------------

check_format_reply <<- NULL
n_user_defined_tsv <<- 0

# make y labels smaller if too big
wrap_string <- function(str, x) {
    # Split the string into chunks of x characters
    split_str <- substring(str, seq(1, nchar(str), by = x), seq(x, nchar(str) + x - 1, by = x))
    # Combine the chunks with a newline separator
    paste(split_str, collapse = "-\n")
}

# Define a function to plot the MPRA data from the TSV file
plot_tsv_file <- function(tsv_file, tsv_name, chr, pos1, pos2, pos) {
    
    tsv_name = wrap_string(tsv_name,30)

    # Load the data from the TSV file
    df <- fread(tsv_file)
        message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", paste0("The ", tsv_name, " format is correct!")))
        # If there are rows to plot, proceed with plot generation
        if (nrow(df) != 0) {
            # Rename columns for easier use in ggplot
            colnames(df) <- c("chr", "pos", "ref", "alt", "ref.mean", "alt.mean", "pval")

            # Make sure the numerical cols are also numerical            
            df$pos = as.numeric(df$pos)
            df$ref.mean = as.numeric(df$ref.mean)
            df$alt.mean = as.numeric(df$alt.mean)
            df$pval = as.numeric(df$pval)

            # Filter the data for the target chromosome and position range
            df <- df[df$chr == chr & df$pos >= pos1 & df$pos <= pos2, ]

            # Calculate allele information based on reference and alternative mean signals
            df$max.allele <- ifelse(df$ref.mean > df$alt.mean,
                "Loss of function\n( REF > ALT )",
                "Gain of function\n( ALT > REF )"
            )

            # Create ggplot object for SuRE signal visualization
            df_plot <- ggplot(df, aes(
                x = pos,
                y = ref.mean,
                xend = pos,
                yend = alt.mean,
                alpha = -log10(pval),
                color = max.allele,
                size = -log10(pval)
            )) +
                geom_segment() + # Use segments to show reference and alternative signals
                guides(size = "none", alpha = "none") + # Remove legends for size and alpha
                scale_color_manual(
                    values = c("#3c5379", "#64BEA5"),
                    breaks = c(
                        "Gain of function\n( ALT > REF )",
                        "Loss of function\n( REF > ALT )"
                    ),
                    guide = guide_legend(title = "")
                ) +
                scale_alpha_continuous(
                    limits = c(0, 1)
                ) +
                theme_bw() + # Use a clean black-and-white theme
                theme(
                    legend.position = "right",
                    axis.title.x = element_blank() # Remove x-axis title
                ) +
                ylab(tsv_name) + # Use file name as plot label
                xlim(c(pos1, pos2)) + # Limit the x-axis to the specified position range
                coord_cartesian(xlim = c(pos1, pos2)) # Ensure proper display within the range
        } else {
            # If no data for the given region, generate an empty plot
            df_plot <- ggplot() +
                geom_line(aes(x = c(pos1, pos1), y = c(0, 1)), linetype = "dashed") +
                labs(fill = "") +
                xlab("") +
                ylab(tsv_name) +
                guides(color = FALSE) +
                coord_cartesian(xlim = c(pos1, pos1)) +
                theme_bw() +
                theme(
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank()
                )
        }
    # Return the plot object
    if(is.na(pos)){
        return(df_plot)
    } else {
        uploaded_query_snps = df[df$chr == query_snps$CHROM & df$pos == query_snps$POS & df$ref == query_snps$REF & df$alt == query_snps$ALT, ]
        if(nrow(uploaded_query_snps)>0){
            df_plot = df_plot +
                geom_segment(aes(x = uploaded_query_snps$pos, xend = uploaded_query_snps$pos, y = uploaded_query_snps$ref.mean, yend = uploaded_query_snps$alt.mean),
                    color = "yellow", size = 0.5, alpha = 1
                ) +
                geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
        } else {
            df_plot = df_plot +
                geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
        }
    }
}
