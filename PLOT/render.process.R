
########################
########################
# User Defined BW+BEDs #
########################
########################

# Initialize variables
user_defined = NULL
user_defined_bw = NULL
user_defined_bed = NULL

# Populate the user_defined list with selected files
user_defined = as.list(selected_files())
if (length(user_defined) > 0) {
names(user_defined) = lapply(user_defined, basename)

# Split the list based on file extensions
user_defined_bw <- user_defined[grepl("\\.bw$", user_defined)]
user_defined_bed <- user_defined[grepl("\\.bed$", user_defined)]

# Check if the lists are empty and set them to NULL if they are
if (length(user_defined_bw) == 0) {
    user_defined_bw <- NULL
}
if (length(user_defined_bed) == 0) {
    user_defined_bed <- NULL
}
}


# Rendering after you know chr pos pos1 pos2
#############
# SuRE plot #
#############

# snp.plot
# static.snp.plot
# suredata
source(file.path(plotDIR, "get.sureplot.R"), local = TRUE)

##############
# Gene plots #
##############

# gene.plot
# gene.plotly
source(file.path(plotDIR, "get.geneplot.R"), local = TRUE)

###################
# Gene expression #
###################
source(file.path(plotDIR, "get.gene.expression.R"), local = TRUE)
genes = unique(genedata$GENE)
gene.expression.plots = get.expression.plot(genes)

###############
# Bigwig plot #
###############

# Bigwigs Quite big
# source(file.path(plotDIR, "bigwigplot.R"))
# bigwig.plot <- wrap_plots(static.snp.plot,
#     gene.plot,
#     plot.bigwigs(chr, pos1, pos2, pos, bigwigs, "#00aeffb3", "common"),
#     plot.bigwigs(chr, pos1, pos2, pos, AC16ATACbw, "#645200", "common"),
#     ncol = 1, heights = c(0.2, 0.1, 0.7, 0.1) 
# )


############################
# User defined Bigwig plot #
############################

source(file.path(plotDIR, "bedplot.R"))
user.defined.plot = NULL
user.defined.plot <- NULL
if (!is.null(user_defined_bw)) {
    n_user_defined_bw = length(user_defined_bw)
    user.defined.plot <- wrap_plots(
        plot.bigwigs(chr, pos1, pos2, pos, user_defined_bw, "#ff0084b3", "common"),
        static.snp.plot,
        gene.plot,
        ncol = 1, heights = c(0.2*n_user_defined_bw, 0.2, 0.1)
    ) 

    if (!is.null(user_defined_bed)) {
        n_user_defined_bw = length(user_defined_bw)
        n_user_defined_bed = length(user_defined_bed)
        user.defined.plot <- wrap_plots(
            plot.bigwigs(chr, pos1, pos2, pos, user_defined_bw, "#ff0084b3", "common"),
            plot.beds(chr, pos1, pos2, pos, user_defined_bed),
            static.snp.plot,
            gene.plot,
            ncol = 1, heights = c(0.2*n_user_defined_bw, 0.1*n_user_defined_bed, 0.2, 0.1) 
        )
    }
} else {
    if (!is.null(user_defined_bed)) {
        n_user_defined_bed = length(user_defined_bed)
        user.defined.plot <- wrap_plots(
            plot.beds(chr, pos1, pos2, pos, user_defined_bed),
            static.snp.plot,
            gene.plot,
            ncol = 1, heights = c(0.1*n_user_defined_bed, 0.3, 0.2) 
        )
    } else {
        # Create an empty Plotly plot
        user.defined.plot <- ggplot() +
            theme_minimal() +
            annotate("text", x = 0.5, y = 0.5, label = "No bigwigs or bed files uploaded", size = 6, hjust = 0.5, vjust = 0.5) +
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.grid = element_blank()
            )
    }
}


###################################
# Selected Variant Overview plots #
###################################

# result_list
# jaspar.data
source(file.path(plotDIR, "selected.variant.overview.R"), local = TRUE)


##############
# SuRE table #
##############

# SNP.info
source(file.path(plotDIR, "get.snp.table.R"), local = TRUE)


##################
# Render Outputs #
##################

# render all outputs
source(file.path(plotDIR, "render.outputs.R"), local = TRUE)

####################
# Download Outputs #
####################

# Make one plot with all tracks
# Make one table with all data

# SNP.info
# list of bigwigs import
# Save the imported data as a new BigWig file
# export.bw(object = bigwigs_data, con = "path_to_new_bigwigs_file.bw")
granges_list <- lapply(bigwigs, function(x) import.bw(x, which = GRanges(seqnames = chr, ranges = IRanges(start = pos1, end = pos2))))
output$downloadData <- downloadHandler(
    filename = function() {
        paste("SuRE_Viz_data_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".zip", sep = "")
    },
    content = function(file) {
        # Create a temporary directory
        temp_dir <- tempdir()

        # Save SNP.info as a CSV
        if (nrow(SNP.info) > 0) {
            csv_path <- file.path(temp_dir, "SNP_info.csv")
            write.csv(SNP.info, csv_path, row.names = FALSE)
        }
        # Bigwigs Quite big
        # Save each GRanges object as a BEDGraph file
        # lapply(names(granges_list), function(gr_name) {
        #     gr_path <- file.path(temp_dir, paste0(gr_name, ".bedGraph"))
        #     export(granges_list[[gr_name]], con = gr_path, format = "bedGraph")
        # })

        # Save Jaspar info if available
        if (!is.null(jaspar.data)) {
            csv_path <- file.path(temp_dir, "JASPAR2020_info.csv")
            write.csv(jaspar.data, csv_path, row.names = FALSE)
        }

        # List of files to include in the zip
        files_to_zip <- c(csv_path, list.files(temp_dir, pattern = "\\.(csv|bedGraph)$", full.names = TRUE))

        # Zip the files
        zip::zipr(file, files = files_to_zip)
    }
)