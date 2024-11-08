#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script handles the processing and rendering of user-defined genomic data files for visualization.
# It initializes and processes user-defined BigWig (.bw) and BED (.bed) files, generates various plots
# including SuRE plots, gene plots, BigWig plots, and user-defined BigWig and BED plots. Additionally,
# it manages the creation of JASPAR plots, ClinVar and gnomAD pages, and SuRE table data. The script
# sources different modules to achieve these tasks and provides functionality to render and download
# output plots and data tables.
#
# Script Workflow:
# 1. **User Defined BW+BEDs**: Initializes variables for user-defined files and processes them based
#    on their extensions (BigWig and BED). Handles empty lists appropriately.
# 2. **Rendering After Knowing chr pos pos1 pos2**:
#    - **SuRE Plot**: Sources and executes the script to generate SuRE plots.
#    - **Gene Plots**: Sources the script for generating gene plots.
#    - **Gene Expression**: Sources the script for gene expression plots and generates the plots.
#    - **BigWig Plot**: Sources the script for BigWig plot generation and combines multiple plots.
# 3. **User Defined BigWig Plot**: Handles the creation of plots for user-defined BigWig and BED files.
#    - Includes conditional checks and dynamic plot generation based on the presence of files.
# 4. **JASPAR Plot**: Sources the script for generating JASPAR plots.
# 5. **ClinVar Page**: Sources the script for displaying ClinVar page information.
# 6. **gnomAD Page**: Sources the script for displaying gnomAD page information.
# 7. **SuRE Table**: Sources the script for generating the SuRE table.
# 8. **Render Outputs**: Sources the script responsible for rendering all outputs.
# 9. **Download Outputs**: Sources the script for handling the download of output data.
#
#-----------------------------------------------------------------------------------------------

########################
########################
# User Defined BW+BEDs #
########################
########################

# Initialize variables for user-defined files
user_defined <- NULL
user_defined_bw <- NULL
user_defined_bed <- NULL

# Populate the user_defined list with selected files
user_defined <- as.list(selected_files())
if (length(user_defined) > 0) {
    names(user_defined) <- lapply(user_defined, basename)

    # Split the list based on file extensions
    user_defined_bw <- user_defined[grepl("\\.bw$", user_defined)]
    user_defined_bw <- user_defined[grepl("\\.bigwig$", user_defined)]
    user_defined_bw <- user_defined[grepl("\\.bigWig$", user_defined)]
    user_defined_bed <- user_defined[grepl("\\.bed$", user_defined)]

    # Check if the lists are empty and set them to NULL if they are
    if (length(user_defined_bw) == 0) {
        user_defined_bw <- NULL
    }
    if (length(user_defined_bed) == 0) {
        user_defined_bed <- NULL
    }
}

# Rendering after knowing chr pos pos1 pos2
#############
# SuRE plot #
#############

# Source the script for SuRE plot generation
source(file.path(appDIR, "modules", "plot", "sureplot.R"), local = TRUE )
log.this.info("SuRE Plot : Plotted ✅")

##############
# Gene plots #
##############

# Source the script for gene plots and plotly visualizations
source(file.path(appDIR, "modules", "plot", "geneplot.R"), local = TRUE )
log.this.info("Gene Plot : Plotted ✅")

###################
# Gene expression #
###################

# Source the script for gene expression plots and generate them
source(file.path(appDIR, "modules", "plot", "gene.expression.R"), local = TRUE )
genes <- unique(genedata$GENE)
gene.expression.plots <- get.expression.plot(genes)
log.this.info("Gene Expression Plot : Plotted ✅")

###############
# Bigwig plot #
###############

# Source the script for BigWig plot generation
source(file.path(appDIR, "modules", "plot", "bigwigplot.R"))

# Combine various plots related to SNPs, genes, and BigWig data
SuREbigwig <<- TRUE
SuREbigwigplots <- plot.bigwigs(chr, pos1, pos2, pos, bigwigs, "#00aeffb3", "common")

SuREbigwig <<- FALSE
bigwig.plot <- wrap_plots(
    static.snp.plot,
    gene.plot,
    SuREbigwigplots,
    plot.bigwigs(chr, pos1, pos2, pos, AC16ATACbw, "#645200", "common"),
    ncol = 1, heights = c(0.2, 0.1, 0.7, 0.1)
)
log.this.info("SuRE Bigwig Tracks : Plotted SuRE bigwig and AC16 ATACseq plot ✅")


############################
# User defined Bigwig plot #
############################

# Source the script for BED plot generation
source(file.path(appDIR, "modules", "plot", "bedplot.R"))
SuREbigwig <<- FALSE
# Initialize the plot for user-defined BigWig and BED files
user.defined.plot <- NULL
if (!is.null(user_defined_bw)) {
    
    log.this.info(paste0("User defined bigwigs to be plotted :\n", paste(user_defined_bw, collapse = "\n")))
    n_user_defined_bw <- length(user_defined_bw)
    user.defined.plot <- wrap_plots(
        plot.bigwigs(chr, pos1, pos2, pos, user_defined_bw, "#ff0084b3", "common"),
        static.snp.plot,
        gene.plot,
        ncol = 1, heights = c(0.2 * n_user_defined_bw, 0.2, 0.1)
    )
    log.this.info("User Defined Plots: Plotted user defined bigwigs ✅")
    
    if (!is.null(user_defined_bed)) {

        log.this.info(paste0("User defined bed to be plotted :\n", paste(user_defined_bed, collapse = "\n")))
        n_user_defined_bed <- length(user_defined_bed)
        user.defined.plot <- wrap_plots(
            plot.bigwigs(chr, pos1, pos2, pos, user_defined_bw, "#ff0084b3", "common"),
            plot.beds(chr, pos1, pos2, pos, user_defined_bed),
            static.snp.plot,
            gene.plot,
            ncol = 1, heights = c(0.2 * n_user_defined_bw, 0.1 * n_user_defined_bed, 0.2, 0.1)
        )
        log.this.info("User Defined Plots : Plotted user defined beds ✅")
        
    }
} else {
    if (!is.null(user_defined_bed)) {

        log.this.info(paste0("User defined bed to be plotted :\n", paste(user_defined_bed, collapse = "\n")))
        n_user_defined_bed <- length(user_defined_bed)
        user.defined.plot <- wrap_plots(
            plot.beds(chr, pos1, pos2, pos, user_defined_bed),
            static.snp.plot,
            gene.plot,
            ncol = 1, heights = c(0.1 * n_user_defined_bed, 0.3, 0.2)
        )
        log.this.info("User Defined Plots : Plotted user defined beds ✅") 
        
    } else {
        # Create an empty plot indicating no user-defined files
        user.defined.plot <- ggplot() +
            theme_minimal() +
            annotate("text", x = 0.5, y = 0.5, label = "No bigwigs or bed files uploaded", size = 6, hjust = 0.5, vjust = 0.5) +
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.grid = element_blank()
            )
        log.this.info("User Defined Plots : User has not defined beds or bigwigs. Plotted an empty plot. ⚠️")
    }
}

###################################
# JASPAR plot #
###################################

# Source the script for JASPAR plot generation
source(file.path(appDIR, "modules", "plot", "jasparplot.R"), local = TRUE )

###################################
# ClinVar page #
###################################

# Source the script for ClinVar page information
source(file.path(appDIR, "modules", "url", "clinvar.R"), local = TRUE )

###################################
# gnomAD page #
###################################

# Source the script for gnomAD page information
source(file.path(appDIR, "modules", "url", "gnomad.R"), local = TRUE )

##############
# SuRE table #
##############

# Source the script for SuRE table data
source(file.path(appDIR, "modules", "table", "get.snp.table.R"), local = TRUE )

##################
# Render Outputs #
##################

# Source the script for rendering all outputs
source(file.path(appDIR, "main", "render.outputs.R"), local = TRUE )

####################
# Download Outputs #
####################
# Source the script for handling the download of output data
source(file.path(appDIR, "main", "downloadData.R"), local = TRUE )
