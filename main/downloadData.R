#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script is responsible for generating a downloadable ZIP file containing various genomic
# data files relevant to variant visualization. The script handles the creation of GRanges objects
# from BigWig files for specific genomic regions, the export of SNP information, and optional JASPAR
# motif data. These data are then packaged into a ZIP file for user download, ensuring that all
# necessary files are included and organized efficiently.
#
# Key Functionalities:
# 1. **GRanges Objects Creation**: Converts BigWig data for a specified chromosome and range into GRanges objects.
# 2. **Download Handler**: Facilitates the download of a ZIP file that includes:
#    - SNP Information: Exported as a CSV file.
#    - BEDGraph Files: Generated from GRanges objects.
#    - JASPAR Motif Data: Included as an optional CSV file, if available.
# 3. **ZIP File Creation**: All relevant files are collected into a temporary directory and compressed
#    into a ZIP file, which is then made available for the user to download.
#-----------------------------------------------------------------------------------------------

# Define the download handler for creating and downloading the ZIP file
if(input_pass_test()){
    # Create a list of GRanges objects from BigWig files for a specific chromosome and range
    granges_list <- lapply(bigwigs, function(x) {
        rtracklayer::import.bw(x, which = GRanges(seqnames = chr, ranges = IRanges(start = pos1, end = pos2)))
    })

    output$downloadData <- downloadHandler(
        # Define the filename of the ZIP file with a timestamp
        filename = function() {
            paste("SuRE_Viz_data_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".zip", sep = "")
        },
        # Define the content of the ZIP file
        content = function(file) {
            # Create a temporary directory for storing intermediate files
            temp_dir <- tempdir()

            # Define the log file path in the temporary directory
            log_path <- file.path(temp_dir, "render.log")

            # Copy the initialized log file to the temporary directory
            writeLines(log.file, log_path)

            # Save SNP information as a CSV file if available
            if (nrow(SNP.info) > 0) {
                csv_path <- file.path(temp_dir, "SNP_info.csv")
                write.csv(SNP.info, csv_path, row.names = FALSE)
            }

            # Save each GRanges object as a BEDGraph file
            lapply(names(granges_list), function(gr_name) {
                gr_path <- file.path(temp_dir, paste0(gr_name, ".bedGraph"))
                rtracklayer::export(granges_list[[gr_name]], con = gr_path, format = "bedGraph")
            })

            # Save JASPAR motif information as a CSV file if available
            if (!is.null(jaspar.data)) {
                csv_path <- file.path(temp_dir, "JASPAR2020_info.csv")
                write.csv(jaspar.data, csv_path, row.names = FALSE)
            }

            # List all files in the temporary directory to be included in the ZIP
            files_to_zip <- c(list.files(temp_dir, pattern = "\\.(csv|bedGraph|log)$", full.names = TRUE))

            # Zip the files into the output ZIP file
            zip::zipr(file, files = files_to_zip)
        }
    )
    }else{
    output$downloadData <- downloadHandler(
        # Define the filename of the error file with a timestamp
        filename = function() {
            paste("error_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".log", sep = "")
        },
        content = function(file) {
            # write log into file
            writeLines(log.file, file)
            })

}