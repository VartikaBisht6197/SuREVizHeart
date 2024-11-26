#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script processes and validates user-uploaded files. It checks if the files are in the
# correct format for TSV, BigWig, and BED file types. If the files are valid, it prints out the
# file details; otherwise, it flags errors and prompts the user with appropriate messages.
# Each file type is validated against specific rules:
# - TSV files: Must have exactly 7 columns with specific data formats.
# - BigWig files: Must be valid binary files.
# - BED files: Must have at least 3 columns and no headers.
#
#-----------------------------------------------------------------------------------------------

# Initialize user input for file types
user_defined_tsv <- NULL
user_defined_bw <- NULL
user_defined_bed <- NULL

# Function to print files of a specific type
print_files <- function(file_list, file_type) {
    if (length(file_list) > 0) {
        message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", paste("The following", file_type, "files have been uploaded:")))
        for (file_key in names(file_list)) {
            message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", paste("  -", file_key, "is a file saved at", file_list[[file_key]])))
        }
    } else {
        message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", paste("No", file_type, "files provided.")))
    }
}

# Iterate through the uploaded files and classify them into corresponding lists
for (file_name in names(user_defined_files$files)) {
    # Extract the file's metadata
    file_info <- user_defined_files$files[[file_name]]

    # Extract file extension and base name (without extension)
    file_ext <- tools::file_ext(file_info$name) # Get file extension
    file_base <- tools::file_path_sans_ext(file_info$name) # Get file name without extension

    # Add file paths to the appropriate list based on extension
    if (file_ext == "tsv") {
        user_defined_tsv[[file_base]] <- file_info$path # Add to TSV list
    } else if (file_ext == "bw") {
        user_defined_bw[[file_base]] <- file_info$path # Add to BigWig (BW) list
    } else if (file_ext == "bed") {
        user_defined_bed[[file_base]] <- file_info$path # Add to BED list
    }
}

# Print all the uploaded files categorized by type
print_files(user_defined_tsv, "TSV")
print_files(user_defined_bed, "BED")
print_files(user_defined_bw, "BigWig")

# Function to handle the validation error and update `input_pass_test`
handle_validation_error <- function(error_message, alert_message) {
    message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", error_message, "❌"))
    shinyalert(
        title = "Invalid input format",
        text = alert_message
    )
    # Use reactive assignment to set the flag for validation failure
    input_pass_test(FALSE)
}

# Validate user inputs
validate_inputs <- function(tsv_list, bw_list, bed_list) {
    # Initialize a flag to track if inputs pass validation
    validation_status <- TRUE

    # Validate TSV files using fread from data.table package
    for (file_name in names(tsv_list)) {
        tryCatch(
            {
                tsv_path <- tsv_list[[file_name]]
                # Using fread to read the TSV file (fast and efficient)
                tsv_data <- fread(tsv_path)

                # Check for exactly 7 columns
                if (ncol(tsv_data) != 7) {
                    stop(paste("TSV file does not have exactly 7 columns. Found", ncol(tsv_data), "columns in file", file_name))
                }

                # Check the format of each column
                if (!all(grepl("^chr[0-9XY]+$", tsv_data[[1]]))) {
                    stop(paste("TSV file has invalid chromosome format in column 1 in file", file_name))
                }
            },
            error = function(e) {
                validation_status <<- FALSE
                handle_validation_error(
                    paste("Error in file", file_name, ":", e$message),
                    e$message
                )
                # Exit the loop after first error
                return()
            }
        )
    }


    # Validate BigWig files
    for (file_name in names(bw_list)) {
        tryCatch(
            {
                bw_path <- bw_list[[file_name]]

                # Attempt to open the BigWig file without specifying a chromosome or range
                bw_data <- import.bw(bw_path)

                # If no error occurs, the BigWig file is valid
                message(paste("BigWig file", file_name, "is valid. ✅"))
            },
            error = function(e) {
                validation_status <<- FALSE
                # Show the error message to the user
                handle_validation_error(
                    paste("Error in file", file_name, ":", e$message),
                    paste("Invalid BigWig file format in file", file_name, ": ", e$message)
                )
                # Exit the loop after first error
                return()
            }
        )
    }


    # Validate BED files
    for (file_name in names(bed_list)) {
        tryCatch(
            {
                bed_path <- bed_list[[file_name]]
                bed_data <- fread(bed_path)

                # Check for at least 3 columns
                if (ncol(bed_data) < 3) {
                    stop("BED file has less than 3 columns.")
                }
                # Check the format of each column
                if (!all(grepl("^chr[0-9XY]+$", bed_data[[1]]))) {
                    stop(paste("BED file has invalid chromosome format in column 1 in file", file_name))
                }
            },
            error = function(e) {
                validation_status <<- FALSE
                handle_validation_error(
                    paste("Error in file", file_name, ":", e$message),
                    e$message
                )
                # Exit the loop after first error
                return()
            }
        )
    }

    # If everything passes, validate the files
    if (validation_status) {
        input_pass_test(TRUE)
        message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "All inputs passed validation. ✅"))
    } else {
        # Explicitly set input_pass_test(FALSE) if there was an error
        input_pass_test(FALSE)
    }
}

# Call the validate function for all file types
validate_inputs(user_defined_tsv, user_defined_bw, user_defined_bed)
