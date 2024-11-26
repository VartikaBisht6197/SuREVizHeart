#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
# Description:
# This Shiny module handles dynamic file upload and removal functionality. Users can upload files,
# and those files are displayed with an option to remove them. The removal and addition of files
# are managed dynamically with Reactivity. To ensure a seamless experience, dynamic observers are
# created and destroyed as needed for each file.
#-----------------------------------------------------------------------------------------------

# Reactive value to store user-uploaded files
user_defined_files <- reactiveValues(files = list())

# Reactive value to store dynamically created observers for file removal
remove_observers <- reactiveValues(observers = list())

# Observe file upload events
observeEvent(input$file, {
    # Check if files are uploaded and handle them
    if (!is.null(input$file) && nrow(input$file) > 0) {
        uploaded_files <- input$file

        # Iterate over uploaded files to process each one
        for (i in 1:nrow(uploaded_files)) {
            file_name <- uploaded_files$name[i] # Get the file name

            # If the file already exists, remove it to avoid duplicates
            if (file_name %in% names(user_defined_files$files)) {
                user_defined_files$files[[file_name]] <- NULL # Remove the file
                # Destroy any existing observer associated with this file
                if (!is.null(remove_observers$observers[[file_name]])) {
                    remove_observers$observers[[file_name]]$destroy()
                    remove_observers$observers[[file_name]] <- NULL
                }
            }

            # Store the uploaded file's metadata in the reactive values
            user_defined_files$files[[file_name]] <- list(
                path = uploaded_files$datapath[i],
                name = uploaded_files$name[i],
                size = uploaded_files$size[i],
                type = uploaded_files$type[i]
            )
        }

        # Dynamically update the UI to reflect the list of uploaded files
        output$selected_files <- renderUI({
            if (length(user_defined_files$files) == 0) {
                return(NULL) # No files to display
            } else {
                # Create a UI row for each uploaded file
                file_ui <- lapply(names(user_defined_files$files), function(file_name) {
                    file_info <- user_defined_files$files[[file_name]]
                    fluidRow(
                        column(8, span(file_info$name)), # Display file name
                        column(2, actionButton(inputId = paste0("remove_file_", file_name), label = "Remove")) # Add "Remove" button
                    )
                })
                do.call(tagList, file_ui) # Combine all UI rows into a single element
            }
        })

        # Create dynamic observers for removing files
        lapply(names(user_defined_files$files), function(file_name) {
            # Only create an observer if it doesn't already exist
            if (is.null(remove_observers$observers[[file_name]])) {
                remove_observers$observers[[file_name]] <- observeEvent(input[[paste0("remove_file_", file_name)]],
                    {
                        # Handle file removal when "Remove" button is clicked
                        user_defined_files$files[[file_name]] <- NULL # Remove file metadata
                        # Clean up reactive values by removing NULL entries
                        user_defined_files$files <- user_defined_files$files[!sapply(user_defined_files$files, is.null)]

                        # Destroy the observer for this file to free up resources
                        remove_observers$observers[[file_name]]$destroy()
                        remove_observers$observers[[file_name]] <- NULL

                        # Re-render the UI to reflect the updated file list
                        output$selected_files <- renderUI({
                            if (length(user_defined_files$files) == 0) {
                                return(NULL) # No files to display
                            } else {
                                # Recreate the UI with the updated list of files
                                file_ui <- lapply(names(user_defined_files$files), function(file_name) {
                                    file_info <- user_defined_files$files[[file_name]]
                                    fluidRow(
                                        column(8, span(file_info$name)), # Display file name
                                        column(2, actionButton(inputId = paste0("remove_file_", file_name), label = "Remove")) # Add "Remove" button
                                    )
                                })
                                do.call(tagList, file_ui) # Combine all UI rows into a single element
                            }
                        })
                    },
                    ignoreInit = TRUE
                ) # Prevent the observer from triggering on creation
            }
        })
    }
})
