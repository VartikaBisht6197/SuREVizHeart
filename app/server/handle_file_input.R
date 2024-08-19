#-----------------------------------------------------------------------------------------------
# File Input Handling and Management
# Author: Vartika Bisht
# Date: 19.08.2024
#-----------------------------------------------------------------------------------------------
# This script manages file inputs in a Shiny application. It allows users to select files using
# a file chooser, displays the selected files, and provides options to remove files from the
# selection. The selected file paths are stored and dynamically updated in the user interface.
#-----------------------------------------------------------------------------------------------

# Define the root directory for file selection
roots <- c(root = "/")

# Initialize the file chooser for input; allows selection of 'bed' and 'bw' file types
shinyFileChoose(input, "file", root = roots, filetypes = c("bed", "bw"))

# Define a reactive value to store the paths of selected files
selected_files <- reactiveVal(character())

# Observe changes in file input and update the reactive value 'selected_files' with the new paths
observeEvent(input$file, {
    # Parse the selected file paths
    file_paths <- parseFilePaths(roots, input$file)$datapath
    # Update 'selected_files' by adding the new file paths and ensuring uniqueness
    selected_files(unique(c(selected_files(), file_paths)))
})

# Render the UI elements to display selected files
output$selected_files <- renderUI({
    # Retrieve the current list of selected file paths
    paths <- selected_files()
    # If no files are selected, return NULL
    if (length(paths) == 0) {
        return(NULL)
    } else {
        # Create UI elements for each selected file
        paths_html <- lapply(paths, function(path) {
            file_name <- basename(path)
            fluidRow(
                column(10, textOutput(paste0("file_", file_name), container = span)),
                column(2, actionButton(inputId = paste0("remove_file_", file_name), label = "Remove"))
            )
        })
        # Combine all UI elements into a single tag list
        do.call(tagList, paths_html)
    }
})

# Render the file names dynamically
observe({
    # Retrieve the current list of selected file paths
    paths <- selected_files()
    # For each file, create a text output for its name
    lapply(paths, function(path) {
        file_name <- basename(path)
        output[[paste0("file_", file_name)]] <- renderText({
            as.character(file_name)
        })
    })
})

# Handle the removal of files dynamically
# Note: ignoreInit = TRUE prevents the event handler from being triggered when initially created
observe({
    # Retrieve the current list of selected file paths
    paths <- selected_files()
    # For each file, set up an event handler for the "Remove" button
    lapply(paths, function(path) {
        file_name <- basename(path)
        observeEvent(input[[paste0("remove_file_", file_name)]],
            {
                # Update the list of selected files by removing the chosen file
                new_paths <- selected_files()
                new_paths <- new_paths[new_paths != path]
                selected_files(new_paths)
            },
            ignoreInit = TRUE
        )
    })
})
