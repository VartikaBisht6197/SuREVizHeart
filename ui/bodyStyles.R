#-----------------------------------------------------------------------------------------------
# Function to Apply Custom CSS Styles to Shiny Application
# Author: Vartika Bisht
# Date: 19.08.2024
#-----------------------------------------------------------------------------------------------
# The 'bodyStyles' function defines custom CSS styles to be applied across the Shiny application.
# These styles are intended to enhance the visual appearance of various UI components, such as
# headers, sidebars, and buttons. The function uses inline CSS within the `tags$style` function
# to ensure that the styles are applied correctly.
#-----------------------------------------------------------------------------------------------

# Custom CSS for enhancing UI appearance
bodyStyles <- function() {
  tags$head(
    tags$style(HTML("
      body {
        font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
        background-color: #f9f9f9;
      }
      .main-header .logo {
        font-size: 1.5em;
      }
      .skin-blue .main-header .logo {
        background-color: #f5f5f5;
        color: black !important;
      }
      .skin-blue .main-header .logo:hover {
        background-color: #e9e9e9;
      }
      .skin-blue .main-header .navbar {
        background-color: #f5f5f5;
      }
      .skin-blue .main-header .navbar .sidebar-toggle {
        color: black !important;
      }
      .skin-blue .main-sidebar {
        background-color: #f5f5f5;
      }
      .skin-blue .main-sidebar .sidebar .sidebar-menu .active a {
        background-color: #e9e9e9;
        border-left-color: #64BEA5;
      }
      .skin-blue .main-sidebar .sidebar .sidebar-menu a {
        color: black;
      }
      .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover {
        background-color: #e9e9e9;
        color: black;
      }
      .content-wrapper, .right-side {
        background-color: #ffffff;
      }
      .box {
        border-top: 2px solid #64BEA5;
      }
      .btn-primary {
        background-color: #64BEA5;
        border-color: #64BEA5;
      }
      .btn-primary:hover {
        background-color: #56a08e;
        border-color: #56a08e;
      }
      .box-title {
        color: black;
      }
      .sidebar-menu .treeview-menu > li > a {
        white-space: normal !important;
        padding-left: 25px !important;
        color: black !important; /* Set non-highlighted text color */
        background-color: #dcdcdc !important;
      }
      .skin-blue .main-sidebar .sidebar .sidebar-menu .active > a {
        background-color: #e9e9e9 !important;
        color: black !important;
      }

      /* Custom Styles for File Input */
      .file-input-label-row {
        display: flex;
        align-items: center;
        gap: 10px;
      }
      .file-input-label-row .btn {
        flex: 1; /* Allow the Browse button to take up remaining space */
      }
      .selected-files {
        margin-top: 15px;
      }
      .file-upload-container {
        border: 2px dashed #64BEA5;
        padding: 20px;
        border-radius: 10px;
        background-color: #f9f9f9;
        transition: background-color 0.3s ease, border-color 0.3s ease;
        cursor: pointer;
        text-align: center;
        margin-top: 30px;
      }
      .file-upload-container:hover {
        background-color: #dff7ef;
        border-color: #56a08e;
      }
      .file-upload-container input[type='file'] {
        font-size: 1.1em;
        color: #555;
        font-weight: bold;
        background-color: transparent;
        border: none;
        box-shadow: none;
        text-align: center;
      }
      .file-upload-container .placeholder-text {
        font-size: 1.1em;
        font-weight: bold;
        color: #777;
      }
      .row-spacing {
          margin-top: 20px;
          }
    "))
  )
}
