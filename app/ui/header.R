#-----------------------------------------------------------------------------------------------
# Define Header Function
# Author: Vartika Bisht
# Date: 19.08.2024
#-----------------------------------------------------------------------------------------------
# This function creates a header for the Shiny dashboard using the `dashboardHeader` function
# from the `shinydashboard` package. It includes a logo image and a title for the application.
# The header is styled to ensure that the logo and text are displayed prominently and align
# properly within the header space.
#-----------------------------------------------------------------------------------------------

header <- function() {
    dashboardHeader(
        # The title of the header is a combination of an image and text
        title = div(
            shiny::tags$img(
                src = "Annogen.logo.png", # Path to the logo image file
                height = 40, # Height of the logo image
                style = "margin-right: 10px;" # Inline style for right margin to space out the logo from text
            ),
            span(
                "SuRE Visualization App", # Application title text
                style = "color: black;" # Inline style to set text color to black
            )
        ),
        titleWidth = 400 # Sets the width of the header to 400 pixels
    )
}
