# SuRE Visualization Shiny App
# chr12:128797635
# esco2
# Install required packages if not already installed
# install.packages(c("formattable", "plotly", "kableExtra", "data.table", "shiny"))

# Load required libraries
library(formattable)
library(plotly)
library(patchwork)
library(ggplotify)
library(universalmotif)
library(kableExtra)
library(data.table)
library(shinyalert)
library(shinycssloaders)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(here)

#options(shiny.fullstacktrace = FALSE)
RscriptDIR <- here()

# Function to convert values to kilo bases
to_kb <- function(x) {
  paste0(x / 1000, "kb")
}

to_number <- function(x) {
  num_str <- gsub("kb", "", x, ignore.case = TRUE)
  num <- as.numeric(num_str)
  num * 1000
}

# Define custom ticks and labels
custom_ticks <- c(1000, 5000, 10000, 20000, 50000, 100000)
custom_tick_labels <- to_kb(custom_ticks)
motifvh <- 1000
motifs <- 1

# Define the UI
ui <- dashboardPage(
  dashboardHeader(
    title = div(
      shiny::tags$img(src = "Annogen.logo.png", height = 40, style = "margin-right: 10px;"),
      span("SuRE Visualization App", style = "color: black;")
    ),
    titleWidth = 400 
  ),
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      # Icon list https://fontawesome.com/icons
      menuItem("Functional Impact Assessment", tabName = "impact_assessment", icon = icon("chart-bar")),
      menuItem("Variant Data Overview", tabName = "variant_summary", icon = icon("table")),
      menuItem("Gene Expression Overview", tabName = "gene_expression_summary", icon = icon("chart-simple")),
      menuItem("Integrating Complementary Datasets", tabName = "sure_profiles", icon = icon("dna")),
      menuItem("Highlighted Variant", tabName = "highlighted_variant", icon = icon("vial"), style = "background-color: #e9e9e9;",
               menuSubItem("Transcription Factor Disruption Analysis", tabName = "jaspar"),
               menuSubItem("Genome Aggregation Database (gnomAD) Viewer", tabName = "gnomad"),
               menuSubItem("ClinVar Viewer", tabName = "clinvar") 
      ),
      menuItem("Uploaded Tracks Viewer", tabName = "comparative_study", icon = icon("project-diagram"))
    )
  ),
  dashboardBody(
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
      "))
    ),
    # Show the error message as a warning
    verbatimTextOutput("warning_output"),
    fluidRow(
      box(
        width = 12, status = "primary", solidHeader = TRUE, title = "Query the human genome",
        fluidRow(
          column(6, textInput("search_query", "Search", value = "chr12:128797635", placeholder = "Enter chr:pos or gene name")),
          column(6, shinyWidgets::sliderTextInput("flank", "Flank", selected = "10kb" , choices = lapply(seq(1000, 100000, 1000), to_kb)) )),
        fluidRow(
          column(12,
            style = "display: flex; align-items: center;",
            actionButton("runButton", "Render", icon = icon("play"), class = "btn-primary", style = "margin-left: 10px;"),
            shinyFilesButton("file", "Upload Files", title = "Select files", multiple = TRUE, class = "btn-primary", icon = icon("upload"), style = "margin-left: 10px;"),
            downloadButton("downloadData", "Download", class = "btn-primary", style = "margin-left: 10px;", icon = icon("download"))),
          column(12, uiOutput("selected_files"), style = "margin-top: 10px;")
        )
      )
    ),
    tabItems(
      tabItem(
        tabName = "impact_assessment",
        box(
          width = 12, status = "primary", solidHeader = TRUE, title = "Functional Impact Assessment",
          includeHTML(file.path(appDIR,"PLOT","SuREViz.html")),
          div(
            style = "overflow: auto;",
            withSpinner(plotlyOutput("plotTab1", height = "800px", width = "100%"))
          )
        )
      ),
      tabItem(
        tabName = "variant_summary",
        box(
          width = 12, status = "primary", solidHeader = TRUE, title = "Variant Data Overview",
          includeHTML(file.path(appDIR,"PLOT","variant.summary.html")),
          div(
            style = "overflow: auto;",
            withSpinner(htmlOutput("tableTab3"))
          )
        )
      ),
      tabItem(
        tabName = "gene_expression_summary",
        box(
          width = 12, status = "primary", solidHeader = TRUE, title = "Gene Expression Overview",
          includeHTML(file.path(appDIR,"PLOT","gene.expression.summary.html")),
          div(
            style = "overflow: auto;",
            withSpinner(plotOutput("plotTab3", width = "100%", height = "1000px"))
          )
        )
      ),
      tabItem(
        tabName = "sure_profiles",
        box(
          width = 12, status = "primary", solidHeader = TRUE, title = "Integrating Complementary Datasets with SuRE Data",
          includeHTML(file.path(appDIR,"PLOT","SuREprofileView.html")),
          div(
            style = "overflow: auto;",
            withSpinner(plotOutput("plotTab2", width = "100%", height = "1000px"))
          )
        )
      ),
      tabItem(
        tabName = "jaspar",
        box(
          width = 12, status = "primary", solidHeader = TRUE, title = "Transcription Factor Disruption Analysis",
          includeHTML(file.path(appDIR,"PLOT","TF.impact.html")),
          div(
            style = "overflow: auto; width: 100%; text-align: center; margin-bottom: 20px;",
            withSpinner(plotOutput("jasparplot", width = "100%", height = paste0(motifvh,"px")))
          ),
          div(
            style = "overflow: auto; width: 100%; margin-top: 10px; text-align: center; margin-top: 20px;",
            withSpinner(tableOutput("fimo_table"))
          )
        )
      )
,
      tabItem(
        tabName = "gnomad",
        box(
          width = 12, status = "primary", solidHeader = TRUE, title = "Genome Aggregation Database (gnomAD) Viewer",
          includeHTML(file.path(appDIR,"PLOT","population.genetics.html")),
          div(
            style = "overflow: auto;",
            withSpinner(htmlOutput("gnomad"))
          )
        )
      ),
      tabItem(
        tabName = "clinvar",
        box(
          width = 12, status = "primary", solidHeader = TRUE, title = "ClinVar Viewer ",
          includeHTML(file.path(appDIR,"PLOT","ClinVar.html")),
          div(
            style = "overflow: auto;",
            withSpinner(htmlOutput("clintable"))
          )
        )
      ),
      tabItem(
        tabName = "comparative_study",
        box(
          width = 12, status = "primary", solidHeader = TRUE, title = "Uploaded Tracks Viewer",
          div(
            style = "overflow: auto;",
            withSpinner(plotOutput("plotTab5", width = "100%", height = "1000px"))
          )
        )
      )
    )
  )
)




# Render memory
# Wipe all memory before render when you press render again
render_mem = 0

# Define server function
server <- function(input, output, session) {
  
  # Hardcoded values
  source(file.path(appDIR,"PLOT","HardcodedPaths.R"))
  
  # Define reactive values
  hover_reactive <- reactiveVal()
  chr_value <- reactiveVal(NULL)
  flank_value <- reactiveVal(NULL)
  pos_value <- reactiveVal(NULL)
  variant_view <- reactiveVal(NULL)
  alert_confirmed <- reactiveVal(FALSE)

  # Initialize the alert_confirmed reactive value
  observe({
    # Reset on app start or when necessary
    alert_confirmed(FALSE)
  })

  ######################
  ######################
  # Observe File Input #
  ######################
  ######################
  
  # Roots
  roots <- c(root = "/")

  shinyFileChoose(input, "file", root = roots, filetypes = c("bed","bw"))
  # shinyFileChoose(input, "file", root = roots)

  # Reactive value to store selected file paths
  selected_files <- reactiveVal(character())

  # Observe file input and update selected_files
  observeEvent(input$file, {
    file_paths <- parseFilePaths(roots, input$file)$datapath
    selected_files(unique(c(selected_files(), file_paths)))
  })

  # Render selected files UI
  output$selected_files <- renderUI({
    paths <- selected_files()
    if (length(paths) == 0) {
      return(NULL)
    } else {
      paths_html <- lapply(paths, function(path) {
        file_name <- basename(path)
        fluidRow(
          column(10, textOutput(paste0("file_", file_name), container = span)),
          column(2, actionButton(inputId = paste0("remove_file_", file_name), label = "Remove"))
        )
      })
      do.call(tagList, paths_html)
    }
  })

  # Render the file names
  observe({
    paths <- selected_files()
    lapply(paths, function(path) {
      file_name <- basename(path)
      output[[paste0("file_", file_name)]] <- renderText({
        as.character(file_name)
      })
    })
  })

  # Dynamically handle the removal of files
  # Add ignoreInit = TRUE in observeEvent to prevent the handlers from being triggered when they are first created.
    observe({
      paths <- selected_files()
      lapply(paths, function(path) {
        file_name <- basename(path)
        observeEvent(input[[paste0("remove_file_", file_name)]],
          {
            new_paths <- selected_files()
            new_paths <- new_paths[new_paths != path]
            selected_files(new_paths)
          },
          ignoreInit = TRUE
        )
      })
    })

#######################  #######################  #######################  #######################  #######################  #######################  #######################
#######################  #######################  #######################    Observe Render button  #######################  #######################  #######################
#######################  #######################  #######################  #######################  #######################  #######################  #######################

  
  observeEvent(input$runButton, {
      
      # Retrieve input values
      search_query <- input$search_query
      flank <- to_number(input$flank)

      print(search_query)

      # reformat given search_query
      source(file.path(plotDIR, "reformat.input.r"), local = TRUE)
  
      # Check if alert was confirmed
      if (alert_confirmed()) {
        # Only update values and continue if the alert was confirmed
        chr_value(chr)
        flank_value(flank)
        pos_value(pos)

        # Source the processing script
        source(file.path(plotDIR, "render.process.R"), local = TRUE)
      } else {
        print("Alert not confirmed. Skipping further processing.")
      }

     })
  
  
#######################  #######################  #######################  #######################  #######################  #######################  #######################
#######################  #######################  #######################      Observe Click        #######################  #######################  #######################
#######################  #######################  #######################  #######################  #######################  #######################  #######################
  
  observe({
    
    click_data_val <- event_data("plotly_click")
    
    # Change render memory to the actual position of the variant
    # If the click is null, or none selected then everything remains
    # the same and the render_mem does not let it 
    # https://stackoverflow.com/a/39907676
    new_render_mem <- ifelse(is.null(click_data_val),"0",click_data_val$x)
    
    if(render_mem!=new_render_mem) {
      
      render_mem <<- new_render_mem 
      
      
      if (!is.null(click_data_val) && render_mem != 0 && !is.integer(click_data_val$y) ) {
  
        # Access chr and flank values
        chr <- chr_value()
        flank <- flank_value()
        
        # Extract the data associated with the clicked point
        new_pos <- click_data_val$x
        
        # Update the pos value
        pos_value(new_pos)
        
        # Reset click_data and user response to NULL
        hover_reactive(NULL)
        click_data_val <- NULL
        
        # Display a confirmation dialog
        shinyalert::shinyalert(
          title = "Plotting the selected SNP",
          text = paste0("Make new plots for selected position \n", chr, ":", new_pos),
          type = "info",
          showCancelButton = TRUE,
          cancelButtonText = "No",
          confirmButtonText = "Yes")
        
      }}
  })

}

# Run the Shiny app
shinyApp(ui, server)
