library(shiny)
# The following packages are required for this app.
required_packages <- c("readxl", "pheatmap", "dplyr", "tibble")

# Check if the packages are installed. If not, install them.
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

library(readxl)
library(pheatmap)
library(dplyr)
library(tibble)

# Increase the maximum file upload size to 10.2 MB (approximately)
options(shiny.maxRequestSize = 10.2 * 1024^2)

ui <- fluidPage(
  titlePanel("Gene Expression Heatmap Generator"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("countFile", "Upload Count Data (.csv or .xlsx)", accept = c(".csv", ".xlsx")),
      selectInput("geneNameColumn", "Select Column Containing Gene Names", choices = NULL), # Added dropdown for gene name column
      textAreaInput("geneList", "Enter Gene List (one gene per line)", ""),
      selectInput("initialSample", "Initial Sample", choices = NULL),
      actionButton("addSample", "Add Sample"),
      uiOutput("additionalSamples"),
      checkboxInput("clusterRows", "Cluster Rows", TRUE),
      checkboxInput("clusterCols", "Cluster Columns", TRUE),
      
      # Added UI for custom color palette
      selectInput("colorPalette", "Color Palette",
                  choices = c("RdBu", "viridis", "magma", "plasma", "inferno", "Custom"),
                  selected = "RdBu"),
      conditionalPanel(
        condition = "input.colorPalette == 'Custom'",
        textInput("negColor", "Color for <0", value = "blue"),
        textInput("zeroColor", "Color for 0", value = "white"),
        textInput("posColor", "Color for >0", value = "red")
      ),
      
      sliderInput("plotWidth", "Plot Width (pixels)", min = 400, max = 1600, value = 800),
      sliderInput("plotHeight", "Plot Height (pixels)", min = 400, max = 1600, value = 800),
      actionButton("generateHeatmap", "Generate Heatmap"),
      downloadButton("downloadHeatmap", "Download Heatmap")
    ),
    
    mainPanel(
      plotOutput("heatmapPlot"),
      textOutput("errorMessage")
    )
  )
)

server <- function(input, output, session) {
  # Update sample choices and gene name column when file is uploaded
  observeEvent(input$countFile, {
    tryCatch({
      file_path <- input$countFile$datapath
      if (tools::file_ext(file_path) == "xlsx") {
        count_data <- read_excel(file_path)
      } else {
        count_data <- read.csv(file_path)
      }
      
      # Get column names for gene name selection and sample selection
      column_names <- names(count_data)
      updateSelectInput(session, "geneNameColumn", choices = column_names) # Update gene name column choices
      sample_choices <- column_names
      updateSelectInput(session, "initialSample", choices = sample_choices)
      
      # Store sample choices for later use
      session$userData$sample_choices <- sample_choices
      session$userData$column_names <- column_names # Store column names
      
    }, error = function(e) {
      output$errorMessage <- renderText(paste("Error reading file:", e$message))
    })
  })
  
  # Reactive value to store added samples
  added_samples <- reactiveVal(NULL)
  
  # Render UI for additional samples
  output$additionalSamples <- renderUI({
    if (!is.null(added_samples())) {
      lapply(1:length(added_samples()), function(i) {
        selectInput(paste0("addedSample", i), paste0("Add Sample ", i),
                    choices = session$userData$sample_choices, selected = added_samples()[i])
      })
    }
  })
  
  # Observe add sample button clicks
  observeEvent(input$addSample, {
    if (is.null(added_samples())) {
      added_samples(session$userData$sample_choices[1])
    } else {
      added_samples(c(added_samples(), session$userData$sample_choices[1]))
    }
  })
  
  heatmap_data <- reactiveVal(NULL) # Reactive value to store heatmap data
  
  observeEvent(input$generateHeatmap, {
    output$errorMessage <- renderText({
      req(input$countFile, input$geneList, input$initialSample, input$geneNameColumn) # Ensure geneNameColumn is required
      
      tryCatch({
        file_path <- input$countFile$datapath
        if (tools::file_ext(file_path) == "xlsx") {
          count_data <- read_excel(file_path)
        } else {
          count_data <- read.csv(file_path)
        }
        
        gene_list <- trimws(unlist(strsplit(input$geneList, "\n")))
        selected_samples <- c(input$initialSample, sapply(1:length(added_samples()), function(i) input[[paste0("addedSample", i)]]))
        gene_name_column <- input$geneNameColumn # Get selected gene name column
        
        if (!gene_name_column %in% colnames(count_data)) {
          return(paste0("Error: Selected gene name column '", gene_name_column, "' not found."))
        }
        
        # Check if all selected samples exist
        if (!all(selected_samples %in% colnames(count_data))) {
          missing_samples <- setdiff(selected_samples, colnames(count_data))
          return(paste("Error: Sample(s) not found:", paste(missing_samples, collapse = ", ")))
        }
        
        plot_data <- count_data %>%
          filter(!!sym(gene_name_column) %in% gene_list) %>% # Use selected column
          select(gene_name_column, all_of(selected_samples)) %>%
          column_to_rownames(gene_name_column)
        
        if (nrow(plot_data) == 0) {
          return("Error: No genes from the list found in the data.")
        }
        
        # Select numeric columns
        numeric_data <- plot_data %>%
          select_if(is.numeric)
        
        heatmap_data(numeric_data) #Store the data for download.
        
        # Color palette logic
        color_palette <- switch(input$colorPalette,
                                "RdBu" = colorRampPalette(c("blue", "white", "red"))(100),
                                "viridis" = viridis::viridis(100),
                                "magma" = viridis::magma(100),
                                "plasma" = viridis::plasma(100),
                                "inferno" = viridis::inferno(100),
                                "Custom" = {
                                  # Create a custom palette based on user input
                                  low_color <- input$negColor
                                  mid_color <- input$zeroColor
                                  high_color <- input$posColor
                                  
                                  # Handle cases where the number of unique values in the matrix is less than 3.
                                  unique_values <- length(unique(as.vector(numeric_data)))
                                  if (unique_values < 3) {
                                    if (unique_values == 1) {
                                      colorRampPalette(c(mid_color))(100)
                                    } else {
                                      colorRampPalette(c(low_color, high_color))(100)
                                    }
                                  } else {
                                    colorRampPalette(c(low_color, mid_color, high_color))(100)
                                  }
                                },
                                colorRampPalette(c("blue", "white", "red"))(100) # Default
        )
        
        output$heatmapPlot <- renderPlot({
          pheatmap(numeric_data, scale = "row",
                   cluster_rows = input$clusterRows,
                   cluster_cols = input$clusterCols,
                   color = color_palette
          )
        }, width = function() input$plotWidth, height = function() input$plotHeight)
        
        return(NULL)
        
      }, error = function(e) {
        return(paste("An error occurred:", e$message))
      })
    })
  })
  
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste("heatmap-", Sys.Date(), ".pdf", sep="")
    },
    content = function(file) {
      if (!is.null(heatmap_data())) {
        
        # Color palette for download
        color_palette_download <- switch(input$colorPalette,
                                         "RdBu" = colorRampPalette(c("blue", "white", "red"))(100),
                                         "viridis" = viridis::viridis(100),
                                         "magma" = viridis::magma(100),
                                         "plasma" = viridis::plasma(100),
                                         "inferno" = viridis::inferno(100),
                                         "Custom" = {
                                           # Create a custom palette based on user input
                                           low_color <- input$negColor
                                           mid_color <- input$zeroColor
                                           high_color <- input$posColor
                                           
                                           # Handle cases where the number of unique values in the matrix is less than 3.
                                           unique_values <- length(unique(as.vector(heatmap_data())))
                                           if (unique_values < 3) {
                                             if (unique_values == 1) {
                                               colorRampPalette(c(mid_color))(100)
                                             } else {
                                               colorRampPalette(c(low_color, high_color))(100)
                                             }
                                           } else {
                                             colorRampPalette(c(low_color, mid_color, high_color))(100)
                                           }
                                         },
                                         colorRampPalette(c("blue", "white", "red"))(100) # Default
        )
        
        pdf(file, width = input$plotWidth/72, height = input$plotHeight/72) #72 points per inch
        pheatmap(heatmap_data(), scale = "row",
                 cluster_rows = input$clusterRows,
                 cluster_cols = input$clusterCols,
                 color = color_palette_download)
        dev.off()
      }
    }
  )
}

shinyApp(ui = ui, server = server)
