library(shiny)
library(shinydashboard)
library(DT)
library(Biostrings)
library(stringr)
library(dplyr)
library(readr)

# Helper functions for primer design
calculate_tm <- function(sequence) {
  # Basic Tm calculation using nearest neighbor method approximation
  # This is a simplified version - in practice, you'd use more sophisticated algorithms
  seq_upper <- toupper(sequence)
  gc_content <- (str_count(seq_upper, "G") + str_count(seq_upper, "C")) / nchar(seq_upper)
  length <- nchar(seq_upper)
  
  # Wallace rule for short primers, adjusted for GC content
  if (length < 14) {
    tm <- (str_count(seq_upper, "A") + str_count(seq_upper, "T")) * 2 + 
      (str_count(seq_upper, "G") + str_count(seq_upper, "C")) * 4
  } else {
    # More accurate formula for longer primers
    tm <- 64.9 + 41 * (gc_content) - (675 / length)
  }
  
  return(round(tm, 1))
}

reverse_complement <- function(sequence) {
  # Convert to DNAString and get reverse complement
  dna_seq <- DNAString(toupper(sequence))
  return(as.character(reverseComplement(dna_seq)))
}

design_sitedirected_primers <- function(template, mutation_pos, new_base, primer_length = 25) {
  template_upper <- toupper(template)
  
  # Ensure mutation position is valid
  if (mutation_pos < 1 || mutation_pos > nchar(template_upper)) {
    return(NULL)
  }
  
  # Design forward primer centered around mutation
  start_pos <- max(1, mutation_pos - floor(primer_length/2))
  end_pos <- min(nchar(template_upper), start_pos + primer_length - 1)
  start_pos <- max(1, end_pos - primer_length + 1)
  
  # Create mutated sequence
  mutated_template <- paste0(
    substr(template_upper, 1, mutation_pos - 1),
    toupper(new_base),
    substr(template_upper, mutation_pos + 1, nchar(template_upper))
  )
  
  forward_primer <- substr(mutated_template, start_pos, end_pos)
  reverse_primer <- reverse_complement(forward_primer)
  
  return(list(
    forward = forward_primer,
    reverse = reverse_primer,
    mutation_pos = mutation_pos,
    original_base = substr(template_upper, mutation_pos, mutation_pos),
    new_base = toupper(new_base)
  ))
}

design_goldengate_primers <- function(sequence, restriction_site = "GGTCTC", overhang = "AATG", primer_length = 25) {
  sequence_upper <- toupper(sequence)
  
  # Forward primer with restriction site and overhang
  core_length <- primer_length - nchar(restriction_site) - nchar(overhang)
  forward_core <- substr(sequence_upper, 1, core_length)
  forward_primer <- paste0(restriction_site, overhang, forward_core)
  
  # Reverse primer with restriction site and overhang
  reverse_core_start <- nchar(sequence_upper) - core_length + 1
  reverse_core <- substr(sequence_upper, reverse_core_start, nchar(sequence_upper))
  reverse_primer <- paste0(restriction_site, overhang, reverse_complement(reverse_core))
  
  return(list(
    forward = forward_primer,
    reverse = reverse_primer,
    restriction_site = restriction_site,
    overhang = overhang
  ))
}

design_cloning_primers <- function(sequence, forward_tail = "GAATTC", reverse_tail = "AAGCTT", primer_length = 20) {
  sequence_upper <- toupper(sequence)
  
  # Forward primer with tail
  forward_core <- substr(sequence_upper, 1, primer_length)
  forward_primer <- paste0(forward_tail, forward_core)
  
  # Reverse primer with tail
  reverse_core_start <- nchar(sequence_upper) - primer_length + 1
  reverse_core <- substr(sequence_upper, reverse_core_start, nchar(sequence_upper))
  reverse_primer <- paste0(reverse_tail, reverse_complement(reverse_core))
  
  return(list(
    forward = forward_primer,
    reverse = reverse_primer,
    forward_tail = forward_tail,
    reverse_tail = reverse_tail
  ))
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Multiplexed DNA Primer Designer"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Site-Directed Mutagenesis", tabName = "sitedirected", icon = icon("dna")),
      menuItem("Golden Gate Assembly", tabName = "goldengate", icon = icon("puzzle-piece")),
      menuItem("DNA Cloning", tabName = "cloning", icon = icon("clone")),
      menuItem("Primer Results", tabName = "results", icon = icon("table")),
      menuItem("Export Data", tabName = "export", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Site-directed mutagenesis tab
      tabItem(
        tabName = "sitedirected",
        fluidRow(
          box(
            title = "Site-Directed Mutagenesis Primer Design", status = "primary", width = 12,
            textAreaInput("sd_template", "Template Sequence:", 
                          placeholder = "Enter your template DNA sequence here...",
                          height = "100px", width = "100%"),
            fluidRow(
              column(4, numericInput("mutation_pos", "Mutation Position:", value = 1, min = 1)),
              column(4, selectInput("new_base", "New Base:", choices = c("A", "T", "G", "C"))),
              column(4, numericInput("sd_primer_length", "Primer Length:", value = 25, min = 15, max = 50))
            ),
            textInput("sd_primer_name", "Primer Name Prefix:", value = "SDM"),
            actionButton("design_sd", "Design Site-Directed Primers", class = "btn-primary")
          )
        )
      ),
      
      # Golden Gate tab
      tabItem(
        tabName = "goldengate",
        fluidRow(
          box(
            title = "Golden Gate Assembly Primer Design", status = "warning", width = 12,
            textAreaInput("gg_sequence", "Target Sequence:", 
                          placeholder = "Enter the sequence to be amplified...",
                          height = "100px", width = "100%"),
            fluidRow(
              column(3, textInput("restriction_site", "Restriction Site:", value = "GGTCTC")),
              column(3, textInput("gg_overhang", "Overhang Sequence:", value = "AATG")),
              column(3, numericInput("gg_primer_length", "Core Primer Length:", value = 25, min = 15, max = 40)),
              column(3, textInput("gg_primer_name", "Primer Name Prefix:", value = "GG"))
            ),
            actionButton("design_gg", "Design Golden Gate Primers", class = "btn-warning")
          )
        )
      ),
      
      # DNA Cloning tab
      tabItem(
        tabName = "cloning",
        fluidRow(
          box(
            title = "DNA Cloning Primer Design", status = "success", width = 12,
            textAreaInput("cloning_sequence", "Target Sequence:", 
                          placeholder = "Enter the sequence to be cloned...",
                          height = "100px", width = "100%"),
            fluidRow(
              column(3, textInput("forward_tail", "Forward Primer Tail:", value = "GAATTC")),
              column(3, textInput("reverse_tail", "Reverse Primer Tail:", value = "AAGCTT")),
              column(3, numericInput("cloning_primer_length", "Core Primer Length:", value = 20, min = 15, max = 35)),
              column(3, textInput("cloning_primer_name", "Primer Name Prefix:", value = "CLONE"))
            ),
            actionButton("design_cloning", "Design Cloning Primers", class = "btn-success")
          )
        )
      ),
      
      # Results tab
      tabItem(
        tabName = "results",
        fluidRow(
          box(
            title = "Primer Design Results", status = "info", width = 12,
            DT::dataTableOutput("primer_table"),
            br(),
            h4("Primer Analysis Summary"),
            verbatimTextOutput("primer_summary")
          )
        )
      ),
      
      # Export tab
      tabItem(
        tabName = "export",
        fluidRow(
          box(
            title = "Export Primer Data", status = "info", width = 12,
            h4("Download Options"),
            br(),
            downloadButton("download_csv", "Download as CSV", class = "btn-primary"),
            br(), br(),
            downloadButton("download_fasta", "Download as FASTA", class = "btn-success"),
            br(), br(),
            h4("Export Preview"),
            DT::dataTableOutput("export_preview")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive value to store all primers
  primer_data <- reactiveVal(data.frame(
    Primer_Name = character(0),
    Sequence = character(0),
    Type = character(0),
    Direction = character(0),
    Tm = numeric(0),
    Length = numeric(0),
    GC_Content = numeric(0),
    Paired_With = character(0),
    Notes = character(0),
    stringsAsFactors = FALSE
  ))
  
  # Site-directed mutagenesis
  observeEvent(input$design_sd, {
    req(input$sd_template, input$mutation_pos, input$new_base)
    
    # Validate sequence
    template_clean <- gsub("[^ATGCatgc]", "", input$sd_template)
    if (nchar(template_clean) == 0) {
      showNotification("Invalid template sequence. Please enter a valid DNA sequence.", type = "error")
      return()
    }
    
    # Update mutation position limits
    updateNumericInput(session, "mutation_pos", max = nchar(template_clean))
    
    if (input$mutation_pos > nchar(template_clean)) {
      showNotification("Mutation position exceeds sequence length.", type = "error")
      return()
    }
    
    primers <- design_sitedirected_primers(template_clean, input$mutation_pos, 
                                           input$new_base, input$sd_primer_length)
    
    if (!is.null(primers)) {
      # Calculate additional metrics
      forward_tm <- calculate_tm(primers$forward)
      reverse_tm <- calculate_tm(primers$reverse)
      forward_gc <- round((str_count(primers$forward, "[GC]") / nchar(primers$forward)) * 100, 1)
      reverse_gc <- round((str_count(primers$reverse, "[GC]") / nchar(primers$reverse)) * 100, 1)
      
      # Create new primer entries
      new_primers <- data.frame(
        Primer_Name = c(paste0(input$sd_primer_name, "_F"), paste0(input$sd_primer_name, "_R")),
        Sequence = c(primers$forward, primers$reverse),
        Type = c("Site-Directed Mutagenesis", "Site-Directed Mutagenesis"),
        Direction = c("Forward", "Reverse"),
        Tm = c(forward_tm, reverse_tm),
        Length = c(nchar(primers$forward), nchar(primers$reverse)),
        GC_Content = c(forward_gc, reverse_gc),
        Paired_With = c(paste0(input$sd_primer_name, "_R"), paste0(input$sd_primer_name, "_F")),
        Notes = c(
          paste0("Mutation: ", primers$original_base, primers$mutation_pos, primers$new_base),
          paste0("Mutation: ", primers$original_base, primers$mutation_pos, primers$new_base)
        ),
        stringsAsFactors = FALSE
      )
      
      # Add to existing data
      current_data <- primer_data()
      primer_data(rbind(current_data, new_primers))
      
      showNotification("Site-directed mutagenesis primers designed successfully!", type = "success")
    } else {
      showNotification("Failed to design primers. Check your inputs.", type = "error")
    }
  })
  
  # Golden Gate assembly
  observeEvent(input$design_gg, {
    req(input$gg_sequence)
    
    sequence_clean <- gsub("[^ATGCatgc]", "", input$gg_sequence)
    if (nchar(sequence_clean) == 0) {
      showNotification("Invalid sequence. Please enter a valid DNA sequence.", type = "error")
      return()
    }
    
    primers <- design_goldengate_primers(sequence_clean, input$restriction_site, 
                                         input$gg_overhang, input$gg_primer_length)
    
    if (!is.null(primers)) {
      forward_tm <- calculate_tm(primers$forward)
      reverse_tm <- calculate_tm(primers$reverse)
      forward_gc <- round((str_count(primers$forward, "[GC]") / nchar(primers$forward)) * 100, 1)
      reverse_gc <- round((str_count(primers$reverse, "[GC]") / nchar(primers$reverse)) * 100, 1)
      
      new_primers <- data.frame(
        Primer_Name = c(paste0(input$gg_primer_name, "_F"), paste0(input$gg_primer_name, "_R")),
        Sequence = c(primers$forward, primers$reverse),
        Type = c("Golden Gate Assembly", "Golden Gate Assembly"),
        Direction = c("Forward", "Reverse"),
        Tm = c(forward_tm, reverse_tm),
        Length = c(nchar(primers$forward), nchar(primers$reverse)),
        GC_Content = c(forward_gc, reverse_gc),
        Paired_With = c(paste0(input$gg_primer_name, "_R"), paste0(input$gg_primer_name, "_F")),
        Notes = c(
          paste0("Restriction site: ", primers$restriction_site, ", Overhang: ", primers$overhang),
          paste0("Restriction site: ", primers$restriction_site, ", Overhang: ", primers$overhang)
        ),
        stringsAsFactors = FALSE
      )
      
      current_data <- primer_data()
      primer_data(rbind(current_data, new_primers))
      
      showNotification("Golden Gate assembly primers designed successfully!", type = "success")
    }
  })
  
  # DNA Cloning
  observeEvent(input$design_cloning, {
    req(input$cloning_sequence)
    
    sequence_clean <- gsub("[^ATGCatgc]", "", input$cloning_sequence)
    if (nchar(sequence_clean) == 0) {
      showNotification("Invalid sequence. Please enter a valid DNA sequence.", type = "error")
      return()
    }
    
    primers <- design_cloning_primers(sequence_clean, input$forward_tail, 
                                      input$reverse_tail, input$cloning_primer_length)
    
    if (!is.null(primers)) {
      forward_tm <- calculate_tm(primers$forward)
      reverse_tm <- calculate_tm(primers$reverse)
      forward_gc <- round((str_count(primers$forward, "[GC]") / nchar(primers$forward)) * 100, 1)
      reverse_gc <- round((str_count(primers$reverse, "[GC]") / nchar(primers$reverse)) * 100, 1)
      
      new_primers <- data.frame(
        Primer_Name = c(paste0(input$cloning_primer_name, "_F"), paste0(input$cloning_primer_name, "_R")),
        Sequence = c(primers$forward, primers$reverse),
        Type = c("DNA Cloning", "DNA Cloning"),
        Direction = c("Forward", "Reverse"),
        Tm = c(forward_tm, reverse_tm),
        Length = c(nchar(primers$forward), nchar(primers$reverse)),
        GC_Content = c(forward_gc, reverse_gc),
        Paired_With = c(paste0(input$cloning_primer_name, "_R"), paste0(input$cloning_primer_name, "_F")),
        Notes = c(
          paste0("Forward tail: ", primers$forward_tail),
          paste0("Reverse tail: ", primers$reverse_tail)
        ),
        stringsAsFactors = FALSE
      )
      
      current_data <- primer_data()
      primer_data(rbind(current_data, new_primers))
      
      showNotification("DNA cloning primers designed successfully!", type = "success")
    }
  })
  
  # Display primer table
  output$primer_table <- DT::renderDataTable({
    primer_data()
  }, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  
  # Primer summary
  output$primer_summary <- renderText({
    data <- primer_data()
    if (nrow(data) == 0) {
      return("No primers designed yet.")
    }
    
    summary_text <- paste0(
      "Total primers: ", nrow(data), "\n",
      "Primer types: ", paste(unique(data$Type), collapse = ", "), "\n",
      "Average Tm: ", round(mean(data$Tm, na.rm = TRUE), 1), "°C\n",
      "Tm range: ", round(min(data$Tm, na.rm = TRUE), 1), "°C - ", round(max(data$Tm, na.rm = TRUE), 1), "°C\n",
      "Average GC content: ", round(mean(data$GC_Content, na.rm = TRUE), 1), "%\n",
      "Length range: ", min(data$Length, na.rm = TRUE), "-", max(data$Length, na.rm = TRUE), " bp"
    )
    
    return(summary_text)
  })
  
  # Export preview
  output$export_preview <- DT::renderDataTable({
    primer_data()
  }, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  
  # Download handlers
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("primer_designs_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(primer_data(), file)
    }
  )
  
  output$download_fasta <- downloadHandler(
    filename = function() {
      paste0("primer_designs_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      data <- primer_data()
      if (nrow(data) > 0) {
        fasta_content <- paste0(
          ">", data$Primer_Name, " | ", data$Type, " | Tm:", data$Tm, "C | GC:", data$GC_Content, "%\n",
          data$Sequence,
          collapse = "\n"
        )
        writeLines(fasta_content, file)
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)