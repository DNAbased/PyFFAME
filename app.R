library(shiny)

ui <- fluidPage(
  
  titlePanel("Hello!"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      actionButton("create_file", 
                   "Create config file", 
                   class = "btn-primary"),
      tags$hr(),
      downloadButton("download_file", 
                   "Download config file", 
                   class = "btn-primary")
      
    ),
    
    mainPanel(
      
      textInput("in_vcf", "VCF-like input file"),
      textInput("in_sequence", "Sequence input file"),
      textInput("in_barcode", "Barcode input file"),
      radioButtons("in_barcode_type", "Barcode input file type",
                  choices=c(".json", ".txt")),
      textInput("db_genome", "Genome input"),
      textInput("de_order", "Design order",
                placeholder="abcde"),
      textInput("de_seq_1", "Added sequence 1"),
      textInput("de_seq_2", "Added sequence 2"),
      textInput("de_seq_3", "Added sequence 3"),
      numericInput("set_feature_size", "Feature size per side",
                   value=85),
      radioButtons("set_all_features", "Create all allelic features?",
                   choices=c("Yes", "No")),
      numericInput("set_indel_max_length", "Discard indels larger than this",
                   value=10),
      radioButtons("set_indel_features", "Only create full-length features",
                   choices=c("Yes", "No"), selected="No"),
      numericInput("set_barcodes_per_feature", "Barcodes added to each feature",
                   value=100),
      radioButtons("set_rev_comp", "Create reverse complementary versions of all features",
                   choices=c("Yes", "No"), selected="No"),
      textInput("enz_file_processed", "Enzyme file (processed)"),
      textInput("enz_file", "Enzyme file (NEB REBASE format)"),
      textInput("enz_used", "All enzymes used",
                placeholder="EcoRI,SbfI"),
      textInput("enz_sites", "OR: all relevant cut sites (DOUBLE QUOTES)",
                placeholder="['GAATTC', 'CCTGCAGG']"),
      numericInput("enz_cumul_cuts", "Total restriction sites expected",
                   value=2),
      numericInput("enz_cumul_cuts_bc", "Total restriction sites expected (in barcode plus the sequence before and after it)",
                   value=2),
      radioButtons("out_format", "Output format",
                   choices=c(".json", ".tsv")),
      textInput("out_output", "Output file path")
      
    )
  )
)


server <- function(input, output) {
  
  observeEvent(input$create_file, {
    #write("test = 'bla'\ntest2 = 'bla2'", "../z_testing.txt")
    in_vcf = ifelse(input$in_vcf == "", "None", paste0("'", input$in_vcf, "'"))
    in_sequence = ifelse(input$in_sequence == "", "None", paste0("'", input$in_sequence, "'"))
    in_barcode = paste0("'", input$in_barcode, "'")
    in_barcode_type = paste0("'", input$in_barcode_type, "'")
    db_genome = paste0("'", input$db_genome, "'")
    de_order = ifelse(input$de_order == "", paste0("'", "abcde", "'"), paste0("'", input$de_order, "'"))
    de_seq_1 = paste0("'", input$de_seq_1, "'")
    de_seq_2 = paste0("'", input$de_seq_2, "'")
    de_seq_3 = paste0("'", input$de_seq_3, "'")
    set_feature_size = input$set_feature_size
    set_all_features = ifelse(input$set_all_features == "Yes", 1, 0)
    set_indel_max_length = input$set_indel_max_length
    set_indel_features = ifelse(input$set_indel_features == "Yes", 1, 0)
    set_barcodes_per_feature = input$set_barcodes_per_feature
    set_rev_comp = ifelse(input$set_rev_comp == "Yes", 1, 0)
    enz_file_processed = ifelse(input$enz_file_processed == "", "None", paste0("'", input$enz_file_processed, "'"))
    enz_file = ifelse(input$enz_file == "", "None", paste0("'", input$enz_file, "'"))
    enz_used = ifelse(input$enz_used == "", "None", paste0("'", input$enz_used, "'"))
    enz_sites = ifelse(input$enz_sites == "", "None", paste0("'", input$enz_sites, "'"))
    enz_cumul_cuts = input$enz_cumul_cuts
    enz_cumul_cuts_bc = input$enz_cumul_cuts_bc
    out_format = paste0("'", input$out_format, "'")
    out_output = paste0("'", input$out_output, "'")
    
    final_output = paste0("in_vcf = ", in_vcf,
                          "\nin_sequence = ", in_sequence,
                          "\nin_barcode = ", in_barcode,
                          "\nin_barcode_type = ", in_barcode_type,
                          "\ndb_genome = ", db_genome,
                          "\nde_order = ", de_order,
                          "\nde_seq_1 = ", de_seq_1,
                          "\nde_seq_2 = ", de_seq_2,
                          "\nde_seq_3 = ", de_seq_3,
                          "\nset_feature_size = ", set_feature_size,
                          "\nset_all_features = ", set_all_features,
                          "\nset_indel_max_length = ", set_indel_max_length,
                          "\nset_indel_features = ", set_indel_features,
                          "\nset_barcodes_per_feature = ", set_barcodes_per_feature,
                          "\nset_rev_comp = ", set_rev_comp,
                          "\nenz_file_processed = ", enz_file_processed,
                          "\nenz_file = ", enz_file,
                          "\nenz_used = ", enz_used,
                          "\nenz_sites = ", enz_sites,
                          "\nenz_cumul_cuts = ", enz_cumul_cuts,
                          "\nenz_cumul_cuts_bc = ", enz_cumul_cuts_bc,
                          "\nout_format = ", out_format,
                          "\nout_output = ", out_output)
    write(final_output, "../config.py")
  })
  output$download_file = downloadHandler(filename=function(){paste0("config.py")},
                                         content=function(filename){file.copy("../config.py", filename)})
  
}

shinyApp(ui, server)

