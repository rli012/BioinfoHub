header=dashboardHeader(title = 'PCa Transcriptomes')

sidebar=dashboardSidebar(
  sidebarMenu(
    style = 'position:fixed; overflow: visible',
    
    tags$div(tags$h4("Gene-Level"),
             style = 'color: white'),
    menuItem("Differential Expression", tabName = 'tab_degene', icon = icon("dna"),
             menuSubItem("Boxplot" , tabName = "tab_boxplot", icon = icon("bar-chart")),
             menuSubItem("Heatmap" , tabName = "heatmap", icon = icon("scatter-chart"))
    ),
    menuItem("Survival Analysis" , tabName = "tab_kmplot", icon = icon("pencil"))
    
  )
)

body=dashboardBody(
  
  tabItems(
    tabItem(tabName="tab_boxplot", tab_boxplot),
    tabItem(tabName="tab_kmplot",tab_kmplot),
    tabItem(tabName="tab_dataset",tab_dataset)
  )
  
)

ui <- dashboardPage(title='PCa Transcriptomes', skin = 'green', header, sidebar, body)

server <- function(input, output, session) {
  output$plot <- renderPlot({})
  outpout$text <- renderText({})
  output$table <- DT::renderDataTable({})
}


shinyApp(
  ui = ui,
  server = server
)



##################################################################################################

gene.default <- 'ENSG00000142515'

gene.expression <- selectizeInput(inputId = "gene.expression", label=h4(strong('Gene')), choices = NULL, selected = gene.default, 
                                  multiple = FALSE, width = 300,
                                  options = list(placeholder = 'Select a gene',
                                                 server = TRUE, selectOnTab=TRUE,
                                                 searchField = c('external_gene_name', 'alias_symbol', 'description', 'ensembl_id', 'entrez_id'),
                                                 labelField = "external_gene_name",
                                                 valueField = "ensembl_id",
                                                 maxOptions = 5,
                                                 render = I("{option: function(item, escape) {
                                  var gene = '<div>' + '<strong>' + escape(item.external_gene_name) + '</strong>:' + '<ul>';
                                  gene = gene + '<li>' + item.alias_symbol + '</li>';
                                  gene = gene + '<li>' + item.description + '</li>';
                                  gene = gene + '<li>' + 'Entrez: ' + item.entrez_id + '</li>';
                                  gene = gene + '<li>' + 'Ensembl: ' + item.ensembl_id + '</li>' + '</ul>' + '</div>';
                                  return gene
                                  }
                                             }")
                                  ))

