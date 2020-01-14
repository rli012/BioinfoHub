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

}


shinyApp(
  ui = ui,
  server = server
)
