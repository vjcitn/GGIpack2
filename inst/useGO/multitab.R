library(shiny)

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   numericInput("ntabs", "ntabs", min=1, max=5, value=1)
  ),
  mainPanel(
   uiOutput("mult")
  )
 )
)

server = function(input, output) {
 z = vector("list", 3)
 lapply(2:1, function(i) {
  dat = iris[1:(i+1),]
  output[[paste0("dat", i)]] = DT::renderDataTable( dat )
  })
 output$dat3 = DT::renderDataTable({
  DT::datatable(mtcars[1:4,])
  })
 output$mult = renderUI({
  #tl = lapply(seq_len(input$ntabs), function(x) tabPanel(x, DT::dataTableOutput(paste0("dat", x))))
  tl = lapply(seq_len(input$ntabs), function(x) tabPanel(x, helpText(x), DT::dataTableOutput(paste0("dat", x))))
  do.call(tabsetPanel, tl)
  })
}

runApp(list(ui=ui, server=server))
