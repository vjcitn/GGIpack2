library(GGIpack)
library(shiny)
data("gloc_hg19", package="GGIpack")
load("ensg.rda")
ensg = ensg[order(names(ensg))]

Sys.setenv("GGI_PARQUET_FOLDER"="/udd/stvjc")

# set up data resources
con = DBI::dbConnect(duckdb::duckdb())
lungpa = try(ggi_gtex_cache("lungpl05.parquet"))
if (inherits(lungpa, "try-error") | nchar(lungpa)==0)
   lungpa = file.path(Sys.getenv("GGI_PARQUET_FOLDER"), "lungpl05.parquet")
lungres = GTExresource(con, tisstag="lung", pfile=lungpa)
wbpa = try(ggi_gtex_cache("wholeblpl05.parquet"))
if (inherits(wbpa, "try-error") | nchar(wbpa)==0)
   wbpa = file.path(Sys.getenv("GGI_PARQUET_FOLDER"), "wholeblpl05.parquet")
whblres = GTExresource(con, tisstag="wholebl", pfile=wbpa)
resl = list(lung=lungres, wholebl=whblres)
fvec = names(resl[[1]]@tbl |> head(2) |> as.data.frame())

#selectizeInput('foo', choices = NULL, ...)
#
## in server
#server <- function(input, output, session) {
#  updateSelectizeInput(session, 'foo', choices = data, server = TRUE)
#}

# simple UI based on selections

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("using gtex eqtl data"),
   checkboxGroupInput("respicks", "resources",
        choices=names(resl), selected=names(resl)[1]),
   numericInput("nrecs", "nrecs", min=5, max=100000, value=10000), 
   radioButtons("focus", "focus", choices=c("chr", "gene", "rsid")),
   conditionalPanel(
    condition = "input.focus == 'chr'",
    radioButtons("chr", "chr", choices=1:22, selected=1, inline=TRUE)
    ),
   conditionalPanel(
    condition = "input.focus == 'gene'",
    selectInput("gene", "gene", choices=NULL)
    ),
   conditionalPanel(
    condition = "input.focus == 'rsid'",
    textInput("snp", "snp")
    ),
   actionButton("stop", "stop app"),
   width=2
   ),
  mainPanel(
   uiOutput("all")
    )
   )
  )

server = function(input, output, session) {
  updateSelectizeInput(session, 'gene', choices = ensg, server = TRUE)
  observeEvent(input$stop, {
    DBI::dbDisconnect(con)
    stopApp()
    })
# prepare output components
# for now very simple processing of tables, later, perform filtering
# based on gene selection, must be reactive (would not use 'head()' but 
# a filter based on locations
#  build_data = reactive({
#    if (input$focus == "chr") {
#    lapply(names(resl), function(x) {
#     resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr)))
#     })
#    }
#    else if (input$focus == "gene") {
#    lapply(names(resl), function(x) {
#     resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene)))
#     })
#    }
#  })
  z = lapply(names(resl), function(x) {
  output[[x]] = 
     DT::renderDataTable({
       if (input$focus == "chr")
         resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
                dplyr::arrange(score) |>
                head(input$nrecs) |> as.data.frame() 
       else if (input$focus == "gene") 
         resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
                as.data.frame() 
       else if (input$focus == "rsid") 
         dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
                as.data.frame() 
     })
    }
    )
  output$theplot = renderPlot({ 
    par(mfrow=c(length(input$respicks), 1))
    for (x in input$respicks) {
       if (input$focus == "chr")
         dat = resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
                head(input$nrecs) |> as.data.frame() 
       else if (input$focus == "gene") 
         dat = resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
                as.data.frame() 
       else if (input$focus == "rsid") 
         dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
                as.data.frame() 
       plot(dat$start, -log10(dat$score), main=x)
       }
    }, height=800L)
# communicate selected components to UI
  output$all = renderUI({
   o = lapply(c(input$respicks, "viz"), function(x) {
            if (x != "viz") tabPanel(x, DT::dataTableOutput(x))
            else tabPanel("viz", shiny::plotOutput("theplot"))
            })
   do.call(tabsetPanel, o)
  })
 }

runApp(list(ui=ui, server=server))
