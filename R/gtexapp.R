#' run the app with GTEx lung and whole blood
#' @examples
#' oldask = options()$ask
#' if (interactive()) {
#'   options(ask=FALSE)
#'   gtexapp()
#' }
#' options(ask=oldask)
#' @export
gtexapp = function() {
 data("gloc_hg19", package="GGIpack2")
 data("newensg", package="GGIpack2")
 ensg = newensg[order(names(newensg))] # preserve old name

 gsyms = names(ensg)
 names(gsyms) = as.character(ensg)
 e2g = function(e) {
   canswap = which(e %in% names(gsyms))
   e[canswap] = gsyms[canswap]
   e
 }
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
  #  prepare output components -- a list of renderDataTable results
    z = lapply(names(resl), function(x) {
      output[[x]] = 
       DT::renderDataTable({
         if (input$focus == "chr")
           dat =resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
                  dplyr::arrange(score) |>
                  head(input$nrecs) |> as.data.frame() 
         else if (input$focus == "gene") 
           dat =resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
                  as.data.frame() 
         else if (input$focus == "rsid") 
           dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
                  as.data.frame() 
         dat$sym = e2g(dat$molecular_trait_id)
         dat
       })
      }
      )
    output$theplot = plotly::renderPlotly({ 
#      par(mfrow=c(length(input$respicks), 1))
      fulldat = NULL
      for (x in input$respicks) {
         if (input$focus == "chr")
           dat = resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |> dplyr::arrange(score) |>
                  head(input$nrecs) |> as.data.frame() 
         else if (input$focus == "gene") 
           dat = resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
                  as.data.frame() 
         else if (input$focus == "rsid") 
           dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
                  as.data.frame() 
         #plot(dat$start, -log10(dat$score), main=x)
         if (is.null(fulldat)) fulldat = dat
         else fulldat = rbind(fulldat, dat)
         }
         fulldat$sym = e2g(fulldat$molecular_trait_id)
         theplot = ggplot2::ggplot(data=fulldat, 
                    ggplot2::aes(x=start, y=-log10(score),
              text=sym)) + ggplot2::geom_point() 
         ntiss = length(unique(fulldat$tissue))
         if (ntiss > 1) theplot = theplot + ggplot2::facet_grid(tissue~.)
         plotly::ggplotly(theplot)
      })
  # communicate selected components to UI
    output$all = renderUI({
     o = lapply(c(input$respicks, "viz"), function(x) {
              if (x != "viz") tabPanel(x, DT::dataTableOutput(x))
              else tabPanel("viz", plotly::plotlyOutput("theplot"))
              })
     o = c(o, list(tabPanel("about", helpText("This is an app demonstrating organization of eQTL data with parquet and duckdb"))))
     names(o) = NULL
     do.call(tabsetPanel, o)
    })
   }
  
  runApp(list(ui=ui, server=server))
}
