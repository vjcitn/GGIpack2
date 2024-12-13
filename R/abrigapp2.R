#' run the app with two abrig resources
#' @examples
#' oldask = options()$ask
#' if (interactive()) {
#'   options(ask=FALSE)
#'   abrigapp2()
#' }
#' options(ask=oldask)
#' @export
abrigapp2 = function() {
 data("gloc_hg19", package="GGIpack2")
 utils::data("geneNames", package = "GGIpack2")

# set up data resources
 con = DBI::dbConnect(duckdb::duckdb())

pfiles =
c("/udd/remcr/abrig/ca_ba.parquet", "/udd/remcr/abrig/ca_brep.parquet",
"/udd/remcr/abrig/cc_s.parquet", "/udd/remcr/abrig/cc_unstim.parquet",
"/udd/remcr/abrig/ch_am.parquet", "/udd/remcr/abrig/comb.parquet")

names(pfiles) = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")


 BALres = ABRIGresource(con, "BAL", pfiles = pfiles)
 BEBres = ABRIGresource(con, "BroncEpiBrush", pfiles = pfiles)
 CD4Stim = ABRIGresource(con, "CD4Stim", pfiles = pfiles)
 CD4Unstim = ABRIGresource(con, "CD4Unstim", pfiles = pfiles)
 AlvMacphage = ABRIGresource(con, "AlvMacphage", pfiles = pfiles)
 PaxRNA = ABRIGresource(con, "PaxRNA", pfiles = pfiles)

 resl = list(BAL=BALres, BroncEpiBr=BEBres,
     CD4Stim=CD4Stim, CD4Unstim=CD4Unstim, AlvMacphage=AlvMacphage, PaxRNA=PaxRNA)

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
    updateSelectizeInput(session, 'gene', choices = geneNames, server = TRUE)
    observeEvent(input$stop, {
      DBI::dbDisconnect(con)
      stopApp()
      })
  #  prepare output components -- a list of renderDataTable results
    z = lapply(names(resl), function(x) {
      output[[x]] = 
       DT::renderDataTable({
         if (input$focus == "chr")
           dat = resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
                  dplyr::arrange(score) |>
                  head(input$nrecs) |> as.data.frame() 
         else if (input$focus == "gene") 
           dat = resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
                  as.data.frame() 
         else if (input$focus == "rsid") 
           dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
                  as.data.frame() 
       })
      }
      )
    output$theplot = plotly::renderPlotly({ 
      thedat = NULL
      for (x in input$respicks) {
         if (input$focus == "chr")
           dat = resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
                  dplyr::arrange(score) |>
                  head(input$nrecs) |> as.data.frame() 
         else if (input$focus == "gene") 
           dat = resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
                  as.data.frame() 
         else if (input$focus == "rsid") 
           dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
                  as.data.frame() 
         dat$tissue = x
         if (is.null(thedat)) thedat = dat
         else thedat = rbind(thedat, dat)
         }
         plotly::ggplotly(ggplot2::ggplot(data=thedat, 
                  ggplot2::aes(x=start, y=-log10(score))) +
              ggplot2::geom_point() + ggplot2::facet_grid(tissue~.))
      })
  # communicate selected components to UI
    output$all = renderUI({
     o = lapply(c(input$respicks, "viz"), function(x) {
              if (x != "viz") tabPanel(x, DT::dataTableOutput(x))
              else tabPanel("viz", plotly::plotlyOutput("theplot"))
              })
     o = c(o, list(tabPanel("about", helpText(sprintf("(GGIpack2 %s) This is an app demonstrating organization of eQTL data with parquet and duckdb", packageVersion("GGIpack2"))))))
     names(o) = NULL

     do.call(tabsetPanel, o)
    })
   }
  
  runApp(list(ui=ui, server=server))
}
