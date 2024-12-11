#' demo app 2
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @rawNamespace import(GenomicRanges, except=c(intersect, union, setdiff))
#' @import igvShiny
#' @import DT
#' @param con a DBI connection
#' @param genelocs a GRanges instance with gene addresses
#' @note Very specialized, just has a few genes, uses specific
#' field from genelocs argument.
#' @examples
#' if (interactive()){
#' utils::data("gloc_hg19", package = "GGIpack2")
#' con = DBI::dbConnect(duckdb::duckdb())
#' abrigapp(con, gloc_hg19)
#' }
#' @export
abrigapp = function(con, genelocs) {
 pfiles <<- ABRIGparquet_paths()
 bpPadding = 30000
 utils::data("geneNames", package = "GGIpack2")
 patquetTableLoc <-system.file("extdata","parquetDataTable.csv", package = "GGIpack2" )
 patquetTable <- utils::read.csv(patquetTableLoc)
 genomeVersion <- "hg19"
 
 #-----------------------------------------------------------------------------------
 if(!dir.exists("tracks"))
   dir.create("tracks")
 shiny::addResourcePath("tracks", "tracks")
 #-------------------------------------------------------------------------------------
 
 
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("GGI demo"),
    selectizeInput("gene", "gene", geneNames) #,
    #actionButton("zoomButton", label = "Zoom")
    ), 
   mainPanel(
    uiOutput("alltabs")
    ) #mainPanel
  ) #sidebarLayout
 )#fluidPage
 
 
 server = function(input, output, session) {
  updateSelectizeInput(session, "gene", choices =sort(geneNames),server = TRUE)
   
  output$stuff = renderPrint({
   mygene = input$gene
   mytiss = input$tiss
   newres = ABRIGresource(con, input$tiss, pfiles = pfiles)
   kk <- filterByRange(newres, genelocs, mygene, ggr_field="gene_name")
   kk@tbl
  }) #outputStuff

#####################################################################
## functions for the datatables
####################################################################
  
  allCelltypes = reactive({
    req(c(input$gene, con, pfiles, genelocs ))
    mygene = input$gene
     allData = GGIpack::allrefs( con =con, gene = mygene, pfiles =pfiles, genelocs = genelocs)
   })#allrefs
  

  
  output$BALstuff = DT::renderDT({
   refs = allCelltypes()
   refs[["BAL"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })#BALstuff
  
  output$BEBstuff = DT::renderDT({
   refs = allCelltypes()
   refs[["BroncEpiBrush"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })#BEBstuff
  
  output$CD4stim = DT::renderDT({
   refs = allCelltypes()
   refs[["CD4Stim"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   }) #CD4Stim
  
   output$CD4Unstim = DT::renderDT({
   refs = allCelltypes()
   refs[["CD4Unstim"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })#CD4Unstim
   
   output$AlvMacphage = DT::renderDT({
   refs = allCelltypes()
   refs[["AlvMacphage"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })#AlvMacphage
   
   output$PaxRNA = DT::renderDT({
   refs = allCelltypes()
   refs[["PaxRNA"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   }) #paxRNA
   
   output$parquet<- DT::renderDT({
     data.frame(patquetTable)
   })#parquet
   ##################################################################################################
   
   
   
   
   ################################################################################################################
   ##igvshiny
   ############################################################################################################### 
   
   output$igvShiny_0 = igvShiny::renderIgvShiny({
     genomeOptions <- igvShiny::parseAndValidateGenomeSpec(genomeName = genomeVersion, initialLocus = "all")
     igvShiny(genomeOptions)
   })#igvShiny_0
   
  
   #dataToGraph = c("BALstuff", "BEBstuff","CD4stim","CD4Unstim", "AlvMacphage", "PaxRNA")
   
   
   #names(dataToGraph) = c("BAL", "BronchEpiBrush", "CD4stim","CD4Unstim", "AlvMacphage", "PaxRNA")
  
   observeEvent(input$gene, {
     refs = allCelltypes()
     print(names(refs))
   #  print(refs[[1]]@tbl)
     for(i in 1:length(refs)){
        GenomicTable = refs[[i]]@tbl |>  as.data.frame() 
       print(table)
       gwasTrack = makeGWASTrack( name=names(GenomicTable), dat = GenomicTable)
       display(gwasTrack, session, id = "igvShiny_0")
     } #for loop
     #tableDn8like =  allfilt[[1]]@tbl |> as.data.frame() 
     #genomicRegion = paste0("chr", min(tableDn8like$CHR),":", formatC(min(tableDn8like$BP)-bpPadding , format="d", big.mark = ","), "-", formatC(max(tableDn8like$BP)+bpPadding, format="d", big.mark = ","), sep ="" )
     #showGenomicRegion(session, "igvShiny_0", genomicRegion)
   }) #observeEvent
   
  # observeEvent(input$zoomButton,{
     
    # tableDn8like = as.data.frame(dataToGraph[[1]])
    # genomicRegion = paste0("chr", min(tableDn8like$CHR),":", formatC(min(tableDn8like$BP)-bpPadding , format="d", big.mark = ","), "-", formatC(max(tableDn8like$BP)+bpPadding, format="d", big.mark = ","), sep ="" )
    # showGenomicRegion(session, "igvShiny_0", genomicRegion)
   #})
   
   
   
#####################################################################################################################   
#output for main panel
#################################################################################################################
    output$alltabs = renderUI({
   tabsetPanel(
     tabPanel(igvShiny::igvShinyOutput("igvShiny_0"),
              #shinyFeedback::useShinyFeedback()
              ),
    tabPanel("BAL",  DT::DTOutput("BALstuff")),
    tabPanel("BronchEpiBrush", DT::DTOutput("BEBstuff")),
    tabPanel("CD4stim", DT::DTOutput("CD4stim")),
    tabPanel("CD4Unstim", DT::DTOutput("CD4Unstim")),
    tabPanel("AlvMacphage", DT::DTOutput("AlvMacphage")),
    tabPanel("PaxRNA", DT::DTOutput("PaxRNA")),
    tabPanel("about", helpText(h3("GGIpack Overview")),
             br(),
             p(sprintf(
               "GGIpack abrigapp version %s.  This app uses parquet files made 
                 from the abrig  data release on 05/15/2023. This new data release merges the 
                 population data for each of the cell types together thus there is no more choice 
                 for population since there is no way to separate the data anymore. The original
                 data can be found in the following  path on the Nantucket server.",
               packageVersion("GGIpack2") ) 
             ),
             br(),
             p("/proj/regeps/regep00/studies/ABRIG/analyses/reahs/cis_eqtl_matrixEqtl.Release.15.05.23/"),
             br(),
             p(" The parquet files were made using the arrow package an R package.  
               Then the parquet files were processed/querried using the  duckdb R package.
               In short these fileswere made by merging all the 23 chromosome files for each of the cell types into one file.
               The table below shows the name of the cell type, the name that it is called in the app, and the name of the actual file."),
             br(),
             DT::DTOutput("parquet"),
             p("The files can be found in the Nantucket server under."),
             br(),
             p("/udd/remcr/abrig/"),
             br(),
             p("The code that was used to create these files can be found in the following 
                   Changit Repository."),
             br(),
             p("https://changit.bwh.harvard.edu/remcr/abrigResources")
       ) #tabPanel for about page
    ) #tabsetPanel
   }) #renderUI
  ##############################################################
   
 }#end of server
 runApp(list(ui=ui, server=server))
 DBI::dbDisconnect(con)
} #end of abrigapp

