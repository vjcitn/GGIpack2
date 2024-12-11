#library(shiny)
#library(JBrowseR)
#library(bslib)
#
#ui <- fluidPage(
#  # Overriding the default bootstrap theme is needed to get proper font size
#  theme = bs_theme(version = 5),
#  titlePanel("JBrowseR Example"),
#  JBrowseROutput("widgetOutput")
#)
#
#server <- function(input, output, session) {
#  # create the assembly configuration
#  assembly <- assembly(
#    "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz",
#    bgzip = TRUE,
#    aliases = c("GRCh37"),
#    refname_aliases = "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/hg19_aliases.txt"
#  )
#  
# #df <- data.frame(
#  #chrom = c('1', '2'),
#   # start = c(123, 456),
#    #end = c(789, 101112),
#    #name = c('feature1', 'feature2')
#  #)
#  
#  path = system.file("extdata/Alveolar_Macrophages_IS.MICA:ILMN_3241692.CAU.meta", package="GGIpack")
#  df_origin = utils::read.table(path, header = TRUE)
#  df <- data.frame(
#    chrom =df_origin$CHR ,
#    start =df_origin$BP,
#    end = df_origin$BP,
#    name = df_origin$SNP
#  )
#  
#  
#  df_track <- JBrowseR::track_data_frame(df, "foo", assembly)
#  
#  # set up the final tracks object to be used
#  tracks <- tracks(
#    df_track
#  )
#  
#  # determine what the browser displays by default
#  default_session <- default_session(
#    assembly,
#    c(df_track),
#    display_assembly = FALSE
#  )
#  
#  
#  output$widgetOutput <- JBrowseR::renderJBrowseR(
#    JBrowseR::JBrowseR("View",
#             assembly = assembly,
#             tracks = tracks,
#             location = "2:456",
#             defaultSession = default_session
#    )
#  )
#}
#
#shinyApp(ui, server)
