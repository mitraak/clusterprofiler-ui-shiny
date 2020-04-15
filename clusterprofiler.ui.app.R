library(ggplot2)
library(AnnotationHub)
library(clusterProfiler)
library(DT)
library(waiter)

plot.select.pathways <- function(df, 
                           pathway.list,
                           alpha=0.05){

    if(is.null(df)) return(NULL)
    idx <- tolower(df$Description) %in% tolower(pathway.list)
    df.subset <- df[idx,]

    # plot horizontal bar chart of selected pathways
    p <- ggplot(df.subset, aes(x=reorder(Description, p.adjust), y=p.adjust, fill=Count)) +
            geom_bar(stat='identity') + coord_flip() + 
            scale_y_continuous(trans='log10') +
            theme_bw() + xlab('Description') + geom_hline(yintercept=alpha, linetype=2)
    return(p)
}

# from lcdb-wf
get.orgdb <- function(species, cache, annotation_key_override=NA){
  
  # Workaround to allow AnnotationHub to use proxy. See
  # https://github.com/Bioconductor/AnnotationHub/issues/4, and thanks
  # Wolfgang!
  proxy <- Sys.getenv('http_proxy')
  if (proxy == ""){
    proxy <- NULL
  }
  
  ah <- AnnotationHub(hub=getAnnotationHubOption('URL'),
                      cache='AnnotationHubCache',
                      proxy=setAnnotationHubOption('PROXY', proxy),
                      localHub=FALSE)
  
  find.annotationhub.name <- function(species.name, override.code) { #autodetect ah names based on loaded database
    if (is.na(override.code)) {
      ah.query <- query(ah, "OrgDb")
      ah.query.speciesmatch <- grepl(paste("^", species.name, "$", sep=""), ah.query$species)
      ah.query.which <- which(ah.query.speciesmatch)
      stopifnot(length(ah.query.which) > 0) #require at least one match
      if (length(ah.query.which) > 1) { #warn of duplicate matches
        print("WARNING: found multiple candidate species in AnnotationHub: ");
        print(ah.query.speciesmatch)
      }
      names(ah.query)[ah.query.which[1]]
    } else {
      override.code
    }
  }
  annotation_key <- find.annotationhub.name(species, annotation_key_override)
  orgdb <- ah[[annotation_key]]
  return(orgdb)
}

# function to perform enrichment
enrich <- function(gene, ont, 
                   orgdb,
                   keytype,
                   species=NULL,
                   universe=NULL,
                   pvalueCutoff=1,
                   qvalueCutoff=1){
  ont <- toupper(ont)
  if(ont %in% c('BP','CC','MF')){
    res <- enrichGO(gene,
                    universe=universe,
                    keyType=keytype,
                    OrgDb=orgdb,
                    ont=ont,
                    pvalueCutoff=pvalueCutoff,
                    qvalueCutoff=qvalueCutoff,
                    readable=TRUE)
  } else if(ont == 'KEGG'){
    org.split <- unlist(strsplit(species, '\\s+'))
    organism <- tolower(paste0(substr(org.split[1],1,1),
                               substr(org.split[2],1,2)))
    
    # get gene symbols for kegg result
    df <- bitr(gene, fromType=keytype, toType=c('SYMBOL','ENTREZID'), OrgDb=orgdb)
    
    res <- enrichKEGG(df$ENTREZID,
                      universe=universe,
                      organism=organism,
                      keyType='ncbi-geneid',
                      pvalueCutoff=pvalueCutoff,
                      qvalueCutoff=qvalueCutoff
                      )
    
    # convert ENTREZID to SYMBOL
    id.vec <- df$SYMBOL
    names(id.vec) <- df$ENTREZID
    
    res.genes <- res@result$geneID

    for(j in 1:length(res.genes)){
      id.split <- strsplit(res.genes[j], "/")[[1]]
      temp <- paste(id.vec[id.split], collapse="/")
      res@result$geneID[j] <- temp
    }
  }
  res
}

ui <- fluidPage(
    # Title of the page
    titlePanel(title='Functional enrichment analysis with clusterProfiler'), 
    
    # main body of the page
    mainPanel(width=12,
      fluidRow(
        column(3,
          # Choose TSV or file with enrichment object
          radioButtons('filetype', label=h4('Type of input'),
                       choices=c('gene list','gene file', 'enrichment object'), 
                       selected='gene list')),
        
          
        column(4,
          # conditional controls
          
          # paste list of genes
          conditionalPanel(
            condition = "input.filetype == 'gene list'",
            textInput('gene.list', label=h4('Genes to analyze'))
          ),
          
          # upload list of genes
          conditionalPanel(
            condition = "input.filetype == 'gene file'",
            fileInput('gene.file', label=h4('Upload file'))
          ),
          
          # upload RDS file
          conditionalPanel(
            condition = "input.filetype == 'enrichment object'",
            fileInput('enrich.file', label=h4('Upload file'))
          )
        ), 
        
        column(4,
               
          # show analysis controls for pasted genes or gene file
          conditionalPanel(
             condition = "input.filetype == 'gene list' || input.filetype == 'gene file'",
             
             h4('Analysis options'),
             # select pathways to analyze
             checkboxGroupInput('pathway', label=h5('Pathways'),
                                c('GO Biological Process'='BP', 
                                  'GO Cellular Component'='CC', 
                                  'GO Molecular Function'='MF', 
                                  'KEGG Pathways'='KEGG')),

             # select species
             selectInput('species', label=h5('Species'),
                         choices=c('Homo sapiens',
                                   'Mus musculus',
                                   'Drosophila melanogaster'),
                         selected='Homo sapiens'),
             # select gene ID source
             selectInput('keytype', label=h5('Gene ID source'),
                         choices=c('ENSEMBL','FLYBASE'),
                         selected='ENSEMBL'),
            
             # submit
             actionButton('go.pathway', label='Submit')
          )
        )
      ),
      hr(),

    # navigation panel with section tabs
    navlistPanel(widths=c(3,9),
      'Visualization',
      id='viz.options',

      # tabPanel to choose enrichment to visualize
      tabPanel(
        'Choose enrichment',
        
        # rds file or after enrichment
        fluidRow(
          use_waiter(),
          column(10, uiOutput('object.contents'))
        ),
        downloadButton('rds_download', label = 'Download results')
      ),
      
      # tabPanel with data table
      tabPanel(
        # title of tabPanel
        'Data table',

        # tabPanel elements
        h2('Enrichment analysis table'),        
        fluidRow(
          column(12, dataTableOutput('contents'))
        )),
        
      # tabPanel with barplot
      tabPanel(
        'Bar plot',

        h3('Selected pathways from table'), br(),
        tags$i('By default we plot all visible rows'), br(),
        plotOutput('barplot', width='100%'),
        downloadButton('barplot_download', label = 'Download plot')
      ),
      
      # tabPanel with enrichment plots
      tabPanel(  
        'Plots',

        # choose plot type
        selectInput('plottype', label=h4('Type of plot'),
                    choices=c('dotplot','emapplot','cnetplot'), selected = 'dotplot'),
        
        # conditional controls

        # dotplot
        conditionalPanel(
          condition = "input.plottype == 'dotplot'",
          fluidRow(
            column(3, selectInput('xvar.dot', label=h6('X-axis variable'),
                                  choices=c('GeneRatio', 'Count'),  selected='GeneRatio')),
            column(3, selectInput('color.dot', label=h6('Color'),
                        choices=c('p.adjust', 'pvalue', 'qvalue'),  selected='p.adjust')),
            column(3, numericInput('numcat.dot', label=h6('# of terms'), value=10)))
        ),

        # emapplot
        conditionalPanel(
         condition = "input.plottype == 'emapplot'",
         fluidRow(
           column(4, numericInput('numcat.emap', label=h6('# of terms'), value=30)),
           column(4, selectInput('color.emap', label=h6('Color'),
                                 choices=c('p.adjust', 'pvalue', 'qvalue'),  selected='p.adjust')))
        ),
        
        # cnetplot
        conditionalPanel(
          condition = "input.plottype == 'cnetplot'",
          fluidRow(
            column(2, numericInput('numcat.cnet', label=h6('# of terms'), value=5)),
            column(3, selectInput('colorEdge.cnet', label=h6('Color edge by terms?'),
                                  choices=c(FALSE, TRUE),  selected=FALSE)),
            column(2, selectInput('circular.cnet', label=h6('Circular layout?'),
                                  choices=c(FALSE, TRUE),  selected=FALSE)),
            column(3, selectInput('nodelabel.cnet', label=h6('Show node label?'),
                                  choices=c(FALSE, TRUE),  selected=TRUE)))
        ),
        
        # output plot
        plotOutput('enrichplot', height='600px', width='100%'),
        downloadButton('enrichplot_download', label = 'Download plot')
        )
    )
    )
    
)

# function to insert commas instead of '/' between gene names
sanitize.gene.names <- function(df, col='geneID'){
    unlist(lapply(strsplit(df[, col], '/'),
                  function(x){ paste(x, collapse=',')}
    ))
}

server <- function(input, output, session){
    # increase upload size
    options(shiny.maxRequestSize=30*1024^2)

    # reactive values to keep track of enrichment object
    enrich.list <- reactiveValues(l = list())
    
    # reactive element for genes
    genes <- reactiveValues(g=NULL)
    
    # reactive values for gene file
    gene.file <- reactiveValues(df=NULL)
    
    # update gene.file if new file is input
    observeEvent(input$gene.file, {
      gene.file$df <- read.table(input$gene.file$datapath, 
                                 sep='\t', header=FALSE, 
                                 stringsAsFactors=FALSE)
    })
    
    # observe event to get gene list
    observeEvent(input$gene.list, { 
      if(!is.null(input$gene.list)){
        genes$g <- input$gene.list
      }
    })
    
    # observe event to get gene file
    observeEvent(gene.file$df, {
      genes$g <- gene.file$df[,1]
    })
    
    # flush enrichment list if new gene file or list is uploaded
    do.flush <- reactive({
      c(input$gene.list, input$gene.file)
    })
    
    observeEvent(do.flush(), {
      enrich.list$l <- list()
    })
    
    # observe event to perform enrichment
    observeEvent(input$go.pathway, {
      show_waiter(spin_chasing_dots())
      
      # find list of pathways that have not been analyzed
      pathway <- input$pathway[!input$pathway %in% names(enrich.list$l)]
      l <- get.enrichment(genes$g, pathway,
                          input$species, input$keytype)
        for(name in names(l)){
          enrich.list$l[[name]] <- l[[name]]
        }
      hide_waiter()
      })
    
    # function to do enrichment
    get.enrichment <- function(gene, pathway, 
                               species, keytype){
      if(is.null(gene)) return(NULL)
      g <- toupper(unlist(strsplit(gene, '\\,*\\s+')))
      
      p <- lapply(pathway, toupper)
      
      orgdb <- get.orgdb(species=species)
      l <- list()
      for(ont in p){
        cat(ont,'\n')
        l[[ont]] <- enrich(gene=g, ont=ont,
                           orgdb=orgdb, keytype=keytype,
                           species=species)
      }
      cat('Enrichment done!\n')
      l
    }
    
    # download enrichment results to RDS file
    output$rds_download <- downloadHandler(
      filename = function(){
        paste0('results.rds')
      },
      content = function(file){
        #ggsave(file, barplot(), width=20, height=10)
        saveRDS(enrich.list$l, file)
      }
    )
    # observe event to read RDS file
    observeEvent(input$enrich.file, {
      l <- readRDS(input$enrich.file$datapath)
      enrich.list$l <- l
    })
    
    # barplot of selected pathways
    barplot <- reactive({
        if(!is.null(enrich.list$l)){
          l <- enrich.list$l
        }
      
        df <- as.data.frame(l[[input$elements]])
        
        if(is.null(input$contents_rows_selected)){
          idx <- input$contents_rows_current
        } else {
          idx <- input$contents_rows_selected
        }
        current.pathways <- df[idx, 'Description']
        p <- plot.select.pathways(df,
                                  current.pathways)
        p + theme(text=element_text(size=20))
    })
    
    # output slot with barplot
    output$barplot <- renderPlot({
      if(is.null(input$enrich.file) & is.null(enrich.list$l)) return(NULL)
      barplot()
    })
    
    # download button for barplot
    output$barplot_download <- downloadHandler(
      filename = function(){
        paste0('barplot-selected.pdf')
      },
      content = function(file){
        ggsave(file, barplot(), width=20, height=10)
      }
    )

    # data table of enrichment results
    output$contents <- renderDataTable({
        if(is.null(input$enrich.file) & is.null(input$gene.list)) return(NULL)
      
        if(!is.null(enrich.list$l)){
          l <- enrich.list$l
        }
      
        df <- as.data.frame(l[[input$elements]]) 
        df$geneID <- sanitize.gene.names(df)
        df %>% dplyr::select(-'ID') %>% 
          dplyr::select(1, 'Count', dplyr::everything()) %>%
          datatable %>% 
          formatSignif(columns=c('p.adjust','pvalue','qvalue'), digits=4) %>%
          formatStyle('Description', 'white-space'='nowrap')
    })
    
    # reactive UI based on uploaded object
    output$object.contents <- renderUI({
      if(!is.null(enrich.list$l)){
        elements <- names(enrich.list$l)
      }
      selectInput('elements', label=h4('Choose enrichment'), choices=elements)
    })
    
    # reactive function for enrichplot
    enrichplot <- reactive({
      if(!is.null(enrich.list$l)){
        l <- enrich.list$l
      }
      
      if(input$plottype == 'dotplot'){
        p <- dotplot(l[[input$elements]],
                     x = input$xvar.dot,
                     color = input$color.dot,
                     showCategory = input$numcat.dot)
        p <- p + theme(axis.title.x=element_text(size=20),
                       axis.title.y=element_text(size=20),
                       axis.text.x=element_text(size=15),
                       axis.text.y=element_text(size=17))
      } else if(input$plottype == 'emapplot'){
        p <- emapplot(l[[input$elements]],
                      showCategory = input$numcat.emap,
                      color = input$color.emap)
      } else if(input$plottype == 'cnetplot'){
        p <- cnetplot(l[[input$elements]],
                      showCategory = input$numcat.cnet,
                      foldChange = input$fc.cnet,
                      colorEdge = as.logical(input$colorEdge.cnet),
                      circular = as.logical(input$circular.cnet),
                      node_label = as.logical(input$nodelabel.cnet))
      }
      p + theme(text=element_text(size=20))
    })
    
    # output slot with enrichplot
    output$enrichplot <- renderPlot({
      if(is.null(input$enrich.file) & is.null(enrich.list$l)) return(NULL)
      enrichplot()
    })

    # download button for enrichplot
    output$enrichplot_download <- downloadHandler(
      filename = function(){
        paste0(input$plottype, '.pdf')
      },
      content = function(file){
        ggsave(file, enrichplot(), width=20, height=10)
      }
    )
}

shinyApp(ui = ui, server = server)
