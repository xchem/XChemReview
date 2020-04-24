# TODO:
# get ngl viewer to read local files?
# Download button?
# Add write-permissions for the app as a whole?

install.packages('/srv/shiny-server/nglShiny', type='source', repos=NULL)
library(nglShiny)
library(shiny)
library(DT)
library(htmlwidgets)
library(devtools)

epochTime <- function() as.integer(Sys.time())
humanTime <- function() format(Sys.time(), "%Y%m%d-%H%M%OS")
nglShiny <- function(options, width = NULL, height = NULL, elementId = NULL)
{
  sprintf("--- ~/github/nglShiny/R/nglShiny ctor");
#  stopifnot("pdbID" %in% names(options))
  
  htmlwidgets::createWidget(
    name = 'nglShiny',
    options,
    width = width,
    height = height,
    # sizingPolicy = htmlwidgets::sizingPolicy(padding=0, browser.fill=TRUE),
    package = 'nglShiny',
    elementId = elementId
  )
} 
defaultPdbID <- "1crn"
nglRepresentations = c('angle', 'axes', 'ball+stick', 'backbone', 'base', 'cartoon', 'contact',
                       'dihedral', 'distance', 'helixorient', 'licorice', 'hyperball', 'label',
                       'line', 'surface', 'point', 'ribbon', 'rocket', 'rope', 'spacefill', 'trace', 'unitcell',
                       'validation')
nglColorSchemes <- c('residueIndex', 'chainIndex', 'entityType', 'entityIndex')
defaultRepresentation <- "cartoon"
defaultColorScheme <- "residueIndex"


# UI Changes
ui <- navbarPage("Staging XChem",           
    # First Page
    tabPanel("Main Page",
        sidebarLayout(
            # Sidebar panel for inputs ----
            sidebarPanel(
                div(
                    id = "form",
                    textInput("name", "Name", ""),
                    selectizeInput('protein', 'Which Protein?', list(), multiple=TRUE),
                    # selectizeInput('site', 'Which Site?', list(), multiple=TRUE),
                    selectizeInput('Xtal', 'Which Structure to flag?', list(), multiple = FALSE),
                    actionButton("download", "Download Data", class = "btn-primary"),
                    selectInput("reason", "Reason for Rejection",
                        c("",  "False Positive", "Needs Refinement", "Just Wrong", 'I don\'t like it')),
                    textOutput('msg'),
                    actionButton("submit", "Submit", class = "btn-primary"),
                    selectizeInput('columns', 'Select Columns to View?', list(), multiple = TRUE)
                    #checkboxInput("check1", "This button does nothing", FALSE),
                ) # div
            ), #sidebarPanel
                            
            mainPanel(
                tabsetPanel(
                    tabPanel("Staged Structures", DT::dataTableOutput("table")),
                    tabPanel("Responses", DT::dataTableOutput("resp")),
                    tabPanel("Help", includeMarkdown("/srv/shiny-server/Pages/include.md"))
                )
            ) # mainPanel
        ) # sidebarLayout
    ), # tabPanel, 
    # End of Page 1.
    # Page 2
    tabPanel('NGL Viewer',
      fluidPage(
      tags$head(
    tags$style("#nglShiny{height:98vh !important;}"),
    tags$link(rel="icon", href="data:;base64,iVBORw0KGgo=")
    ),  

        sidebarLayout(
            sidebarPanel(
                actionButton("fitButton", "Fit"),
                actionButton("defaultViewButton", "Defaults"),
                actionButton("clearRepresentationsButton", "Clear Representations"),
                #shinyFilesButton('pdbSelector', label='File select', title='Which Structure to view?', multiple=FALSE),
                selectizeInput('pdbSelector', 'Which Structure to view?', list(), multiple = FALSE),
                #selectInput("pdbSelector", "", pdbIDs, selected=defaultPdbID),
                selectInput("representationSelector", "", nglRepresentations, selected=defaultRepresentation),
                selectInput("colorSchemeSelector", "", nglColorSchemes, selected=defaultColorScheme),
                hr(),
                width=2
            ), # sidebarPanel
            mainPanel(
                nglShinyOutput('nglShiny'),
                width=10
            ) # mainPanel
        ) # sidebarlayout
    )#, 
    )# tabPanel
    # End of Page 2
) # End of UI

server <- function(input, output, session) {

# Things in Global Scope
fieldsAll <- c("name", 'Xtal', "reason")
responsesDir <- file.path("/srv/shiny-server/Responses")
dataDir <- file.path('/srv/shiny-server/Data/')
pdbIDs <- c("1crn",  # crambin refined against 0.945-A x-ray diffraction data.
            "2UWS",  # photosynthetic reaction center from Rb. sphaeroides, pH 6.5, charge-separated state
            "1IZL",  # Crystal structure of oxygen-evolving photosystem II from Thermosynechococcus vulcanus at 3.7-A resolution
            #dir(dataDir, pattern='.pdb', full=T), # My structure (will need to add a mtz...)
            '6TNU')
            #'/foo/Data/refine_16.pdb',
            #'/foo')



    # Stuff people shouldn't see.
    #options <- list(pdbID="1pcr")
    #options <- list(pdbID="3kvk")
    options <- list(pdbID="1crn")
    #options <- list(pdbID="1rqk")
    output$nglShiny <- renderNglShiny(
        nglShiny(list(), 300, 300)
    )
  
    loadData <- function() {
        files <- list.files(file.path(responsesDir), full.names = TRUE)
        data <- lapply(files, read.csv, stringsAsFactors = FALSE)
        data <- do.call(rbind, data)
        data
    }

    formData <- reactive({
        data <- sapply(fieldsAll, function(x) input[[x]])
        data <- c(data, timestamp = epochTime())
        data <- t(data)
        data
    })
  
    saveData <- function(data) {
        fileName <- sprintf("%s_%s.csv",
                            humanTime(),
                            digest::digest(data))
    
        write.csv(x = data, file = file.path(responsesDir, fileName),
                  row.names = FALSE, quote = TRUE)
    }
  
    # Outputs:
  
    # Main Table handler
    db <- read.csv(paste(dataDir,'mock.csv', sep='/'), stringsAsFactors=F, row.names=1)
    dedupe <- duplicated(db[,1])
    db <- db[!dedupe,]
    rownames(db) <- db[,1]
    db <- db[,-1]
    inputData <- reactive({db})
    r1 <- reactive({
        # Subset data
        if(is.null(input$protein) & is.null(input$columns)) inputData()
        else if(is.null(input$columns) & !is.null(input$protein)) inputData()[inputData()$Protein %in% input$protein, ]
        else if(!is.null(input$columns) & is.null(input$protein)) inputData()[ ,input$columns]
        else inputData()[inputData()$Protein %in% input$protein, input$columns]
    })
  
    output$table <- DT::renderDataTable({r1()})

    # Response Table
    output$resp <- DT::renderDataTable(
        loadData(),
        rownames = FALSE,
        options = list(searching = FALSE, lengthChange = FALSE)
    ) 
  
    output$msg <- renderText({'Please click once, \n refresh page to see response being updated'})  
    output$value <- renderPrint({ input$action})
  
    # Events
    observeEvent(input$submit, {
        print(formData())
        saveData(formData())
    })
  
    observe({
        updateSelectizeInput(session, 'columns', choices = colnames(inputData()))
        updateSelectizeInput(session, 'protein', choices = sort(unique(inputData()$Protein)))
        updateSelectizeInput(session, "Xtal", choices = sort(rownames( inputData() )))
        updateSelectizeInput(session, 'pdbSelector', choices=sort(pdbIDs))
    })
    

    observeEvent(input$fitButton, {
        session$sendCustomMessage(type="fit", message=list())
    })
  
    observeEvent(input$defaultViewButton, {
        session$sendCustomMessage(type="removeAllRepresentations", message=list())
        session$sendCustomMessage(type="setRepresentation", message=list(defaultRepresentation))
        session$sendCustomMessage(type="setColorScheme", message=list(defaultColorScheme))
        session$sendCustomMessage(type="fit", message=list())
    })
  
    observeEvent(input$clearRepresentationsButton, {
        session$sendCustomMessage(type="removeAllRepresentations", message=list())
        #updateSelectInput(session, "representationSelector", label=NULL, choices=NULL,  selected=defaultRepresentation)
        #updateSelectInput(session, "colorSchemeSelector", label=NULL, choices=NULL,  selected=defaultColorScheme)
    })
  
    #shinyFileChoose(input, 'pdbSelector', root=c(root='.'), filetypes=c('pdb'))

    observeEvent(input$pdbSelector, {
        choice = input$pdbSelector
        #choice = paste0(c('.', unlist(choice)[2:3]), collapse='/')
        message(sprintf("pdb: %s", choice))
        session$sendCustomMessage(type="setPDB", message=list(choice))
        #updateSelectInput(session, "pdbSelector", label=NULL, choices=NULL,  selected=choice)
    })
  
    observeEvent(input$representationSelector, {
        choice = input$representationSelector;
        message(sprintf("rep: %s", choice))
        session$sendCustomMessage(type="setRepresentation", message=list(choice))
        updateSelectInput(session, "representationSelector", label=NULL, choices=NULL,  selected=choice)
    })
  
    observeEvent(input$colorSchemeSelector, {
        choice = input$colorSchemeSelector;
        message(sprintf("colorScheme: %s", choice))
        session$sendCustomMessage(type="setColorScheme", message=list(choice))
        updateSelectInput(session, "colorSchemeSelector", label=NULL, choices=NULL,  selected=choice)
    })
  
    # Trigger download for input$Xtal
    output$downloadBtn <- downloadHandler(
        filename = function() { 
            'Disintegrate.png'
        },
        content = function(file) {
            file.copy(file, sprintf('%s',file))
        }
    )
} # Server

# Run the application 
app <- shinyApp(ui = ui, server = server)
runApp(app, host ="0.0.0.0", port = 3838, launch.browser = FALSE)
