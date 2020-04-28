# TODO:
# Get ngl viewer to read local files?
# Get Download button to work?
# Consider modularising bits.
# Add write-permissions for the app as a whole?
#gpath <- '.'
gpath <- '/srv/shiny-server/'
install.packages(sprintf('%s/%s', gpath, 'nglShiny'), type='source', repos=NULL)
library(devtools)
library(shiny)
library(DT)
library(htmlwidgets)
library(nglShiny)
library(caTools)
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

nglRepresentations = c(
    'angle', 
    'axes', 
    'ball+stick', 
    'backbone', 
    'base', 
    'cartoon', 
    'contact',
    'dihedral', 
    'distance', 
    'helixorient', 
    'licorice', 
    'hyperball', 
    'label',
    'line', 
    'surface', 
    'point', 
    'ribbon', 
    'rocket', 
    'rope', 
    'spacefill', 
    'trace', 
    'unitcell',
    'validation'
    )

#nglColorSchemes <- c('residueIndex', 'chainIndex', 'entityType', 'entityIndex')
nglColorSchemes <-  c(
    'atomindex',
    'bfactor',
    'chainid',
    'chainindex',
    'chainname',
    'densityfit',
    'electrostatic',
    'element',
    'entityindex',
    'entitytype',
    'geoquality',
    'hydrophobicity',
    'modelindex',
    'moleculetype',
    'occupancy',
    'random',
    'residueindex',
    'resname',
    'sstruc',
    'uniform',
    'value',
    'volume'
    )

defaultRepresentation <- "ball+stick"
defaultColorScheme <- "chainIndex"
possRes <- c("",  "Needs refinement", "Unconvincing ligand density (I\'ve checked event maps)", "Ligand restraints issue", "Low resolution/poor data quality", 'Other')

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
                    # selectizeInput('site', 'Which Site?', list(), multiple=TRUE),
                    selectizeInput('Xtal', 'Which Structure to flag?', list(), multiple = FALSE),
                    #actionButton("download", "Download Data", class = "btn-primary"),
                    selectInput("reason", "Reason for Rejection", possRes),
                    textOutput('msg'),
                    actionButton("submit", "Submit", class = "btn-primary"),
                    selectizeInput('protein', 'Select Protein Rows', list(), multiple=TRUE),
                    selectizeInput('columns', 'Select Columns to View?', list(), multiple = TRUE)
                    #checkboxInput("check1", "This button does nothing", FALSE),
                ), width=2 # div
            ), #sidebarPanel
                            
            mainPanel(
                tabsetPanel(
                    tabPanel("Staged Structures", DT::dataTableOutput("table")),
                    tabPanel("Responses", DT::dataTableOutput("resp")),
                    tabPanel("Help", includeMarkdown(sprintf('%s/%s', gpath, "Pages/include.md")))
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
                textInput("name2", "Name", ""),
                selectizeInput('Xtal2', 'Which Structure to View/Flag?', list('Mpro-x0104'), multiple = FALSE),
                selectInput("reason2", "Reason for Rejection", possRes),
                textOutput('stats'),
                hr(),
                textOutput('msg2'),
                actionButton("submit2", "Submit", class = "btn-primary"),
                hr(),
                actionButton("fitButton", "Fit"),
                actionButton("defaultViewButton", "Defaults"),
                actionButton("clearRepresentationsButton", "Clear Representations"),

                #selectizeInput('pdbSelector', 'Which Structure to view?', list("rcsb://1crn"), multiple = FALSE), # Replace with Xtal2?

                selectInput("representationSelector", "", nglRepresentations, selected=defaultRepresentation),
                selectInput("colorSchemeSelector", "", nglColorSchemes, selected=defaultColorScheme),
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
defaultPdbID <- ""
# Things in Global Scope
responsesDir <- file.path(sprintf('%s/%s', gpath, "Responses"))
dataDir <- file.path(sprintf('%s/%s', gpath, "Data"))
pdbIDs <- dir(dataDir, pattern='.pdb', full = TRUE, rec=TRUE)
mapIDs <- dir(dataDir, pattern='.ccp4', full = TRUE, rec=TRUE)
    # Stuff people shouldn't see.
    #options <- list(pdbID="1pcr")
    #options <- list(pdbID="3kvk")
    options <- list(pdbID="")
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

    fieldsAll <- c("name", 'Xtal', "reason")
    formData <- reactive({
        data <- sapply(fieldsAll, function(x) input[[x]])
        data <- c(data, timestamp = epochTime())
        data <- t(data)
        data
    })

    fieldsAll2 <- c("name2", 'Xtal2', "reason2")
    formData2 <- reactive({
        data <- sapply(fieldsAll2, function(x) input[[x]])
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
  
    r2 <- reactive({
        dat <- inputData()[input$Xtal2, ]
        sn <- names(dat)[2:8]
        val <- dat[2:8]
        values <- sprintf('%s: %s \n %s: %s \n %s: %s \n %s: %s \n %s: %s \n %s: %s \n',
            sn[1], val[1], sn[2], val[2], sn[3], val[3], sn[4], val[4], sn[5], val[5], sn[6], val[6] 
            )
        values
    })

    output$stats <- renderText({r2()})

    output$table <- DT::renderDataTable({r1()})

    # Response Table
    output$resp <- DT::renderDataTable(
        loadData(),
        rownames = FALSE,
        options = list(searching = FALSE, lengthChange = FALSE)
    ) 
  
    output$msg <- renderText({'Please click once, \n refresh page to see response being updated'})  
    output$msg2 <- renderText({'Please click once, \n refresh page to see response being updated'})  
    output$value <- renderPrint({ input$action})
  
    # Events
    observeEvent(input$submit, {
        print(formData())
        saveData(formData())
    })
    
    observeEvent(input$submit2, {
        print(formData2())
        saveData(formData2())
    })

    observe({
        updateSelectizeInput(session, 'columns', choices = colnames(inputData()))
        updateSelectizeInput(session, 'protein', choices = sort(unique(inputData()$Protein)))
        updateSelectizeInput(session, "Xtal", choices = sort(rownames( inputData() )))
        #updateSelectizeInput(session, "Xtal2", choices = sort(rownames( inputData() )))
        updateSelectizeInput(session, "Xtal2", choices = c('Mpro-x0104', 'Mpro-x0161') )
        #updateSelectizeInput(session, 'pdbSelector', choices=sort(pdbIDs))
    })

    observeEvent(input$fitButton, {
        session$sendCustomMessage(type="fit", message=list())
    })
  
    observeEvent(input$clearRepresentationsButton, {
        session$sendCustomMessage(type="removeAllRepresentations", message=list())
        #updateSelectInput(session, "representationSelector", label=NULL, choices=NULL,  selected=defaultRepresentation)
        #updateSelectInput(session, "colorSchemeSelector", label=NULL, choices=NULL,  selected=defaultColorScheme)
    })
  
    observeEvent(input$Xtal2, {
        fail <- try({
        choice = input$Xtal2
        if(grepl('rcsb', choice)){
            message(sprintf("pdb: %s", choice))
            session$sendCustomMessage(type="setPDB", message=list(choice))
        } else {
            # From Xtal2, look in folder for .pdb and then ligand centered Map.
            # If pdb is not on pdb... Do things.
            message(sprintf("pdb: %s", choice))
            syscall <- sprintf('cat %s', dir(sprintf('%s/%s',dataDir, choice), pattern = 'pdb', full.names=T))
            message(syscall)
            pdbstrings <- system(syscall, intern = TRUE)
            fname <- dir(sprintf('%s/%s',dataDir, choice), pattern = 'ccp4', full.names=T)
            message(fname)
            choice <- paste0(pdbstrings, collapse='\n')
            defaultPdbID <<- choice
            event <-  readBin(fname, what = 'raw', file.info(fname)$size)
            event <- base64encode(event, size=NA, endian=.Platform$endian)
            defaultShell <<- event
            session$sendCustomMessage(type="setPDB2", message=list(choice))
            session$sendCustomMessage(type="addEvent", message=list(event))
        }
        }, silent=T)
        if(inherits(fail, 'try-error')) session$sendCustomMessage(type="removeAllRepresentations", message=list())
        #updateSelectInput(session, "pdbSelector", label=NULL, choices=NULL,  selected=choice)
    })

    observeEvent(input$defaultViewButton, {
        session$sendCustomMessage(type="removeAllRepresentations", message=list())
        session$sendCustomMessage(type="setPDB2", message=list(defaultPdbID))
        session$sendCustomMessage(type="addEvent", message=list(defaultShell))
        #session$sendCustomMessage(type="setRepresentation", message=list(defaultRepresentation))
        #session$sendCustomMessage(type="setColorScheme", message=list(defaultColorScheme))
        #session$sendCustomMessage(type="fit", message=list())
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
