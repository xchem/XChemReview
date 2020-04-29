#################################################################################
# Structure is horrible, Time for some organisation:
# Libraries and function definitions
#################################################################################
debug = TRUE
# Set Path: May need to add something later for files on /dls
gpath <- '.'
#gpath <- '/srv/shiny-server/'

# Load Required packages:
# Installing home-brewed version of nglShiny Package as we some source changes.
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
  htmlwidgets::createWidget(
    name = 'nglShiny',
    options,
    width = width,
    height = height,
    package = 'nglShiny',
    elementId = elementId
  )
} 

#################################################################################
# Global Variables used in server/UI
#################################################################################
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

possDec <- c("", "Release", "Release (notify)", "More Work", "Reject")
possAns <- possAns2 <- c('Select Decision')

possDec <- c("", "Release", "Release (notify)", "More Work", "Reject")
possRes <- list('Release' = c('Everything is Wonderful'),
                'Release (notify)' = c(
                    'Alternate binding conformation',
                    'Incomplete Density',
                    'Weak Density',
                    'Low Resolution',
                    'Poor Data quality'
                    ),
                'More Work' = c(
                    'Repeat Experiment',
                    'Check Geometry',
                    'Check Conformation',
                    'Check Refinement'
                    ),
                'Reject' = c(
                    'Density too weak',
                    'Insubstantial Evidence',
                    'Bad coordination',
                    'Incomplete Density'
                    )
                )

possAns <- possAns2 <- c('Select Decision')
defaultPdbID <- ""

#################################################################################
# UI Code
#################################################################################

ui <- navbarPage("Staging XChem", id='beep',          
    # First Page
    tabPanel("Main Page",
        sidebarLayout(
            # Sidebar panel for inputs ----
            sidebarPanel(
                div(
                    id = "form",
                    textInput("name", "Name", ""),
                    # selectizeInput('site', 'Which Site?', list(), multiple=TRUE),
                    selectizeInput('Xtal', 'Which Structure?', list(), multiple = FALSE),
                    #actionButton("download", "Download Data", class = "btn-primary"),
                    selectInput("decision", "Decision", possDec),
                    selectizeInput("reason", "Reason(s)", list(), multiple=TRUE),
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
                    selectizeInput('Xtal2', 'Which Structure?', list('Mpro-x0104'), multiple = FALSE),
                    selectInput("decision2", "Decision", possDec),
                    selectizeInput("reason2", "Reason(s)", list(), multiple=TRUE),
                    textOutput('stats'),
                    hr(),
                    textOutput('msg2'),
                    actionButton("submit2", "Submit", class = "btn-primary"),
                    actionButton("Back", "Back", class = "btn-primary"),
                    hr(),
                    actionButton("fitButton", "Fit"),
                    actionButton("defaultViewButton", "Defaults"),
                    actionButton("clearRepresentationsButton", "Clear Representations"),
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


#################################################################################
# Server Code
#################################################################################
server <- function(input, output, session) {
    sessionTime <- epochTime()
    # Things in Global Scope
    responsesDir <- file.path(sprintf('%s/%s', gpath, "Responses"))
    dataDir <- file.path(sprintf('%s/%s', gpath, "Data"))
    pdbIDs <- dir(dataDir, pattern='.pdb', full = TRUE, rec=TRUE)
    mapIDs <- dir(dataDir, pattern='.ccp4', full = TRUE, rec=TRUE)
    options <- list(pdbID="")

    # Functions
    # Read-in Responses Form Data
    loadData <- function() {
        files <- list.files(file.path(responsesDir), full.names = TRUE)
        data <- lapply(files, read.csv, stringsAsFactors = FALSE)
        data <- do.call(rbind, data)
        data
    }
    # Save Responses.
    saveData <- function(data) {
        fileName <- sprintf("%s_%s.csv",
                            humanTime(),
                            digest::digest(data))
        write.csv(x = data, file = file.path(responsesDir, fileName),
                  row.names = FALSE, quote = TRUE)
    }

    # Reactives
    # Page 1 Form Handler
    fieldsAll <- c("name", 'Xtal', "decision", "reason")
    formData <- reactive({
        data <- sapply(fieldsAll, function(x) paste0(input[[x]], collapse='; '))
        data <- c(data, timestamp = epochTime())
        data <- t(data)
        data
    })

    # Page 2 Form Handler
    fieldsAll2 <- c("name2", 'Xtal2', "decision2", "reason2")
    formData2 <- reactive({
        data <- sapply(fieldsAll2, function(x) paste0(input[[x]], collapse='; '))
        names(data) <- gsub('2','',names(data))
        data <- c(data, timestamp = epochTime())
        data <- t(data)
        data
    })


    # Outputs 

    # Main Table Output Handler
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
  
    output$table <- DT::renderDataTable({r1()}, selection = 'single')

    # NGL viewer side panel stats
    r2 <- reactive({
        dat <- inputData()[input$Xtal2, ]
        sn <- names(dat)[2:8]
        val <- dat[2:8]
        values <- sprintf('%s: %s \n %s: %s \n %s: %s \n %s: %s \n %s: %s \n %s: %s \n %s: %s \n',
            sn[1], val[1], sn[2], val[2], sn[3], val[3], sn[4], val[4], sn[5], val[5], sn[6], val[6], sn[7], val[7] 
            )
        values
    })

    output$stats <- renderText({r2()})

    # NGL Viewer
    output$nglShiny <- renderNglShiny(
        nglShiny(list(), 300, 300)
    )

    # Response Table
    output$resp <- DT::renderDataTable(
        loadData(),
        rownames = FALSE,
        options = list(searching = FALSE, lengthChange = FALSE)
    ) 
  
    # Generic Output Messages.
    output$msg <- renderText({'Please click once'})  
    output$msg2 <- renderText({'Please click once'})  


    # Observers, behaviour will be described as best as possible

    # Upon Row Click
    observeEvent(input$table_rows_selected, {
        # Check if Row has been updated since session began, ensure that loadData()[,] # will also get relevant xtal data?
        rdat <- r1()[input$table_rows_selected,]
        if(sessionTime > max( loadData()[,'timestamp']) ){ 
            # Update Form window (weird bug with changing decision reupdates form...)
            updateSelectizeInput(session, "Xtal", selected = rownames(rdat), choices = sort(rownames( inputData() )))
            updateSelectizeInput(session, "Xtal2", selected = rownames(rdat), choices = sort(rownames( inputData() )))
            # Move to NGL viewer Page
            updateTabsetPanel(session, "beep", selected = 'NGL Viewer')

        } else {
            # Show Dialog that things have changed, allow user to restart session (OK) or cancel out and look at something else
            showModal(modalDialog(title = "Someone has already reviewed this crystal", 
                "Someone has recently reviewed this structure. Restarting the session to capture their reponse, you can then review the structure or choose another."
                , easyClose=TRUE, footer = tagList( modalButton("Cancel"), actionButton("ok", "Restart Session"))
            ))
        }
    })
    observeEvent(input$ok, {
        session$reload()
    })
  
    resetForm <- function(){
        updateSelectizeInput(session, "Xtal", selected = '', choices = sort(rownames( inputData() )))
        updateSelectizeInput(session, "Xtal2", selected = '', choices = sort(rownames( inputData() )))
        updateTabsetPanel(session, "beep", selected = 'Main Page') 
        session$reload()
    }

    # Upon Main Page Submit
    observeEvent(input$submit, {
        # Add check for recent updates?
        if(debug) print(formData())
        saveData(formData())
        resetForm()
    })
    
    # Upon NGL Viewer Page Submit
    observeEvent(input$submit2, {
        # Add check for recent updates
        if(debug) print(formData2())
        saveData(formData2())
        resetForm()
    })

    # Change reasons based on decisions
    # Main Page
    observeEvent(input$decision,{
        possAns <- possRes[[input$decision]]
    })

    # NGL page
    observeEvent(input$decision2,{
        possAns2 <- possRes[[input$decision2]]
    })

    # Upon pressing Fit, Fit structure in window
    observeEvent(input$fitButton, {
        session$sendCustomMessage(type="fit", message=list())
    })

    # Upon pressing Clear, Remove Everything
    observeEvent(input$clearRepresentationsButton, {
        session$sendCustomMessage(type="removeAllRepresentations", message=list())
        #updateSelectInput(session, "representationSelector", label=NULL, choices=NULL,  selected=defaultRepresentation)
        #updateSelectInput(session, "colorSchemeSelector", label=NULL, choices=NULL,  selected=defaultColorScheme)
    })
  
    # When input$Xtal2 (structure dropdown on NGL viewer page) is updated do things
    # Load structure and event to NGL stage!
    observeEvent(input$Xtal2, {
        fail <- try({
            choice = input$Xtal2
            if(grepl('rcsb', choice)){ # If a pdb structure do regular stuff
                if(debug) message(sprintf("pdb: %s", choice))
                session$sendCustomMessage(type="setPDB", message=list(choice))
            } else {
                # If pdb is not on pdb... Do things.
                if(debug) message(sprintf("pdb: %s", choice))
                syscall <- sprintf('cat %s', dir(sprintf('%s/%s',dataDir, choice), pattern = 'pdb', full.names=T))
                if(debug) message(syscall)
                pdbstrings <- system(syscall, intern = TRUE)
                fname <- dir(sprintf('%s/%s',dataDir, choice), pattern = 'ccp4', full.names=T)
                if(debug) message(fname)
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
    })

    # Go back to main Panel, do a refresh for good measure.
    observeEvent(input$Back, {
        resetForm()
    })

    # When pressed re-create original xtal ngl view...
    observeEvent(input$defaultViewButton, {
        session$sendCustomMessage(type="removeAllRepresentations", message=list())
        session$sendCustomMessage(type="setPDB2", message=list(defaultPdbID))
        session$sendCustomMessage(type="addEvent", message=list(defaultShell))
        #session$sendCustomMessage(type="setRepresentation", message=list(defaultRepresentation))
        #session$sendCustomMessage(type="setColorScheme", message=list(defaultColorScheme))
        #session$sendCustomMessage(type="fit", message=list())
    })
    
    # Add defaults
    observeEvent(input$representationSelector, {
        choice = input$representationSelector;
        message(sprintf("rep: %s", choice))
        session$sendCustomMessage(type="setRepresentation", message=list(choice))
        updateSelectInput(session, "representationSelector", label=NULL, choices=NULL,  selected=choice)
    })
  
    # Add colours
    observeEvent(input$colorSchemeSelector, {
        choice = input$colorSchemeSelector;
        message(sprintf("colorScheme: %s", choice))
        session$sendCustomMessage(type="setColorScheme", message=list(choice))
        updateSelectInput(session, "colorSchemeSelector", label=NULL, choices=NULL,  selected=choice)
    })
  
    # Placeholder download button!
    output$downloadBtn <- downloadHandler(
        filename = function() { 
            'Disintegrate.png'
        },
        content = function(file) {
            file.copy(file, sprintf('%s',file))
        }
    )

        # Generic Observers?   
    observe({
        updateSelectizeInput(session, 'columns', choices = colnames(inputData()))
        updateSelectizeInput(session, 'protein', choices = sort(unique(inputData()$Protein)))
        updateSelectizeInput(session, "Xtal", choices = sort(rownames( inputData() )))
        updateSelectizeInput(session, "Xtal2", choices = c('Mpro-x0104', 'Mpro-x0161') )
        updateSelectizeInput(session, 'reason', choices = possRes[[input$decision]])
        updateSelectizeInput(session, 'reason2', choices = possRes[[input$decision2]])
    })
} # Server

#################################################################################
# Runtime
#################################################################################
app <- shinyApp(ui = ui, server = server)
runApp(app, host ="0.0.0.0", port = 3838, launch.browser = FALSE)
