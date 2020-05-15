#################################################################################
# Structure is horrible, Time for some organisation:
# Libraries and function definitions
#################################################################################
{
rm(list=ls())
message <- function (..., domain = NULL, appendLF = TRUE) {
    args <- list(...)
    cond <- if (length(args) == 1L && inherits(args[[1L]], "condition")) {
        if (nargs() > 1L) 
            warning("additional arguments ignored in message()")
        args[[1L]]
    }
    else {
        msg <- .makeMessage(..., domain = domain, appendLF = appendLF)
        call <- sys.call()
        simpleMessage(msg, call)
    }
    defaultHandler <- function(c) {
        cat(conditionMessage(c), file = stdout(), sep = "")
    }
    withRestarts({
        signalCondition(cond)
        defaultHandler(cond)
    }, muffleMessage = function() NULL)
    invisible()
}
debug = TRUE
local = FALSE
# Set Path: May need to add something later for files on /dls

# Server Bindings
gpath <- '/srv/shiny-server/'
responsesDir <- '/dls/science/users/mly94721/xchemreview/Responses/' 
dataDir <- '/dls/science/users/mly94721/xchemreview/Data/'
library(devtools)

if(local){
 gpath <- '.'
 responsesDir <-file.path(sprintf('%s/%s', gpath, "Responses"))
 dataDir <- file.path(sprintf('%s/%s', gpath, "Data"))
 devtools::install_github('tjgorrie/nglShiny')
}
# Load Required packages:
# Installing home-brewed version of nglShiny Package as we some source changes.
#install.packages(sprintf('%s/%s', gpath, 'nglShiny'), type='source', repos=NULL)


library(shiny)
library(DT)
library(htmlwidgets)
library(nglShiny)
library(caTools)
library(DBI)

# Who doesn't love unclosed loops.
#source('./db_config.R')
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

getRootFP <- function(pdbpath){
    splits <- strsplit(pdbpath, split ='/')[[1]]
    n <- length(splits)
    paste0(c(splits[1:(n-2)],''), collapse='/')
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

defaultPdbID <- ""
defaultshell <- ""

possDec_int <- 1:4
names(possDec_int) <- c("Release", "Release (notify)", "More Work", "Reject")

}
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
                    textInput("name", "FedID", ""),
                    # selectizeInput('site', 'Which Site?', list(), multiple=TRUE),
                    selectizeInput('Xtal', 'Which Structure?', list(), multiple = FALSE),
                    #actionButton("download", "Download Data", class = "btn-primary"),
                    selectInput("decision", "Decision", possDec),
                    selectizeInput("reason", "Reason(s)", list(), multiple=TRUE),
                    textOutput('msg'),
                    actionButton("submit", "Submit", class = "btn-primary"),
                    selectizeInput('protein', 'Select Specific Protein', list(), multiple=TRUE),
                    selectizeInput('columns', 'Select Columns to View? (delete/add values as needed)', list(), multiple = TRUE)
                    #checkboxInput("check1", "This button does nothing", FALSE),
                ), width=2 # div
            ), #sidebarPanel
                            
            mainPanel(
                tabsetPanel(
                    tabPanel("Staged Structures", DT::dataTableOutput("table")),
                    #tabPanel("Responses", DT::dataTableOutput("resp")),
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
                    textInput("name2", "FedID", ""),
                    selectizeInput('Xtal2', 'Which Structure?', list(), multiple = FALSE),
                    selectInput("decision2", "Decision", possDec),
                    selectizeInput("reason2", "Reason(s)", list(), multiple=TRUE),
                    textOutput('stats'),
                    hr(),
                    textOutput('msg2'),
                    actionButton("submit2", "Submit", class = "btn-primary"),
                    actionButton("Back", "Back", class = "btn-primary"),
                    hr(),
                    textOutput('msg3'),
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
    if(debug) message('Server Init')
    sessionTime <- epochTime()
    # Things in Global Scope
    #pdbIDs <- dir(dataDir, pattern='.pdb', full = TRUE, rec=TRUE)
    #mapIDs <- dir(dataDir, pattern='.ccp4', full = TRUE, rec=TRUE)
    options <- list(pdbID="")

    # Functions
    # Read-in Responses Form Data
    loadData <- function() {
        files <- list.files(file.path(responsesDir), full.names = TRUE)
        data <- lapply(files, read.csv, stringsAsFactors = FALSE)
        data <- do.call(rbind, data)
        # Now only contains ID instead of xtal name...
        data
    }
    # Save Responses.
    saveData <- function(data) {
        fileName <- sprintf("%s_%s.csv",
                            humanTime(),
                            digest::digest(data))
        if(!data[,'fedid'] %in% c('', ' ')){ # Create Modal that prevent empty data from being submitted!!
            write.csv(x = data, file = file.path(responsesDir, fileName),
                  row.names = FALSE, quote = TRUE)
            con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
            dbAppendTable(con, 'review_responses', value = data, row.names=NULL)
            dbDisconnect(con)
        }
    }

    # Reactives
    # Page 1 Form Handler
    fieldsAll <- c("name", 'Xtal', "decision", "reason")
    formData <- reactive({
        data <- sapply(fieldsAll, function(x) paste0(input[[x]], collapse='; '))
        # Get Crystal ID
        data <- c(dbdat[data[2], 'Id'], data[1], possDec_int[data[3]] ,data[3:4], timestamp = epochTime())
        data <- data.frame(t(data), stringsAsFactors=F)
        # Force Coercion
        data[,1] <- as.integer(data[,1])
        data[,3] <- as.integer(data[,3])
        data[,6] <- as.integer(data[,6])
        colnames(data) <- c('crystal_id', 'fedid', 'decision_int', 'decision_str', 'reason', 'time_submitted')
        data
    })

    # Page 2 Form Handler
    fieldsAll2 <- c("name2", 'Xtal2', "decision2", "reason2")
    formData2 <- reactive({
        data <- sapply(fieldsAll2, function(x) paste0(input[[x]], collapse='; '))
        names(data) <- gsub('2','',names(data))
        # Get Crystal ID
        data <- c(dbdat[data[2], 'Id'], data[1], possDec_int[data[3]] ,data[3:4], timestamp = epochTime())
        data <- data.frame(t(data), stringsAsFactors=F)
        # Force Coercion
        data[,1] <- as.integer(data[,1])
        data[,3] <- as.integer(data[,3])
        data[,6] <- as.integer(data[,6])
        colnames(data) <- c('crystal_id', 'fedid', 'decision_int', 'decision_str', 'reason', 'time_submitted')
        data
    })


    # Outputs, invest in putting in postgres instead of slurping everything into memory
    # Otherwise extremely slow to shove this in the front end? and on session loading...
    # Can use it for
    source('/dls/science/users/mly94721/xchemreview/db_config.R') # Config file...
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    refinement_data <- dbGetQuery(con, "SELECT id, crystal_name_id, r_free, ramachandran_outliers, res, rmsd_angles, rmsd_bonds, lig_confidence_string, cif, pdb_latest, mtz_latest FROM refinement WHERE outcome=4")
    crystal_data <- dbGetQuery(con, sprintf("SELECT id, crystal_name, compound_id, target_id FROM crystal WHERE id IN (%s)", paste(refinement_data[,'crystal_name_id'], collapse=',')))
    target_data <- dbGetQuery(con, sprintf("SELECT * FROM target WHERE id IN (%s)", paste(crystal_data[,'target_id'], collapse=',')))
    compound_data <- dbGetQuery(con, sprintf("SELECT * FROM compounds WHERE id IN (%s)", paste(crystal_data[,'compound_id'], collapse=',')))
    # Sort Responses...
    response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
    dbDisconnect(con)

    comps <- compound_data[,2]
    names(comps) <- as.character(compound_data[,1])
    targs <- target_data[,2]
    names(targs) <- as.character(target_data[,1])

    # Collapse into DF
    jd <- cbind(refinement_data[match(crystal_data[,1], refinement_data[,2]), ], crystal_data[,-1])
    jd$Smiles <- comps[as.character(jd$compound_id)]
    jd$Protein <- targs[as.character(jd$target_id)]

    dbdat <- jd
    colnames(dbdat) <- c('Id', 'xId', 'RFree', 'Ramachandran.Outliers', 'Resolution', 'RMSD_Angles', 'RMSD_bonds', 'lig_confidence', 'CIF', 'Latest.PDB', 'Latest.MTZ', 'Xtal', 'cId', 'tID', 'Smiles', 'Protein')
    gc()
    
    # Main Table Output Handler
    # dbdat <- read.csv(paste(dataDir,'mock.csv', sep='/'), stringsAsFactors=F, row.names=1)
    dedupe <- duplicated(dbdat[,'Xtal'])# duplicated(dbdat[,1])
    if(any(dedupe)) dbdat <- dbdat[!dedupe,]
    #dbdat <- dbdat[,-1]
    dbdat$Decision <- ''
    dbdat$Reason <- ''
    #currentRes <- loadData()
    # Get most Recent Response per xtal
    rownames(dbdat) <- as.character(dbdat$Id)
    if(nrow(response_data) > 0){
        tofill <- as.data.frame(t(sapply(split(response_data, response_data$crystal_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
        rownames(tofill) <- as.character(tofill$crystal_id)
        rninter <- intersect(rownames(tofill), rownames(dbdat))
        dbdat[rninter, 'Decision'] <- unlist(tofill[rninter, 'decision_str'])
        dbdat[rninter, 'Reason'] <- unlist(tofill[rninter, 'reason'])

    }
    rownames(dbdat) <- dbdat[,'Xtal']
    # Sort Data 
    dbdat <- do.call('rbind', 
        lapply(c('Release (notify)', '', 'More Work', 'Release', 'Reject'), function(dec){
            dbdat[ dbdat[ , 'Decision'] == dec , ]
        })
    )
    if(debug) print('Data Loaded')
    inputData <- reactive({dbdat})

    if(debug) print('Data Reactivised')
    r1 <- reactive({
        if(debug) print('Subsetting Table')
        # Subset data
        if(is.null(input$protein) & is.null(input$columns)) inputData()
        else if(is.null(input$columns) & !is.null(input$protein)) inputData()[inputData()$Protein %in% input$protein, ]
        else if(!is.null(input$columns) & is.null(input$protein)) inputData()[ ,input$columns]
        else inputData()[inputData()$Protein %in% input$protein, input$columns]
    })
  
    output$table <- DT::renderDataTable({r1()}, selection = 'single')
    if(debug) print('Response Table Loading')
    # Response Table
    output$resp <- DT::renderDataTable(
        loadData(),
        rownames = FALSE,
        options = list(searching = FALSE, lengthChange = FALSE)
    ) 
    if(debug) print('Response Table Success')
    # NGL viewer side panel stats
    if(debug) print('NGL viewer Tab')
    r2 <- reactive({
        if(debug) print('Update NGL Viewer Side Panel')
        dat <- inputData()[input$Xtal2, ]
        # 
        sn <- names(dat)[3:9]
        val <- dat[3:11]
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

    # Generic Output Messages.
    output$msg <- renderText({'Please click once'})  
    output$msg2 <- renderText({'Please click once'})  
    output$msg3 <- renderText({'If crystal isn\'t showing up, press defaults'})

    sessionGreaterThanMostRecentResponse <- function(id, sessionTime){
        con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
        response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
        dbDisconnect(con)
        mostrecent <- as.data.frame(t(sapply(split(response_data, response_data$crystal_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
        rownames(mostrecent) <- as.character(mostrecent$crystal_id)
        t0 <- mostrecent[as.character(id), 'time_submitted'][[1]]
        output <- ifelse(is.null(t0), TRUE, sessionTime > t0)
        return(output)
    }

    displayModalWhoUpdated <- function(id){
                con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
                response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
                dbDisconnect(con)

                mostrecent <- as.data.frame(t(sapply(split(response_data, response_data$crystal_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
                rownames(mostrecent) <- as.character(mostrecent$crystal_id)
                user <- mostrecent[as.character(id), 'fedid'][[1]]

                showModal(modalDialog(title = "Someone has recently reviewed this crystal", 
                    sprintf("A User (%s) has recently reviewed this structure. Restarting the session to update their response. If you disagree with the current response, please submit another response or select another crystal.", user)
                    , easyClose=TRUE, footer = tagList( modalButton("Cancel"), actionButton("ok", "Restart Session"))
                ))
    }


    # Observers, behaviour will be described as best as possible

    # Upon Row Click
    observeEvent(input$table_rows_selected, {
        if(debug) print('Row Click')
        # Check if Row has been updated since session began, ensure that loadData()[,] # will also get relevant xtal data?
        # Connect to DB and get most recent time...        

        rdat <- r1()[input$table_rows_selected,,drop=FALSE]
        selrow <- rownames(rdat) 
        cId <- dbdat[selrow, 'Id']
        
        #if(sessionTime > max( loadData()[,'timestamp']) ){ 
        if(sessionGreaterThanMostRecentResponse(id=cId, sessionTime=sessionTime)){
            # Update Form window (weird bug with changing decision reupdates form...)
            updateSelectizeInput(session, "Xtal", selected = rownames(rdat), choices = sort(rownames( inputData() )))
            updateSelectizeInput(session, "Xtal2", selected = rownames(rdat), choices = sort(rownames( inputData() )))
            # Move to NGL viewer Page
            updateTabsetPanel(session, "beep", selected = 'NGL Viewer')
        } else {
            displayModalWhoUpdated(id=cId)
        }
    })

    observeEvent(input$ok, {
        if(debug) print('Reload Session')
        session$reload()
    })
  
    resetForm <- function(){
        if(debug) print('Reset Form')
        updateSelectizeInput(session, "Xtal", selected = '', choices = sort(rownames( inputData() )))
        updateSelectizeInput(session, "Xtal2", selected = '', choices = sort(rownames( inputData() )))
        updateTabsetPanel(session, "beep", selected = 'Main Page') 
        session$reload()
    }

    # Upon Main Page Submit
    observeEvent(input$submit, {

        fData <- formData()
        if(debug) print(fData)
        if(any(fData%in%c('', ' '))) {
            showModal(modalDialog(title = "Please fill all fields in the form", 
                "One or more fields have been left empty. Please provide your FedID, a decision and reason(s) before clicking submit."
                , easyClose=TRUE, footer = tagList(modalButton("Cancel"))
            ))
        } else {
             # Get ID...
            cId <- fData[ ,'crystal_id']
            if(sessionGreaterThanMostRecentResponse(id=cId, sessionTime=sessionTime)){
                saveData(fData)
                resetForm()
            } else {
                displayModalWhoUpdated(id=cId)
            }
        }

    })
    
    # Upon NGL Viewer Page Submit
    observeEvent(input$submit2, {

        fData <- formData2()
        if(debug) print(fData)
        if(any(fData%in%c('', ' '))){
            showModal(modalDialog(title = "Please fill all fields in the form", 
                "One or more fields have been left empty. Please provide your FedID, a decision and reason(s) before clicking submit."
                , easyClose=TRUE, footer = tagList(modalButton("Cancel"))
            ))
        } else {
             # Get ID...
            cId <- fData[ ,'crystal_id']
            if(sessionGreaterThanMostRecentResponse(id=cId, sessionTime=sessionTime)){
                saveData(fData)
                resetForm()
            } else {
                displayModalWhoUpdated(id=cId)
            }
        }
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
    # Really need to sort this logic ball out...
    observeEvent(input$Xtal2, {
        # Retry everything to ensure that view loads after stage load...
        session$sendCustomMessage(type="removeAllRepresentations", message=list())
        session$sendCustomMessage(type="setPDB", message=list(''))
        choice = input$Xtal2
        # If pdb is not on pdb... Do things.
        if(debug) message(sprintf("pdb: %s", choice))
        filepath <- dbdat[choice,'Latest.PDB']
        XtalRoot <- try(getRootFP(filepath), silent=T)
        #syscall <- sprintf('cat %s', dir(sprintf('%s/%s',dataDir, choice), pattern = 'pdb', full.names=T))
        syscall <- sprintf('cat %s', filepath)
        if(debug) message(syscall)
        tryAddPDB <- try({
            pdbstrings <- system(syscall, intern = TRUE)
            choice <- paste0(pdbstrings, collapse='\n')
            defaultPdbID <<- choice
            session$sendCustomMessage(type="setPDB2", message=list(choice))
            }, silent = TRUE)
        if(inherits(tryAddPDB, 'try-error')){
            defaultPdbID <<- ''
            session$sendCustomMessage(type="removeAllRepresentations", message=list())
        } else {
            session$sendCustomMessage(type="setPDB2", message=list(choice))
            if(!inherits(XtalRoot, 'try-error')){
                fname <- dir(XtalRoot, pattern = '_event.ccp4', full.names=T)
                #fname <- dir(sprintf('%s/%s',dataDir, choice), pattern = 'ccp4', full.names=T)
                if(debug) message(sprintf('%s: %s', 'eMap:', fname))

                tryAddEvent <- try({
                    event <- readBin(fname, what = 'raw', file.info(fname)$size)
                    event <- base64encode(event, size=NA, endian=.Platform$endian)
                    defaultShell <<- event
                    session$sendCustomMessage(type="addEvent", message=list(event))
                }, silent=T)
                if(inherits(tryAddEvent, 'try-error')){
                    defaultShell <<- ''
                    session$sendCustomMessage(type="removeAllRepresentations", message=list())
                } else {
                    session$sendCustomMessage(type="addEvent", message=list(event))  
                }
            }
        }
    })

    # Go back to main Panel, do a refresh for good measure.
    observeEvent(input$Back, {
        resetForm()
    })

    # When pressed re-create original xtal ngl view...
    observeEvent(input$defaultViewButton, {
        try({session$sendCustomMessage(type="removeAllRepresentations", message=list())}, silent = T)
        try({session$sendCustomMessage(type="setPDB2", message=list(defaultPdbID))}, silent = T)
        try({session$sendCustomMessage(type="addEvent", message=list(defaultShell))}, silent = T)
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
    # Default Coloumn Order... May add more...
    defOrder <- c(
        'Protein', 
        'Smiles',  
        'Decision',
        'Reason',
        'Resolution', 
        'RFree', 
        'lig_confidence', 
        'RMSD_Angles', 
        'RMSD_bonds',  
        'Ramachandran.Outliers', 
        'CIF',
        'Latest.PDB',
        'Latest.MTZ'
    )

    observe({
        updateSelectizeInput(session, 'columns', selected = defOrder, choices = colnames(inputData()))
        updateSelectizeInput(session, 'protein', choices = sort(unique(inputData()$Protein)))
        updateSelectizeInput(session, "Xtal", selected = input$Xtal, choices = sort(rownames( inputData() )))
        updateSelectizeInput(session, "Xtal2", selected = input$Xtal2, choices = sort(rownames( inputData() )))
        #updateSelectizeInput(session, "Xtal2", choices = c('Mpro-x0104', 'Mpro-x0161'))
        updateSelectizeInput(session, 'reason', choices = possRes[[input$decision]])
        updateSelectizeInput(session, 'reason2', choices = possRes[[input$decision2]])
    })
} # Server

#################################################################################
# Runtime
#################################################################################
{
app <- shinyApp(ui = ui, server = server)
cmd <- commandArgs(T)
ip <- cmd[1]
port <- cmd[2]
# ip = '0.0.0.0'
# port = 3838

runApp(app, host=ip, port = as.numeric(port), launch.browser = FALSE)
}
