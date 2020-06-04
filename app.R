# App Refactor
{
debug = TRUE
local = FALSE
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
# Set Path: May need to add something later for files on /dls

# Server Bindings
gpath <- '/srv/shiny-server/'
responsesDir <- '/dls/science/users/mly94721/xchemreview/Responses/' 
dataDir <- '/dls/science/users/mly94721/xchemreview/Data/'
library(devtools)
library(shiny)
library(DT)
library(htmlwidgets)
library(caTools)
library(DBI)

if(local){
 gpath <- '.'
 responsesDir <-file.path(sprintf('%s/%s', gpath, "Responses"))
 source('./db_config.R')
 install.packages('~/nglshiny', repos=NULL, type='source')
 library(nglShiny)
} else {
    install.packages("/dls/science/users/mly94721/xchemreview/nglshiny", repos=NULL, type='source', lib="/dls/science/users/mly94721/R/")
    library(nglShiny, lib.loc = "/dls/science/users/mly94721/R/")
    # Move this to docker...
    install.packages('sendmailR', repos = 'http://cran.rstudio.com/' ,lib ="/dls/science/users/mly94721/R/")
    library(sendmailR, lib.loc = "/dls/science/users/mly94721/R/")
}

# Who doesn't love unclosed loops.
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

defOrder <- c( 
    'Smiles',  
    'Decision',
    'Reason',
    'Resolution', 
    'RFree',
    'Rwork', 
    'lig_confidence', 
    'RMSD_Angles', 
    'RMSD_bonds',  
    'Ramachandran.Outliers',
    'Protein'
)

colss <- c( 
    'Smiles',  
    'Decision',
    'Reason',
    'Resolution', 
    'RFree',
    'Rwork', 
    'lig_confidence', 
    'RMSD_Angles', 
    'RMSD_bonds',  
    'Ramachandran.Outliers', 
    'CIF',
    'Latest.PDB',
    'Latest.MTZ',
    'Protein'
)

if(!local) source('/dls/science/users/mly94721/xchemreview/db_config.R') # Config file...
con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
refinement_data <- dbGetQuery(con, "SELECT id, crystal_name_id, r_free, rcryst, ramachandran_outliers, res, rmsd_angles, rmsd_bonds, lig_confidence_string, cif, pdb_latest, mtz_latest FROM refinement WHERE outcome=4 OR outcome=5")
crystal_data <- dbGetQuery(con, sprintf("SELECT id, crystal_name, compound_id, target_id FROM crystal WHERE id IN (%s)", paste(refinement_data[,'crystal_name_id'], collapse=',')))
target_data <- dbGetQuery(con, sprintf("SELECT * FROM target WHERE id IN (%s)", paste(crystal_data[,'target_id'], collapse=',')))
compound_data <- dbGetQuery(con, sprintf("SELECT * FROM compounds WHERE id IN (%s)", paste(crystal_data[,'compound_id'], collapse=',')))
comps <- compound_data[,2]
names(comps) <- as.character(compound_data[,1])
targs <- target_data[,2]
names(targs) <- as.character(target_data[,1])
# Collapse into DF
jd <- cbind(refinement_data[match(crystal_data[,1], refinement_data[,2]), ], crystal_data[,-1])
jd$Smiles <- comps[as.character(jd$compound_id)]
jd$Protein <- targs[as.character(jd$target_id)]
colnames(jd) <- c('Id', 'xId', 'RFree', 'Rwork', 'Ramachandran.Outliers', 'Resolution', 'RMSD_Angles', 'RMSD_bonds', 'lig_confidence', 'CIF', 'Latest.PDB', 'Latest.MTZ', 'Xtal', 'cId', 'tID', 'Smiles', 'Protein')
response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
dbDisconnect(con)

proteinList <- sort(unique(jd$Protein))
# Source mailing list from file...
source('/dls/science/users/mly94721/xchemreview/mailing_list.R')

xtalList <-  sort(unique(jd[,'Xtal']))

rm(refinement_data, crystal_data, target_data, compound_data, jd, targs, comps)
gc()


defaultRepresentation <- "ball+stick"
defaultColorScheme <- "chainIndex"

possDec <- c("", "Release", "Release (notify)", "More Work", "Reject")
possAns <- possAns2 <- c('Select Decision')

possRes <- tapply(X=response_data$reason, INDEX=response_data$decision_str,
                    function(x){
                        unique(unlist(strsplit(x, '; ')))
                        })

possRes[['Release']] <- c(possRes[['Release']], 'Everything is Wonderful')
possRes[['Release (notify)']] <- c(possRes[['Release (notify)']], 'Alternate binding conformation','Incomplete Density','Weak Density','Low Resolution','Poor Data quality')
possRes[['More Work']] <- c(possRes[['More Work']], 'Cannot View Density', 'Repeat Experiment', 'Check Geometry', 'Check Conformation', 'Check Refinement')
possRes[['Reject']] <- c(possRes[['Reject']], 'Density too weak', 'Insubstantial Evidence','Bad coordination','Incomplete Density')
possDec_int <- 1:4
names(possDec_int) <- c("Release", "Release (notify)", "More Work", "Reject")
}

# UI Code
ui <- navbarPage("XChem Review", id='beep',
	# First Page
	tabPanel('Main',
		fluidRow(
			column(2,
                uiOutput('proteinselect'),
                div(
                    id = "form",
                    textInput("name", "FedID", ""),
                    uiOutput('xtalselect'),
                    #selectizeInput('Xtal', 'Which Structure?', choices = xtalList, multiple = FALSE),
                    selectInput("decision", "Decision", choices = possDec),
                    selectizeInput("reason", "Reason(s)", list(), multiple=TRUE, options(create=TRUE)),
                    textOutput('msg'),
                    actionButton("submit", "Submit", class = "btn-primary"),
                    actionButton('clear', 'Clear', class = 'btn-primary'),
                    #selectInput('protein', 'Select Specific Protein', choices = proteinList, selected= uiOutput("inVar"), multiple=TRUE),
                    selectInput('columns', 'Select Columns to View? (delete/add more values as needed)', choices=colss, selected= defOrder, multiple = TRUE)
                )
			),
			column(8,
				nglShinyOutput('nglShiny', height = '400px'),
                div(style="height: 50px;", br()),
            	DT::dataTableOutput("table")
			),
			column(2,
                textOutput('msg3'),
                #actionButton("fitButton", "Fit"),
                actionButton("defaultViewButton", "View/Update Parameters"),
                actionButton("clearRepresentationsButton", "Clear Representations"),
                #actionButton("updateView", "Update Parameters"),
                hr(),          
                checkboxInput('eventMap', 'Event map', value = TRUE),
                uiOutput('isoEventSlider'), 
                #sliderInput("isoEvent", "Event ISO",
                #    min = 0, max = 10,
                #    value = 1, step = 0.1),
                checkboxInput('twofofcMap', '2fofc map', value = TRUE),
                uiOutput('iso2fofcSlider'),
                #sliderInput("iso2fofc", "2fofc ISO",
                #    min = 0, max = 10,
                #    value = 1.5, step = 0.1),
                checkboxInput('fofcMap', 'fofc Map', value = TRUE),
                uiOutput('isofofcSlider'),
                #sliderInput("isofofc", "fofc ISO",
                #    min = 0, max = 10,
                #    value = 3, step = 0.1),
                hr(),
                numericInput("boxsize", 'Box Size', value = 10, min = 0, max = 100),
                numericInput("clipDist", "Clipping Distance", value=5, min = 0, max = 100),
                sliderInput("fogging", "Fogging:",
                  min = 0, max = 100,
                  value = c(45,58)),
                sliderInput("clipping", "Clipping:",
                  min = 0, max = 100,
                  value = c(47,100)),  
                hr()#,      
                #selectInput("representationSelector", "", nglRepresentations, selected=defaultRepresentation, width=0),
                #selectInput("colorSchemeSelector", "", nglColorSchemes, selected=defaultColorScheme,width=0)
			)
		) # Fluid row
	), # Tab Panel
	tabPanel('Help',
		includeMarkdown(sprintf('%s/%s', gpath, "Pages/include.md"))
	) # Tab Panel
) # Nav Bar Page
# End of UI

#################################################################################
# Server Code
#################################################################################
server <- function(input, output, session) {
    defaultPdbID <- ""
    defaultShell <- ""
    if(debug) message('Server Init')
    session$allowReconnect('force')
	sessionDisconnect <- function() message('User Disconnected')
    session$onSessionEnded(sessionDisconnect)
    sessionTime <- epochTime()
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

    sendEmail <- function(structure, user, decision, reason){
        protein <- gsub('-x[0-9]+', '', structure)
        sendmailR::sendmail(
            from = '<XChemStructureReview@diamond.ac.uk>',
            to = sort(unique(emailListperStructure[[protein]])),#'<tyler.gorrie-stone@diamond.ac.uk>', #emailListperStructure[[structure]],
            subject = sprintf('%s has been labelled as %s', structure, decision),
            msg = sprintf(
'%s has been labelled as %s by %s for the following reason(s): %s.

If you wish to review this change please go to xchemreview.diamond.ac.uk while 
connected to the diamond VPN or via NX.

Direct Link (must be connected to diamond VPN): https://xchemreview.diamond.ac.uk/?xtal=%s&protein=%s

If you disagree with this decision please discuss and change the outcome by submitting a new response.

This email was automatically sent by The XChem Review app

If you believe you have been sent this message in error, please email tyler.gorrie-stone@diamond.ac.uk',
            structure, decision, user, reason, structure, protein),
            control = list(
                smtpServer = 'exchsmtp.stfc.ac.uk',
                smtpPort = 25
                )
        )
    }

    # Save Responses.
    saveData <- function(data, xtaln) {
        fileName <- sprintf("%s_%s.csv",
                            humanTime(),
                            digest::digest(data))
        if(!data[,'fedid'] %in% c('', ' ')){ # Create Modal that prevent empty data from being submitted!!
            write.csv(x = data, file = file.path(responsesDir, fileName),
                  row.names = FALSE, quote = TRUE)
            con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
            dbAppendTable(con, 'review_responses', value = data, row.names=NULL)
            dbDisconnect(con)
            sendEmail(xtaln, data[,'fedid'], data[,'decision_str'], data[,'reason'])
        }
    }

    sessionGreaterThanMostRecentResponse <- function(id, sessionTime){
        con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
        response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
        dbDisconnect(con)
        if(nrow(response_data) > 0){
            mostrecent <- as.data.frame(t(sapply(split(response_data, response_data$crystal_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
            rownames(mostrecent) <- as.character(mostrecent$crystal_id)
            t0 <- mostrecent[as.character(id), 'time_submitted'][[1]]
            output <- ifelse(is.null(t0), TRUE, sessionTime > t0)
        } else {
            output <- TRUE
        }
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

    # Reactives?
    # Form Handler

    fieldsAll <- c("name", 'Xtal', "decision", "reason")
    formData <- reactive({
        data <- sapply(fieldsAll, function(x) paste0(input[[x]], collapse='; '))
        # Get Crystal ID
        xtalname <- data[2]
        data <- c(dbdat[data[2], 'Id'], data[1], possDec_int[data[3]] ,data[3:4], timestamp = epochTime())
        data <- data.frame(t(data), stringsAsFactors=F)
        # Force Coercion
        data[,1] <- as.integer(data[,1])
        data[,3] <- as.integer(data[,3])
        data[,6] <- as.integer(data[,6])
        colnames(data) <- c('crystal_id', 'fedid', 'decision_int', 'decision_str', 'reason', 'time_submitted')
        list(data=data, xtalname=xtalname)
    })

    # Outputs, invest in putting in postgres instead of slurping everything into memory
    # Otherwise extremely slow to shove this in the front end? and on session loading...
    # Can use it for
    if(!local) source('/dls/science/users/mly94721/xchemreview/db_config.R') # Config file...

    getData <- function(db, host_db, db_port, db_user, db_password){
        con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
        refinement_data <- dbGetQuery(con, "SELECT id, crystal_name_id, r_free, rcryst, ramachandran_outliers, res, rmsd_angles, rmsd_bonds, lig_confidence_string, cif, pdb_latest, mtz_latest FROM refinement WHERE outcome=4 OR outcome=5")
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
        colnames(dbdat) <- c('Id', 'xId', 'RFree', 'Rwork', 'Ramachandran.Outliers', 'Resolution', 'RMSD_Angles', 'RMSD_bonds', 'lig_confidence', 'CIF', 'Latest.PDB', 'Latest.MTZ', 'Xtal', 'cId', 'tID', 'Smiles', 'Protein')
        gc()
    
        # Main Table Output Handler
        # dbdat <- read.csv(paste(dataDir,'mock.csv', sep='/'), stringsAsFactors=F, row.names=1)
        dedupe <- duplicated(dbdat[,'Xtal']) # duplicated(dbdat[,1])
        if(any(dedupe)) dbdat <- dbdat[!dedupe,]
        #dbdat <- dbdat[,-1]
        dbdat$Decision <- ''
        dbdat$Reason <- ''

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
        return(dbdat)
    }

    dbdat <- getData(db=db, host_db=host_db, db_port=db_port, 
                    db_user=db_user, db_password=db_password)
    if(debug) print('Data Loaded')

    inputData <- reactive({dbdat})

    # NGL Viewer
    output$nglShiny <- renderNglShiny(
        nglShiny(list(), 300, 300)
    )

    observeEvent(input$decision,{
        if(debug) message('Updating Decision')
        updateSelectizeInput(session, 'reason', choices = possRes[[input$decision]])
    })

    if(debug) print('Data Reactivised')
    r1 <- reactive({
        if(debug) print('Subsetting Table') # Based on input$protein and input$columns?
        # Subset data
        if(is.null(input$protein) & is.null(input$columns)) inputData()
        else if(is.null(input$columns) & !is.null(input$protein)) inputData()[inputData()$Protein %in% input$protein, ]
        else if(!is.null(input$columns) & is.null(input$protein)) inputData()[ ,input$columns]
        else inputData()[inputData()$Protein %in% input$protein, input$columns]
    })
  
    output$table <- DT::renderDataTable({r1()}, options= list(pageLength=20),selection = 'single')

    # Generic Output Messages.
    output$msg <- renderText({'Please click once'})  
    output$msg3 <- renderText({'NGL Viewer Controls'})

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
        } else {
            displayModalWhoUpdated(id=cId)
        }
    })

    observeEvent(input$ok, {
        if(debug) print('Reload Session')
        session$reload()
    })
  
    observeEvent(input$updateView,{
        session$sendCustomMessage(type="fit", message=list())
    })

    resetForm <- function(){
        if(debug) print('Reset Form')
        updateSelectizeInput(session, "Xtal", selected = '', choices = sort(rownames( inputData() )))
        session$reload()
        #dbdat <- getData(db=db, host_db=host_db, db_port=db_port, 
        #                db_user=db_user, db_password=db_password)
        #inputData <- reactive({dbdat})
        #sessionTime <- epochTime()
    }

    # Upon Main Page Submit
    observeEvent(input$submit, {

        fData <- formData()[[1]]
        xtaln <- formData()[[2]]
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
                saveData(fData, xtaln)
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

    # Upon pressing Fit, Fit structure in window
    observeEvent(input$fitButton, {
        session$sendCustomMessage(type="fit", message=list())
    })

    observeEvent(input$updateView,{
        session$sendCustomMessage(type="updateParams", message=list(input$clipDist, 
            input$clipping[1], input$clipping[2], input$fogging[1], input$fogging[2]))
    })

    # Upon pressing Clear, Remove Everything
    observeEvent(input$clearRepresentationsButton, {
        session$sendCustomMessage(type="removeAllRepresentations", message=list())
        #updateSelectInput(session, "representationSelector", label=NULL, choices=NULL,  selected=defaultRepresentation)
        #updateSelectInput(session, "colorSchemeSelector", label=NULL, choices=NULL,  selected=defaultColorScheme)
    })
  
    # Load structure and event to NGL stage!
    uploadPDB <- function(filepath, input){
        syscall <- sprintf('cat %s', filepath)
        if(debug) message(syscall)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse='\n')
        session$sendCustomMessage(
            type="setPDB2", 
            message=list(choice, 
                as.character(input$clipDist), 
                as.character(input$clipping[1]),
                as.character(input$clipping[2]),
                as.character(input$fogging[1]),
                as.character(input$fogging[2])
            )
        )
    }

    getExt <- function(x) sapply(strsplit(x, '[.]'), tail, 1)

    findFirstMatchingFile <- function(x, fp){
        results <- lapply(x, function(y, fp){
            dir(fp, pattern=y, full.names=T)
        }, fp=fp)
        filesFound <- sapply(results, length)>0
        if(any(filesFound)){
            # If any files are found, take the first one, otherwise return NA...
            first <- which(sapply(results, length)>0)[1] # Take First one...
            output <- results[[first]][1]
        } else {
            output <- NA
        }
        return(output)
    }

    findFiles <- function(fp){
        eventmapStrings <- c('_event.ccp4', 'event_map.map')
        fofc2Strings <- c('_2fofc.cpp4', '^2fofc.map')
        fofcStrings <- c('_fofc.ccp4', '^fofc.map')

        # first look for .ccp4 files
        eventmapFile <- findFirstMatchingFile(eventmapStrings, fp=fp)
        fofc2File <- findFirstMatchingFile(fofc2Strings, fp=fp)
        fofcFile <- findFirstMatchingFile(fofcStrings, fp=fp)
        return(c(eventmapFile, fofc2File, fofcFile))
    }

    uploadEMaps <- function(XtalRoot, input){
        if(debug) print(dir(XtalRoot))
        theFiles <- findFiles(XtalRoot) # row 1 is: 1 = event map, 2 = 2fofc and 3 = fofc
        #fname <- dir(XtalRoot, pattern = '_event.ccp4', full.names=T)[1]
        #fname <- dir(XtalRoot, pattern = 'event', full.names=T)[1]
        fname <- theFiles[1]
        if(debug) message(sprintf('%s: %s', 'event Map', fname))

        if(input$isoEvent){   
            event <- readBin(fname, what = 'raw', file.info(fname)$size)
            event <- base64encode(event, size=NA, endian=.Platform$endian)
            # addEvent requires:
            # filepath as event blob (base64string)
            # desired iso level

            session$sendCustomMessage(type="addEvent", 
                message=list(
                    event, 
                    as.character(input$isoEvent), 
                    as.character('orange'), 
                    as.character('false'), 
                    as.character(getExt(fname)),
                    as.character(input$boxsize)
                )
            )
        }

        if(input$iso2fofc){
            #fname <- dir(XtalRoot, pattern = '_2fofc.ccp4', full.names=T)
            #fname <- dir(XtalRoot, pattern = '2fofc.map', full.names=T)[1]
            fname <- theFiles[2]
            if(debug) message(sprintf('%s: %s', '2fofc', fname))
            event <- readBin(fname, what = 'raw', file.info(fname)$size)
            event <- base64encode(event, size=NA, endian=.Platform$endian)
            session$sendCustomMessage(type="addEvent", 
                message=list(
                    event, 
                    as.character(input$iso2fofc), 
                    as.character('blue'), 
                    as.character('false'), 
                    as.character(getExt(fname)),
                    as.character(input$boxsize)
                )
            )
        }
        if(input$isofofc){
            #fname <- dir(XtalRoot, pattern = '_fofc.ccp4', full.names=T)[1]
            #fname <- dir(XtalRoot, pattern = '^fofc.map', full.names=T)[1]
            fname <- theFiles[3]
            if(debug) message(sprintf('%s: %s', 'fofc', fname))
            event <- readBin(fname, what = 'raw', file.info(fname)$size)
            event <- base64encode(event, size=NA, endian=.Platform$endian)
            session$sendCustomMessage(type="addEvent", 
                message=list(
                    event, 
                    as.character(input$isofofc), 
                    as.character('lightgreen'), 
                    as.character('false'), 
                    as.character(getExt(fname)),
                    as.character(input$boxsize)
                )
            )
            session$sendCustomMessage(type="addEvent", 
                message=list(
                    event, 
                    as.character(input$isofofc), 
                    as.character('tomato'), 
                    as.character('true'), 
                    as.character(getExt(fname)),
                    as.character(input$boxsize)
                )
            )
        }
    }

    # Really need to sort this logic ball out...
    observeEvent(input$Xtal, {
        # Retry everything to ensure that view loads after stage load...
        choice = input$Xtal
        filepath <- dbdat[choice,'Latest.PDB']
        XtalRoot <- try(getRootFP(filepath), silent=T)
        defaultPdbID <- filepath
        defaultShell <- XtalRoot
        tryAddPDB <- try(uploadPDB(filepath=defaultPdbID, input=input), silent=T)
        if(inherits(tryAddPDB, 'try-error')){
            defaultPdbID <- ''
            defaultShell <- ''
            session$sendCustomMessage(type="removeAllRepresentations", message=list())
        } else {
            if(!inherits(XtalRoot, 'try-error')){
                tryAddEvent <- try(uploadEMaps(XtalRoot=defaultShell, input=input), silent=T)
                if(inherits(tryAddEvent, 'try-error')){
                    defaultShell <- ''
                    session$sendCustomMessage(type="removeAllRepresentations", message=list())
                }
            }
        }
    })

    # Go back to main Panel, do a refresh for good measure.
    observeEvent(input$clear, {
        resetForm()
    })

    # When pressed re-create original xtal ngl view...
    observeEvent(input$defaultViewButton, {
        #try({session$sendCustomMessage(type="removeAllComponents", message=list())}, silent=T)
        choice = input$Xtal
        filepath <- dbdat[choice,'Latest.PDB']
        XtalRoot <- try(getRootFP(filepath), silent=T)
        defaultPdbID <- filepath
        defaultShell <- XtalRoot
        tryAddPDB <- try(uploadPDB(filepath=defaultPdbID, input=input), silent=T)
        if(inherits(tryAddPDB, 'try-error')){
            defaultPdbID <- ''
            defaultShell <- ''
            session$sendCustomMessage(type="removeAllRepresentations", message=list())
        } else {
            if(!inherits(XtalRoot, 'try-error')){
                tryAddEvent <- try(uploadEMaps(XtalRoot=defaultShell, input=input), silent=T)
                if(inherits(tryAddEvent, 'try-error')){
                    defaultShell <- ''
                    session$sendCustomMessage(type="removeAllRepresentations", message=list())
                }
            }
        }
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

    # Ticker
    autoInvalidate <- reactiveTimer(10000)
  	observe({
    	autoInvalidate()
    	cat("")
  	})

    output$proteinselect <- renderUI({
        query <- parseQueryString(session$clientData$url_search)
        if (!is.null(query[['protein']])) {
            selectInput('protein', 'Select Specific Protein', choices = proteinList, selected=query[['protein']], multiple=TRUE) 
        } else {
            selectInput('protein', 'Select Specific Protein', choices = proteinList, selected=list(), multiple=TRUE)                  
        }
    })

    output$xtalselect <- renderUI({
        query <- parseQueryString(session$clientData$url_search)
        if(!is.null(query[['xtal']])){
            selectizeInput('Xtal', 'Which Structure?', choices = xtalList, selected = query[['xtal']], multiple = FALSE)
        } else {
            selectizeInput('Xtal', 'Which Structure?', choices = xtalList, multiple = FALSE)
        }
    })

    output$isoEventSlider <- renderUI({
        if(isoEvent){
        sliderInput("isoEvent", "",
                    min = 0, max = 10,
                    value = 1, step = 0.1)
        } else {
            NULL
        }
    })

    output$iso2fofcSlider <- renderUI({
        if(iso2fofc){
        sliderInput("iso2fofc", "",
                    min = 0, max = 10,
                    value = 1.5, step = 0.1)
        } else {
            NULL
        }
    })

    output$isofofcSlider <- renderUI({
        if(isofofc){
            sliderInput("isofofc", "",
                min = 0, max = 10,
                value = 3, step = 0.1)
        } else {
            NULL
        }
    })

} # Server

#################################################################################
# Runtime
#################################################################################
{
app <- shinyApp(ui = ui, server = server)
if(local){
ip <- '0.0.0.0'
port <- '3838'
} else {
cmd <- commandArgs(T)
ip <- cmd[1]
port <- cmd[2]
}
runApp(app, host=ip, port = as.numeric(port), launch.browser = FALSE)
}
