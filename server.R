# Probably should break this up too...
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

    sendEmail <- function(structure, user, decision, reason, comments){
        protein <- gsub('-x[0-9]+', '', structure)
        sendmailR::sendmail(
            from = '<XChemStructureReview@diamond.ac.uk>',
            to = sort(unique(emailListperStructure[[protein]])),#'<tyler.gorrie-stone@diamond.ac.uk>', #emailListperStructure[[structure]],
            subject = sprintf('%s has been labelled as %s', structure, decision),
            msg = sprintf(
'%s has been labelled as %s by %s for the following reason(s): %s.

With these additional comments:

%s

-------------------------------
If you wish to review this change please go to xchemreview.diamond.ac.uk while 
connected to the diamond VPN or via NX.

Direct Link (must be connected to diamond VPN): https://xchemreview.diamond.ac.uk/?xtal=%s&protein=%s

If you disagree with this decision please discuss and change the outcome by submitting a new response.

This email was automatically sent by The XChem Review app

If you believe you have been sent this message in error, please email tyler.gorrie-stone@diamond.ac.uk',
            structure, decision, user, reason, comments, structure, protein),
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
            sendEmail(xtaln, data[,'fedid'], data[,'decision_str'], data[,'reason'], input$comments)
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
        data <- c(dbdat[data[2], 'xId'], data[1], possDec_int[data[3]] ,data[3:4], timestamp = epochTime())
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
        rownames(dbdat) <- as.character(dbdat$xId)
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

    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
    dbDisconnect(con)

    possRes <- tapply(X=response_data$reason, INDEX=response_data$decision_str,
                        function(x){
                            unique(unlist(strsplit(x, '; ')))
                        })
    possRes[['Release']] <- unique(c(possRes[['Release']], 'Everything is Wonderful'))
    possRes[['Release (notify)']] <- unique(c(possRes[['Release (notify)']], 'Alternate binding conformation','Incomplete Density','Weak Density','Low Resolution','Poor Data quality'))
    possRes[['More Work']] <- unique(c(possRes[['More Work']], 'Cannot View Density', 'Repeat Experiment', 'Check Geometry', 'Check Conformation', 'Check Refinement'))
    possRes[['Reject']] <- unique(c(possRes[['Reject']], 'Density too weak', 'Insubstantial Evidence','Bad coordination','Incomplete Density'))
    possDec_int <- 1:4
    names(possDec_int) <- c("Release", "Release (notify)", "More Work", "Reject")


    inputData <- reactive({dbdat})

    # NGL Viewer
    output$nglShiny <- renderNglShiny(
        nglShiny(list(), width=NULL, height=100)
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
  
    output$table <- DT::renderDataTable({r1()},
                                        selection = 'single', 
                                        options = list(
                                            pageLength = 20, 
                                            drawCallback = I("function( settings ) {document.getElementById('table').style.width = '600px';}")
                                            )
                                        )

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
        xId <- dbdat[selrow, 'xId']        
        #if(sessionTime > max( loadData()[,'timestamp']) ){ 
        if(sessionGreaterThanMostRecentResponse(id=xId, sessionTime=sessionTime)){
            # Update Form window (weird bug with changing decision reupdates form...)
            updateSelectizeInput(session, "Xtal", selected = rownames(rdat), choices = sort(rownames( inputData() )))
        } else {
            displayModalWhoUpdated(id=xId)
        }
    })

    observeEvent(input$ok, {
        if(debug) print('Reload Session')
        session$reload()
    })
  
    observeEvent(input$updateView,{
        session$sendCustomMessage(type="updateParams", message=list())
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
            xId <- fData[ ,'crystal_id']
            if(sessionGreaterThanMostRecentResponse(id=xId, sessionTime=sessionTime)){
                saveData(fData, xtaln)
                resetForm()
            } else {
                displayModalWhoUpdated(id=xId)
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
        session$sendCustomMessage(type="ligfit", message=list())
    })

    #observeEvent(input$updateView,{
    #    session$sendCustomMessage(type="updateParams", message=list(input$clipDist, 
    #        input$clipping[1], input$clipping[2], input$fogging[1], input$fogging[2]))
    #})

    observeEvent(input$clipping, {
        session$sendCustomMessage(type="updateParams", message=list(input$clipDist, 
            input$clipping[1], input$clipping[2], input$fogging[1], input$fogging[2]))
    })

    observeEvent(input$fogging, {
        session$sendCustomMessage(type="updateParams", message=list(input$clipDist, 
            input$clipping[1], input$clipping[2], input$fogging[1], input$fogging[2]))
    })

    observeEvent(input$clipDist, {
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
#        withProgress(message = 'Uploading PDB', style='notification', detail = 'Finding File', value = 0, {
            syscall <- sprintf('cat %s', filepath)
            if(debug) message(syscall)
#            incProgress(.25, detail = 'Reading File')
            pdbstrings <- system(syscall, intern = TRUE)
            choice <- paste0(pdbstrings, collapse='\n')
#            incProgress(.25, detail = 'Sending File to stage')
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
#            setProgress(1)
#        })
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

    updateVisabilities <- function(event, twofofc, fofc){
        session$sendCustomMessage(type='updateVisabilities',
            list(
                tolower(as.character(as.logical(event))),
                tolower(as.character(as.logical(twofofc))),
                tolower(as.character(as.logical(fofc))),
                tolower(as.character(as.logical(fofc)))
                )
        )
    }

    uploadEMaps <- function(XtalRoot, input){
#        withProgress(message = 'Loading maps', detail = 'Finding Files', style='notification', value=0, {
            if(debug) print(dir(XtalRoot))
            theFiles <- findFiles(XtalRoot) # row 1 is: 1 = event map, 2 = 2fofc and 3 = fofc
#            incProgress(.1, details='Load Event Map')
            #fname <- dir(XtalRoot, pattern = '_event.ccp4', full.names=T)[1]
            #fname <- dir(XtalRoot, pattern = 'event', full.names=T)[1]
            fname <- theFiles[1]
            if(debug) message(sprintf('%s: %s', 'event Map', fname))

            if(TRUE){   
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

            updateVisabilities(event=input$eventMap, twofofc=input$twofofcMap, fofc=input$fofcMap)

#            incProgress(.3, details='Load 2fofc Map')
            if(TRUE){
                #fname <- dir(XtalRoot, pattern = '_2fofc.ccp4', full.names=T)
                #fname <- dir(XtalRoot, pattern = '2fofc.map', full.names=T)[1]
                fname <- theFiles[2]
                if(debug) message(sprintf('%s: %s', '2fofc', fname))
                event <- readBin(fname, what = 'raw', file.info(fname)$size)
                event <- base64encode(event, size=NA, endian=.Platform$endian)
                session$sendCustomMessage(type="add2fofc", 
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
            updateVisabilities(event=input$eventMap, twofofc=input$twofofcMap, fofc=input$fofcMap)  
#            incProgress(.3, details='Load fofc Map')
            if(TRUE){
                #fname <- dir(XtalRoot, pattern = '_fofc.ccp4', full.names=T)[1]
                #fname <- dir(XtalRoot, pattern = '^fofc.map', full.names=T)[1]
                fname <- theFiles[3]
                if(debug) message(sprintf('%s: %s', 'fofc', fname))
                event <- readBin(fname, what = 'raw', file.info(fname)$size)
                event <- base64encode(event, size=NA, endian=.Platform$endian)
                session$sendCustomMessage(type="addfofc", 
                    message=list(
                        event, 
                        as.character(input$isofofc), 
                        as.character('lightgreen'), 
                        as.character('false'), 
                        as.character(getExt(fname)),
                        as.character(input$boxsize)
                    )
                )
                session$sendCustomMessage(type="addfofc_negative", 
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
            updateVisabilities(event=input$eventMap, twofofc=input$twofofcMap, fofc=input$fofcMap)
    }



    observeEvent(input$eventMap, {
        updateVisabilities(event=input$eventMap, twofofc=input$twofofcMap, fofc=input$fofcMap)
    })
    observeEvent(input$twofofcMap, {
        updateVisabilities(event=input$eventMap, twofofc=input$twofofcMap, fofc=input$fofcMap)
    })
    observeEvent(input$fofcMap, {
        updateVisabilities(event=input$eventMap, twofofc=input$twofofcMap, fofc=input$fofcMap)
    })

    observeEvent(input$boxsize,{
        if(input$eventMap)  session$sendCustomMessage(type='twiddleEvent',          list(as.character(input$isoEvent),as.character(input$boxsize)))
        if(input$twofofcMap)session$sendCustomMessage(type='twiddle2fofc',          list(as.character(input$iso2fofc),as.character(input$boxsize)))
        if(input$fofcMap)   session$sendCustomMessage(type='twiddlefofc',           list(as.character(input$isofofc),as.character(input$boxsize)))
        if(input$fofcMap)   session$sendCustomMessage(type='twiddlefofc_negative',  list(as.character(input$isofofc),as.character(input$boxsize)))
    })

    observeEvent(input$isoEvent,{
        session$sendCustomMessage(type='twiddleEvent',          list(as.character(input$isoEvent),as.character(input$boxsize)))
    })

    observeEvent(input$iso2fofc,{
        session$sendCustomMessage(type='twiddle2fofc',          list(as.character(input$iso2fofc),as.character(input$boxsize)))
    })

    observeEvent(input$isofofc,{
        session$sendCustomMessage(type='twiddlefofc',           list(as.character(input$isofofc),as.character(input$boxsize)))
        session$sendCustomMessage(type='twiddlefofc_negative',  list(as.character(input$isofofc),as.character(input$boxsize)))
    })

    # Really need to sort this logic ball out...
    observeEvent(input$Xtal, {
#        withProgress(message = 'Loading Crystal', style='notification', value=.1,{
            # Retry everything to ensure that view loads after stage load...
            choice = input$Xtal
            filepath <- dbdat[choice,'Latest.PDB']
            XtalRoot <- try(getRootFP(filepath), silent=T)
            defaultPdbID <- filepath
            defaultShell <- XtalRoot

            spfile <- tail(dir(XtalRoot, pattern='A-1101.png', full.names=T, rec=T),1)
            output$spiderPlot <- renderImage({
                if(length(spfile) == 1){
                    list(src = spfile,
                    contentType = 'image/png',
                    width=200,
                    height=200)
                } else { 
                    list(src = '',
                    contentType = 'image/png',
                    width=200,
                    height=200)
                }
            }, deleteFile=FALSE)

            tryAddPDB <- try(uploadPDB(filepath=defaultPdbID, input=input), silent=T)
#            incProgress(.5, detail = 'Attempting to load maps')
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
#            setProgress(1)
#        })

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

    #output$spiderplot <- renderUI({
    #    
    #    })

    output$xtalselect <- renderUI({
        query <- parseQueryString(session$clientData$url_search)
        if(!is.null(query[['xtal']])){
            selectizeInput('Xtal', 'Which Structure?', choices = xtalList, selected = query[['xtal']], multiple = FALSE)
        } else {
            selectizeInput('Xtal', 'Which Structure?', choices = xtalList, selected = list(), multiple = FALSE)
        }
    })

    output$isoEventSlider <- renderUI({
        if(input$eventMap){
            sliderInput("isoEvent", "",
                    min = 0, max = 10,
                    value = 1, step = 0.1)
        } else {
            NULL
        }
    })

    output$iso2fofcSlider <- renderUI({
        if(input$twofofcMap){
            sliderInput("iso2fofc", "",
                    min = 0, max = 10,
                    value = 1.5, step = 0.1)
        } else {
            NULL
        }
    })

    output$isofofcSlider <- renderUI({
        if(input$fofcMap){
            sliderInput("isofofc", "",
                min = 0, max = 10,
                value = 3, step = 0.1)
        } else {
            NULL
        }
    })

    controlPanel = TRUE
    output$controlRow <- renderUI({
        if(TRUE){
            fluidRow(
                    column(6, numericInput("boxsize", 'Box Size', value = 10, min = 0, max = 100, width='100px')),
                    column(6, numericInput("clipDist", "Clip Dist", value=5, min = 0, max = 100, width='100px'))
                    )
        } else {
            NULL
        }
    })

    output$controlFog <- renderUI({
        if(controlPanel){
            sliderInput("fogging", "Fogging:", min = 0, max = 100, value = c(45,58))
        } else {
            NULL
        }
    })


    output$controlClip <- renderUI({
        if(controlPanel){
            sliderInput("clipping", "Clipping:", min = 0, max = 100, value = c(47,100))
        } else {
            NULL
        }
    })


} # Server


                    