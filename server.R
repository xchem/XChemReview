# Probably should break this up too...
#################################################################################
# Server Code
#################################################################################
server <- function(input, output, session) {
    sID <- sample(1:100000, 1)
    debugMessage <- function(sID, text){
        message(sprintf('sid: %s | %s | %s', sID, text, Sys.time()))
    }
    defaultPdbID <- ""
    defaultShell <- ""
    if(debug) debugMessage(sID=sID, sprintf('Session init'))
    session$allowReconnect(FALSE)
	sessionDisconnect <- function() debugMessage(sID=sID, 'Disconnected')
    session$onSessionEnded(sessionDisconnect)
    sessionTime <- epochTime()
    options <- list(pdbID="")

    # Page 1

    # Functions
    # Read-in Responses Form Data
    loadData <- function() {
        files <- list.files(file.path(responsesDir), full.names = TRUE)
        data <- lapply(files, read.csv, stringsAsFactors = FALSE)
        data <- do.call(rbind, data)
        # Now only contains ID instead of xtal name...
        data
    }

    slackFilters <- c('has joined the channel')
    getChannelList <- function(){
        channellist <- content(POST(url='https://slack.com/api/conversations.list', body=list(token=api)))
        channels <- sapply(c('name', 'id'), function(x) sapply(channellist[[2]], '[[', x))
        rownames(channels) <- channels[,1]
        channels <- channels[!rownames(channels)%in%c('welcome', 'team', 'project'),]
        return(channels)
    }

    getChannel <- function(structure, channels){
        channelname <- tolower(gsub('[^[:alnum:]]', '', structure))
        ifelse(channelname %in% rownames(channels), channels[channelname, 2], NA)
    }

    parseConversation <- function(channel){
        history <- content(POST(url='https://slack.com/api/conversations.history', 
                body = list(token = api, channel = channel)))
        convoblock <- do.call('rbind', parseMessageContent(conversations_history = history[[2]], channel=channel))
        return(convoblock[nrow(convoblock):1,])
    }

    parseMessageContent <- function(conversations_history, channel){
        messageContent <- lapply(conversations_history, parseIndividualMessage, channel=channel)
        return(messageContent)
    }

    parseReply <- function(x, content, ts){
        user <- content[[x]]$user
        text <- content[[x]]$text
        out <- c(user, text, ts)
        return(out)
    }

    parseThread <- function(channel, ts, api){
        thread <- content(POST(url='https://slack.com/api/conversations.replies',
                            body = list(token = api, channel = channel, ts = ts)))
        content <- thread[[1]]
        contentid <- length(content):1
        out <- t(sapply(contentid, parseReply, content=content, ts=ts))
        return(out)
    }

    parseIndividualMessage <- function(message, channel){
        user <- message$user
        text <- message$text
        ts <- message$thread_ts
        broadcastedreply <- any(grepl('root', names(message)))
        if(broadcastedreply) return(NULL)
        if(!is.null(ts)){
            out <- parseThread(channel = channel, ts = ts, api = api)
        } else {
            out <- t(c(user, text, message$ts))
        }
        return(out)
    }

    userList <- function(){
        users <- content(POST(url='https://slack.com/api/users.list', 
            body = list(token = api)))
        userlist <- t(sapply(users$members, function(x){
            pull <- c(x$id, x$real_name)
            if(length(pull)<2) pull <- c(pull, NA)
            return(pull)
        }))
        users <- userlist[,2]
        names(users) <- userlist[,1]
        return(users)
    }

    refreshChat <- function(channel){
        users <- userList()
        convo <- try(parseConversation(channel=channel), silent=T)
        if(!inherits(convo, 'try-error')){
            convo <- convo[!convo[,2] == '', ,drop=F]
            convo <- convo[!grepl(slackFilters, convo[,2]), ,drop=F]
            convo[,1] <- users[convo[,1]]
            stamps <- duplicated(convo[,3])
            textdump <- paste(mapply(X=1:nrow(convo), Y=stamps, function(X,Y){
                date <- as.POSIXct(as.numeric(convo[X, 3]), origin="1970-01-01")
                if(Y) sprintf('\t- <%s>: %s', convo[X,1], convo[X,2])
                else sprintf('[%s] - <%s>: %s ',  date, convo[X,1], convo[X,2])
            }), collapse='\n')
            output$chat <- renderText({textdump})
        } else {
            output$chat <- renderText({'Select A Channel'})
        }
    }

    createChannel <- function(structure){
        # lowercase and despecial
        channelname <- tolower(gsub('[^[:alnum:]]', '', structure))
        resp <- POST(url='https://slack.com/api/conversations.create', body=list(token=api, name =channelname))
        return(content(resp)$channel$id)
    }

    sendMessageToSlack <- function(channel, message, name){
        content(POST(url='https://slack.com/api/chat.postMessage', 
            body=list(token=apiuser, 
                channel=channel, 
                text= message, 
                as_user='true', 
                username=name)))
    }

    sendEmail <- function(structure, user, decision, reason, comments){
        if(debug) debugMessage(sID=sID, sprintf('Communicating to Slack...'))
        channels <- getChannelList()
        channelID <- getChannel(structure, channels)
        if(is.na(channelID)){
            # Create Channel, then get ID immediately
            channelID <- createChannel(structure=structure)
        }
        # Post comment...
        sendMessageToSlack(channel=channelID, 
                            message= sprintf('%s has been labelled as %s by %s for the following reason(s): %s.

With these additional comments:

%s', structure, decision, user, reason, comments), 
                            name = 'xchemreview-bot')

        if(debug) debugMessage(sID=sID, sprintf('Sending Email'))
        protein <- gsub('-x[0-9]+', '', structure)
        sendmailR::sendmail(
            from = '<XChemStructureReview@diamond.ac.uk>',
            to = sort(unique(emailListperStructure[[protein]])),#'<tyler.gorrie-stone@diamond.ac.uk>', #emailListperStructure[[structure]],
            #to = '<tyler.gorrie-stone@diamond.ac.uk>',
            subject = sprintf('%s has been labelled as %s', structure, decision),
            msg = sprintf(
'%s has been labelled as %s by %s for the following reason(s): %s.

With these additional comments:

%s

-------------------------------
If you wish to review this change please go to xchemreview.diamond.ac.uk while 
connected to the diamond VPN or via NX.

Direct Link (must be connected to diamond VPN): https://xchemreview.diamond.ac.uk/?xtal=%s&protein=%s

If you disagree with this decision please discuss at the in the slack channel: (https://xchemreview.slack.com/archives/%s) or change the outcome by submitting a new response.

This email was automatically sent by The XChem Review app

If you have trouble joining the slack channel please use this invitation link: https://join.slack.com/t/xchemreview/shared_invite/zt-fpocaf6e-JQp~U6rcbGrre33E~7~faw

If you believe you have been sent this message in error, please email tyler.gorrie-stone@diamond.ac.uk',
            structure, decision, user, reason, comments, structure, protein, channelID),
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
        refinement_data <- dbGetQuery(con, "SELECT id, crystal_name_id, r_free, rcryst, ramachandran_outliers, res, rmsd_angles, rmsd_bonds, lig_confidence_string, spacegroup, outcome, cif, pdb_latest, mtz_latest FROM refinement WHERE outcome=4 OR outcome=5 OR outcome=6")
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
        colnames(dbdat) <- c('Id', 'xId', 'RFree', 'Rwork', 'Ramachandran.Outliers', 'Resolution', 'RMSD_Angles', 
            'RMSD_bonds', 'lig_confidence', 'Space_Group', 'XCEoutcome', 'CIF', 'Latest.PDB', 'Latest.MTZ', 'Xtal', 
            'cId', 'tID', 'Smiles', 'Protein')
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
        lapply(c("","Release", "More Refinement", "More Experiments", "Reject"), function(dec){
            dbdat[ dbdat[ , 'Decision'] == dec , ]
        })
        )
        return(dbdat)
    }

    dbdat <- getData(db=db, host_db=host_db, db_port=db_port, 
                    db_user=db_user, db_password=db_password)
    if(debug) debugMessage(sID=sID, sprintf('Data Loaded'))

    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
    dbDisconnect(con)

    possDec <- c("", "Release", "More Refinement", "More Experiments", "Reject")
    possAns <- possAns2 <- c('Select Decision')
    possRes <- list()
    possRes[['Release']] <- c('High Confidence', 'Clear Density, Unexpected Ligand', 'Correct Ligand, Weak Density', 'Low Confidence', 'No Ligand Present')
    possRes[["More Refinement"]] <- c('Check ligand conformation','Check sidechain rotamers','Check Rfactors','Check that refinement converged','Improve water model','Build alternate conformations','Fix geometry','Trim down ligand','Density did not load','Other')
    possRes[["More Experiments"]] <- c('Get better density','Get better resolution','Confirm ligand identity','Check if ligand changed','Other')
    possRes[["Reject"]] <- c('Density not convincing','Too few interactions','Binding site too noisy','Other')
    possDec_int <- 1:4
    names(possDec_int) <- c("Release", "More Refinement", "More Experiments", "Reject")


    inputData <- reactive({dbdat})
    output$progtext <- renderText({'Table Loaded'}) 

    # NGL Viewer
    output$nglShiny <- renderNglShiny(
        nglShiny(name = 'nglShiny', list(), width=NULL, height=100)
    )

    # Control Panel listener
    callModule(modalModule, "foo")

    # Update behaviour for these...
    observeEvent(input$clickedAtoms, {
        print(input$clickedAtoms)
    })

    output$selAtoms <- renderText({'Selected Atoms'})

    observeEvent(input$clickNames, {
        print(input$clickNames)
        output$selAtoms <- renderText({paste(input$clickNames, collapse = ';')})
    })

    observeEvent(input$decision,{
        if(debug) debugMessage(sID=sID, sprintf('Update Decision'))
        updateSelectizeInput(session, 'reason', choices = possRes[[input$decision]])
    })

    if(debug) debugMessage(sID=sID, sprintf('Data Reactivised'))

    r1 <- reactive({
        rowidx <- rep(FALSE, nrow(inputData()))
        outcome <- inputData()$XCEoutcome
        if(input$out4) rowidx[outcome==4] <- TRUE
        if(input$out5) rowidx[outcome==5] <- TRUE
        if(input$out6) rowidx[outcome==6] <- TRUE
        if(debug) debugMessage(sID=sID, sprintf('Subsetting Table')) # Based on input$protein and input$columns?
        # Subset data
        if(is.null(input$protein) & is.null(input$columns)) inputData()[rowidx, ]
        else if(is.null(input$columns) & !is.null(input$protein)) inputData()[rowidx & inputData()$Protein %in% input$protein, ]
        else if(!is.null(input$columns) & is.null(input$protein)) inputData()[rowidx,input$columns]
        else inputData()[rowidx & inputData()$Protein %in% input$protein, input$columns]
    })
  
    output$table <- DT::renderDataTable(
        #{r1()},
        {datatable(r1(), selection = 'single', options = list(
            pageLength = 20
        )) %>% formatStyle(
        'Decision',
        target = 'row',
        backgroundColor = styleEqual(c('Release', 'More Refinement', 'More Experiments', 'Reject'), c('#648FFF', '#FFB000', '#FE6100', '#DC267F'))
        )}  
    )

    # Generic Output Messages.
    output$msg <- renderText({'Please click once'}) 
    output$missingFiles <- renderText({''}) 
    output$msg3 <- renderText({'NGL Viewer Controls'})
    output$progtext <- renderText({''}) # User Feedback...
    # Observers, behaviour will be described as best as possible
    # Upon Row Click
    observeEvent(input$table_rows_selected, {
        # Check if Row has been updated since session began, ensure that loadData()[,] # will also get relevant xtal data?
        # Connect to DB and get most recent time...        
        rdat <- r1()[input$table_rows_selected,,drop=FALSE]
        selrow <- rownames(rdat)
        if(debug) debugMessage(sID=sID, sprintf('Selecting: %s', selrow)) 
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
        if(debug) debugMessage(sID=sID, sprintf('Reloading Session'))
        session$reload()
    })
  
    observeEvent(input$updateView,{
        session$sendCustomMessage(type="updateParams", message=list())
    })

    resetForm <- function(){
        if(debug) debugMessage(sID=sID, sprintf('Resetting Form'))
        updateSelectizeInput(session, "Xtal", selected = '', choices = sort(rownames( inputData() )))
        session$reload()
    }

    # Upon Main Page Submit
    observeEvent(input$submit, {
        fData <- formData()[[1]]
        xtaln <- formData()[[2]]
        if(debug) debugMessage(sID=sID, sprintf('Submitting Form'))
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
    })
  
    # Load structure and event to NGL stage!
    uploadPDB <- function(filepath, input){
#        withProgress(message = 'Uploading PDB', style='notification', detail = 'Finding File', value = 0, {
            syscall <- sprintf('cat %s', filepath)
            if(debug) debugMessage(sID=sID, sprintf('Executing: %s', syscall))
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
        eventmapStrings <- c('_event_cut.ccp4', '_event.ccp4', 'event_map.map', '_map.native.ccp4')
        fofc2Strings <- c('_2fofc_cut.cpp4', '_2fofc.cpp4', '^2fofc.map')
        fofcStrings <- c('_fofc_cut.cpp4', '_fofc.ccp4', '^fofc.map')

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
            theFiles <- findFiles(XtalRoot) # row 1 is: 1 = event map, 2 = 2fofc and 3 = fofc
            if(any(is.na(theFiles))){
                maptype <- c('event', '2fofc', 'fofc')
                missfiles <- paste0(maptype[which(is.na(theFiles))], collapse = ' and ')
                if(debug) print(missfiles)
                output$missingFiles <- renderText({
                    sprintf('Unable to find %s maps for this structure. Please check folder for files ending in .map or .cpp4', missfiles)
                }) 
            } else {
                output$missingFiles <- renderText({sprintf('Using: event: %s, 2fofc: %s, fofc: %s', basename(theFiles[1]), basename(theFiles[2]), basename(theFiles[3]))})
            }
            fname <- theFiles[1]
            if(debug) debugMessage(sID=sID, sprintf('Event Map: %s', fname))
            output$progtext <- renderText({'Uploading map files... event map...'}) 
            if(TRUE){   
                event <- readBin(fname, what = 'raw', file.info(fname)$size)
                event <- base64encode(event, size=NA, endian=.Platform$endian)
                session$sendCustomMessage(type="addEvent", 
                    message=list(
                        as.character(event), 
                        as.character(input$isoEvent), 
                        as.character('orange'), 
                        as.character('false'), 
                        as.character(getExt(fname)),
                        as.character(input$boxsize),
                        tolower(as.character(as.logical(input$eventMap)))

                    )
                )
            }

            output$progtext <- renderText({'Uploading map files... 2fofc map...'}) 
            if(TRUE){
                fname <- theFiles[2]
                if(debug) debugMessage(sID=sID, sprintf('2fofc Map: %s', fname))
                event <- readBin(fname, what = 'raw', file.info(fname)$size)
                event <- base64encode(event, size=NA, endian=.Platform$endian)
                session$sendCustomMessage(type="add2fofc", 
                    message=list(
                        event, 
                        as.character(input$iso2fofc), 
                        as.character('blue'), 
                        as.character('false'), 
                        as.character(getExt(fname)),
                        as.character(input$boxsize),
                        tolower(as.character(as.logical(input$twofofcMap)))
                    )
                )
            }

            output$progtext <- renderText({'Uploading map files... fofc map...'}) 
            if(TRUE){
                #fname <- dir(XtalRoot, pattern = '_fofc.ccp4', full.names=T)[1]
                #fname <- dir(XtalRoot, pattern = '^fofc.map', full.names=T)[1]
                fname <- theFiles[3]
                if(debug) debugMessage(sID=sID, sprintf('fofc Map: %s', fname))
                event <- readBin(fname, what = 'raw', file.info(fname)$size)
                event <- base64encode(event, size=NA, endian=.Platform$endian)
                session$sendCustomMessage(type="addfofc_positive", 
                    message=list(
                        event, 
                        as.character(input$isofofc), 
                        as.character('lightgreen'), 
                        as.character('false'), 
                        as.character(getExt(fname)),
                        as.character(input$boxsize),
                        tolower(as.character(as.logical(input$fofcMap)))
                    )
                )
                session$sendCustomMessage(type="addfofc_negative", 
                    message=list(
                        event, 
                        as.character(input$isofofc), 
                        as.character('tomato'), 
                        as.character('true'), 
                        as.character(getExt(fname)),
                        as.character(input$boxsize),
                        tolower(as.character(as.logical(input$fofcMap)))
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
        if(input$fofcMap)   session$sendCustomMessage(type='twiddlefofc_positive',  list(as.character(input$isofofc), as.character(input$boxsize)))
        if(input$fofcMap)   session$sendCustomMessage(type='twiddlefofc_negative',  list(as.character(input$isofofc), as.character(input$boxsize)))
    })

    observeEvent(input$isoEvent,{
        session$sendCustomMessage(type='twiddleEvent',          list(as.character(input$isoEvent),as.character(input$boxsize)))
    })

    observeEvent(input$iso2fofc,{
        session$sendCustomMessage(type='twiddle2fofc',          list(as.character(input$iso2fofc),as.character(input$boxsize)))
    })

    observeEvent(input$isofofc,{
        session$sendCustomMessage(type='twiddlefofc_positive',  list(as.character(input$isofofc),as.character(input$boxsize)))
        session$sendCustomMessage(type='twiddlefofc_negative',  list(as.character(input$isofofc),as.character(input$boxsize)))
    })

    observeEvent(input$assembly2,{
        session$sendCustomMessage(type='updateAssembly', list(as.character(input$assembly2)))
    })

    # Really need to sort this logic ball out...
    observeEvent(input$Xtal, {
            session$sendCustomMessage(type = 'setup', message=list())
            starttime <- Sys.time()
            choice = input$Xtal
            filepath <- dbdat[choice,'Latest.PDB']
            XtalRoot <- try(getRootFP(filepath), silent=T)
            defaultPdbID <- filepath
            defaultShell <- XtalRoot
            output$progtext <- renderText({'Uploading PDB File...'}) 
            tryAddPDB <- try(uploadPDB(filepath=defaultPdbID, input=input), silent=T)
            output$progtext <- renderText({'Uploading PDB File... Done'}) 
            if(inherits(tryAddPDB, 'try-error')){
                defaultPdbID <- ''
                defaultShell <- ''
                session$sendCustomMessage(type="removeAllRepresentations", message=list())
            } else {
                if(!inherits(XtalRoot, 'try-error')){
                    output$progtext <- renderText({'Uploading map files... '}) 
                    tryAddEvent <- try(uploadEMaps(XtalRoot=defaultShell, input=input), silent=T)
                    if(inherits(tryAddEvent, 'try-error')){
                        defaultShell <- ''
                        session$sendCustomMessage(type="removeAllRepresentations", message=list())
                    }
                    endtime <- Sys.time()
                    output$progtext <- renderText({sprintf('Currently Viewing: %s (TimeTaken: %s seconds)', input$Xtal, signif(endtime-starttime, 3))}) 
                }
            }

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

            ligfile <- tail(dir(sprintf('%s/compound', XtalRoot), pattern = '.png', full.names=T),1)
            output$ligimage <- renderImage({
                if(length(ligfile) == 1){
                    list(src = ligfile,
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
            selectizeInput('Xtal', 'Which Structure?', choices = c('', xtalList), selected = query[['xtal']], multiple = FALSE)
        } else {
            selectizeInput('Xtal', 'Which Structure?', choices = c('', xtalList), multiple = FALSE)
        }
    })

    output$isoEventSlider <- renderUI({
        if(input$eventMap){
            sliderInput("isoEvent", "",
                    min = 0, max = 3,
                    value = 1, step = 0.1)
        } else {
            NULL
        }
    })

    output$iso2fofcSlider <- renderUI({
        if(input$twofofcMap){
            sliderInput("iso2fofc", "",
                    min = 0, max = 3,
                    value = 1.5, step = 0.1)
        } else {
            NULL
        }
    })

    output$isofofcSlider <- renderUI({
        if(input$fofcMap){
            sliderInput("isofofc", "",
                min = 0, max = 3,
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

    # Frag View
    output$FragViewnglShiny <- renderNglShiny(
        nglShiny(name = 'nglShiny', list(), width=NULL, height=100)
    )

    uploadPDB2 <- function(filepath){
            syscall <- sprintf('cat %s', filepath)
            if(debug) debugMessage(sID=sID, sprintf('Executing: %s', syscall))
            pdbstrings <- system(syscall, intern = TRUE)
            choice <- paste0(pdbstrings, collapse='\n')
            session$sendCustomMessage(
                type="setapoPDB", 
                message=list(choice))
    }

    getAlignedStagingFolder <- function(){
        stagefold <- getStagingFolder()
        if(length(dir('stagefold', pattern='aligned')) > 0){
            out <- file.path('/dls/science/groups/i04-1/fragprep/staging', input$fragSelect, 'aligned')
        } else {
            out <- stagefold
        }
        return(out)
    }

    getStagingFolder <- function(){
        file.path('/dls/science/groups/i04-1/fragprep/staging', input$fragSelect)
    }

    uploadMol <- function(filepath){
            syscall <- sprintf('cat %s', filepath)
            if(debug) debugMessage(sID=sID, sprintf('Executing: %s', syscall))
            pdbstrings <- system(syscall, intern = TRUE)
            choice <- paste0(pdbstrings, collapse='\n')
            session$sendCustomMessage(
                type="addMol", 
                message=list(choice))
    }

    uploadMol2 <- function(filepath){
            syscall <- sprintf('cat %s', filepath)
            if(debug) debugMessage(sID=sID, sprintf('Executing: %s', syscall))
            pdbstrings <- system(syscall, intern = TRUE)
            choice <- paste0(pdbstrings, collapse='\n')
            session$sendCustomMessage(
                type="addMolandfocus", 
                message=list(choice))
    }
 
    getMolFiles <- function(folderName){
        folderPath <- getAlignedStagingFolder()
        molfiles <- dir(folderPath, rec=T, pattern='.mol', full.names=TRUE)
        names(molfiles) <- basename(molfiles)
        return(molfiles)
    }

    getMetaFiles <- function(folderName){
        folderPath <- getAlignedStagingFolder()
        files <- dir(folderPath, rec=T, pattern='meta.csv', full.names=TRUE)
        return(files)
    }

    showCurrentMetaData <- function(){
        if(debug) debugMessage(sID=sID, sprintf('Fetching Data'))
        files <- getMetaFiles(input$fragSelect)


        try({
            metaoutfile <- file.path(getStagingFolder(), 'metadata.csv')
            if(debug) debugMessage(sID=sID, sprintf('Saving Data to %s', metaoutfile))
            output <- do.call('rbind', lapply(files, function(x) read.csv(x, row.names=NULL, stringsAsFactors=F, header=F)))
            colnames(output) <- c('','crystal_name', 'RealCrystalName', 'smiles','new_smiles','alternate_name','site_name','pdb_entry')
            output2 <- output[!output$site_name=='IGNORE',]
            write.csv(output2[,-1], file=metaoutfile, quote=F)
        }, silent = T)

        try(colnames(output) <- c('', 'Fragalysis Label', 'Crystal Name', 'Smiles', 'New Smiles', 'Alt Fragalysis Label', 'Site Label', 'PDB Entry'), silent=T)
        rownames(output) <- output[,2]
        return(output[,-1])
    }

    updateTable <- function(){
        #gsearch <- input$therow_search
        #gsearch <- ifelse(is.null(gsearch),'',gsearch)
        #rowsel <- input$therow_rows_selected
        #rowsel <- ifelse((is.null(rowsel) | rowsel==''), 'single', list(mode='single', selected=rowsel, target='row'))
        #csearch <- input$therow_search_columns
        #csearch <- ifelse(is.null(csearch), list(), csearch)
        #print(gsearch)
        #print(rowsel)
        #print(csearch)
        if(debug) debugMessage(sID=sID, sprintf('Updating Table'))
        dummy <- rbind(c(1,1), c(1,1))
        mdr <- reactiveValues(x=showCurrentMetaData())
        output$therow <-  DT::renderDataTable({datatable(dummy, selection = 'single', options = list(pageLength = 200))})
        output$therow <-  DT::renderDataTable({
            datatable(mdr$x, 
                filter='top', 
                selection = 'single',#rowsel, 
                options = list(
                    #searchCols = csearch,
                    scrollX=TRUE,
                    #search = list(
                    #    regex = FALSE, 
                    #    caseInsensitive = TRUE, 
                    #    search = gsearch),
                    pageLength = 200
                )
            )
        })
        updateSelectizeInput(session, 'site_name', choices = sort(mdr$x[,6]))
        updateSelectizeInput(session, 'site_name2', choices = sort(mdr$x[,6]))
        if(debug) debugMessage(sID=sID, sprintf('Updated Table'))
    }

    output$writeButton <- renderUI({
        if(input$desync) {
            actionButton('write', 'Write metadata to table', style="background-color: #FF0000")
        } else {
            actionButton('write', 'Write metadata to table')
        }

    })

    observeEvent(input$fragSelect,{
        if(debug) debugMessage(sID=sID, sprintf('Selecting: %s', input$fragSelect))
        folderPath <- getAlignedStagingFolder()
        apofile <- tail(dir(folderPath, rec =T, pattern = 'apo.pdb', full.names=TRUE),1)
        molfiles <- getMolFiles(input$fragSelect)
        molfil <- names(molfiles)
        updateSelectInput(session, 'goto', choices = molfil)
        updateTable()     
        tryAddPDB <- try(uploadPDB2(filepath=apofile), silent=T)
        molout <- try(sapply(molfiles, uploadMol), silent=T)
        
    })

    observeEvent(input$gonext, {
        molfiles <- getMolFiles(input$fragSelect)
	    molbase <- names(molfiles)
        nmol <- length(molfiles)
        id <- which(molbase == input$goto)
        next_id <- id + 1
        if(next_id > nmol) next_id <- 1 # Overflow back to start of list
        # Cycle along to next ligand in molfil
        if(debug) debugMessage(sID=sID, sprintf('Switching to: %s', molbase[next_id]))
        updateSelectInput(session, 'goto', selected = molbase[next_id], choices=molbase)
    })

    observeEvent(input$goback, {
        molfiles <- getMolFiles(input$fragSelect)
	    molbase <- names(molfiles)
        nmol <- length(molfiles)
        id <- which(molbase == input$goto)
        next_id <- id - 1
        if(next_id < 1) next_id <- nmol # Underflow to end of list
        if(debug) debugMessage(sID=sID, sprintf('Switching to: %s', molbase[next_id]))
        updateSelectInput(session, 'goto', selected = molbase[next_id], choices=molbase)
    })

    observeEvent(input$goto, {
        output$metastatus <- renderText({'STATUS: Pending...'}) 
        molfiles <- getMolFiles(input$fragSelect)
        folder <- dirname(molfiles[input$goto])
        # Fill as it is seen:
        files <- dir(folder, pattern='.csv', full.names=T)
        fn <- strsplit(input$goto, split='[.]')[[1]][1]
        updateTextInput(session, 'newCrystalName', value = fn)
        if(length(files) > 0){
            dat <- read.csv(files, row.names=NULL, stringsAsFactors=F, header=F)
            updateTextInput(session, 'crysname', value = dat[1, 3])
            updateTextInput(session, 'smiles', value = dat[1, 4])
            updateTextInput(session, 'new_smiles', value = dat[1, 5])
            updateTextInput(session, 'alternate_name', value = dat[1, 6])
            updateSelectizeInput(session, 'site_name', select = dat[1, 7])
            updateTextInput(session, 'pdb_entry', value = dat[1, 8])
        } else {
            inputfolder <- file.path('/dls/science/groups/i04-1/fragprep/staging', input$fragSelect)
            smilesfn <- strsplit(input$goto, split='[.]')[[1]][1]
            smilestr <- system(sprintf('cat %s/%s/%s_smiles.txt', inputfolder,smilesfn,smilesfn), intern=T)
            updateTextInput(session, 'crysname', value = '')
            updateTextInput(session, 'smiles', value = smilestr)
            updateTextInput(session, 'new_smiles', value = '')
            updateTextInput(session, 'alternate_name', value = '')
            updateSelectizeInput(session, 'site_name', selected = '')
            updateTextInput(session, 'pdb_entry', value = '')
        }

        # Go to specific ligand do not edit go next loop
        if(debug) debugMessage(sID=sID, sprintf('Selected: %s', input$goto))
        if(debug) debugMessage(sID=sID, sprintf('trying to view: %s', molfiles[input$goto]))
        gogogo <- try(uploadMol2(molfiles[input$goto]), silent=T)
    })

    observeEvent(input$write, {
        if(debug) debugMessage(sID=sID, sprintf('Writing a row'))
        output$metastatus <- renderText({'STATUS: Written!'}) 
        molfiles <- getMolFiles(input$fragSelect)
        folder <- dirname(molfiles[input$goto])
        if(debug) debugMessage(sID=sID, sprintf('Writing %s to %s', input$site_name, input$goto))
        fragname <- gsub('.mol', '', input$goto)
        fn <- file.path(folder, sprintf('%s_meta.csv', fragname))
        output <- t(c(fragname, 
                    input$crysname,
                    input$smiles, 
                    input$new_smiles, 
                    input$alternate_name, 
                    input$site_name, 
                    input$pdb_entry
                    ))
        write.table(output, file = fn, quote = F, col.names=F, sep=',')
        if(!input$desync) updateTable()
    })

    # On Table Rowclick # Potentially slow? Unneeded? # Go back to
    observeEvent(input$therow_rows_selected, {
        if(debug) debugMessage(sID=sID, sprintf('Selecting Row'))
        molfiles <- getMolFiles(input$fragSelect)
        molbase <- names(molfiles)
        mdr <- reactiveValues(x=showCurrentMetaData())
        choice <- mdr$x[input$therow_rows_selected,1]
        updateSelectizeInput(session, 'goto', selected = sprintf('%s.mol', choice), choices=molbase)
    })

    output$massChange <- renderText({'Update Site Labels'})

    observeEvent(input$mcl, {
        if(debug) debugMessage(sID=sID, sprintf('Changing Labels'))
        oldlabel <- input$site_name2
        newlabel <- input$new_label
        # Perhaps open a modal
        folderPath <- getAlignedStagingFolder()
        metacsv <- dir(folderPath, pattern='meta.csv', rec=T, full.names=T)
        for(afile in metacsv){
            message(afile)
            afile_data <- read.csv(afile, row.names=NULL, stringsAsFactors=F, header=F)
            currentlabel <- afile_data[1,7]
            if(!is.na(currentlabel)){
                if(afile_data[1,7] == oldlabel){
                    afile_data[1,7] <- newlabel
                    write.table(afile_data, file=afile, quote=F, col.names=F, row.names=F, sep=',')
                }
            }
        }
        if(!input$desync) updateTable()
    })

    observeEvent(input$restartViewer, {
        try({session$sendCustomMessage(type="removeAllComponents", message=list())}, silent=T)
    })

    observeEvent(input$changeName, {
        if(debug) debugMessage(sID=sID, sprintf('Changing Names'))
        oldname <- strsplit(input$goto, split='[.]')[[1]][1]
        newname <- input$newCrystalName
        molfiles <- getMolFiles(input$fragSelect)
        folder <- dirname(molfiles[input$goto])
        newfolder <- gsub(oldname, newname, folder)

        # sanity check newname for / . etc Reject with modal.

        # Have a modal with accept...
        oldfiles <- dir(folder, rec=T, full.names=T)
        newfiles <- gsub(oldname, newname, oldfiles)
        # Change folder name

        # check if name exists
        existingNames <- names(molfiles)
        dupe <- newname %in% existingNames
        if(dupe){
            showModal(modalDialog(title = "A fragment with this name already exists",
                    sprintf('The name %s already exists, please use a different name', newname) 
                    , easyClose=TRUE))
        } else {
            # If newname has none of these characters, legal = TRUE
            illegal <- '[./~\\!\"\'£$%^&*()@<>?|+`¬]'
            legal <- !grepl(illegal, newname)
            if(legal){
                file.rename(folder, newfolder)
                # Reget Old Files
                oldfiles <- dir(newfolder, rec=T, full.names=T)
                file.rename(oldfiles, newfiles)
                # Edit label in metadata...
                metacsv <- dir(newfolder, pattern='_meta.csv', rec=T, full.names=T)
                if(length(metacsv) > 0){
                    metadat <- read.csv(metacsv, row.names=NULL, stringsAsFactors=F, header=F)
                    metadat[1,2] <- newname
                    write.table(metadat, file=metacsv, quote=F, col.names=F, row.names=F, sep=',')
                } else {
                showModal(modalDialog(title = "You used some naughty characters",
                    sprintf('Please use a file name that does not include any of these symbols: 
                        %s', illegal) 
                    , easyClose=TRUE))
                }
            }
        }
        # Mandatory Update
        if(!input$desync) updateTable()
        molfiles <- getMolFiles(input$fragSelect)
        molbase <- names(molfiles)
        updateSelectizeInput(session, 'goto', selected = sprintf('%s.mol', newname), choices=molbase)
    })

    observeEvent(input$updateTable,{
        updateTable()
    })

    # Frag Chat
    output$scrollDialog <- renderText({'Scroll for More...'})

    observeEvent(input$updateSlackChannels,{
        channels <- getChannelList()
        channelSelect <- channels[,2]
        names(channelSelect) <- channels[,1]
        updateSelectizeInput(session, "channelSelect", select = input$channelSelect, choices = names(channelSelect))
        refreshChat(channel = channelSelect[input$channelSelect]) 
    })

    observeEvent(input$channelSelect, {
        message(input$channelSelect)
        channels <- getChannelList()
        channelSelect <- channels[,2]
        names(channelSelect) <- channels[,1]
        refreshChat(channel = channelSelect[input$channelSelect]) 
        output$chatURL <- renderText({sprintf('https://xchemreview.slack.com/archives/%s', channelSelect[input$channelSelect])})
    })

    observeEvent(input$slackSubmit, {
        if(input$slackUser == '' | input$TextInput == ''){
            showModal(modalDialog(title = "Cannot send empty messages.",
                    'Please add some text to the Name and Text Fields.', easyClose=TRUE))
        } else {
        sendMessageToSlack(channel = channelSelect[input$channelSelect], 
                            message = sprintf('on behalf of: %s. \n %s', input$slackUser, input$TextInput), 
                            name=input$slackUser)
        }
        refreshChat(channel = channelSelect[input$channelSelect]) 
    })
} # Server
