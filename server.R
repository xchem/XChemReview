server <- function(input, output, session){
    fedid <- Sys.getenv(x = 'SHINYPROXY_USERNAME', unset = "", names = NA)
    sID <- fedid
    debug = TRUE
    if(debug) debugMessage(sID=sID, sprintf('Session init'))
    session$allowReconnect(FALSE)
    sessionDisconnect <- function(){
        # Clean-up:
        clnup <- file.path(configuration$xcr_path, 'www', sprintf('report_%s.pdf', sID))
        try(file.remove(clnup), silent=T)
        debugMessage(sID=sID, 'Disconnected')
    }
    session$onSessionEnded(sessionDisconnect)
    query <- parseQueryString(isolate(session$clientData$url_search))
    if(!is.null(query[['key']])){
        debugMessage(sID=sID, sprintf('Decoding Key: %s', query[['key']]))
        target_list2 <- sort(c(decrypt(query[['key']])))
        target_list <- target_list2
        fragfolders <- c('', target_list2)
    }
    sessionTime <- reactive({epochTime()})
    possDec <- c("", "Release", "More Refinement", "More Experiments", "Reject")
    possAns <- possAns2 <- c('Select Decision')
    possRes <- list()
    possRes[['Release']] <- c('High Confidence', 'Clear Density, Unexpected Ligand', 'Correct Ligand, Weak Density', 'Low Confidence', 'No Ligand Present')
    possRes[["More Refinement"]] <- c('Check ligand conformation',
        'Check sidechain rotamers',
        'Check Rfactors',
        'Check that refinement converged',
        'Improve water model',
        'Incomplete Event Density',
        'Build alternate conformations',
        'Fix geometry',
        'Trim down ligand',
        'Density did not load',
        'Other')
    possRes[["More Experiments"]] <- c('Get better density',
        'Get better resolution',
        'Confirm ligand identity',
        'Check if ligand changed',
        'Other')
    possRes[["Reject"]] <- c(
        'Density not convincing',
        'Too few interactions',
        'Binding site too noisy',
        'Not the ligand',
        'Duplicate Pose',
        'Other')
    possDec_int <- 1:4
    names(possDec_int) <- c("Release", "More Refinement", "More Experiments", "Reject")

    sessionlist <- reactiveValues()
    sessionlist$current_emaps <- ''
    sessionlist$lig_name <- ''
    sessionlist$apo_file <- ''
    sessionlist$mol_file <- ''
    sessionlist$event <- ''
    sessionlist$twofofc_file <- ''
    sessionlist$fofc_file <- ''
    sessionlist$target_name <- ''
    sessionlist$resolution <- ''
    sessionlist$r_free <- ''
    sessionlist$r_cryst <- ''
    sessionlist$ramachandran_outliers <- ''
    sessionlist$rmsd_angles <- ''
    sessionlist$rmsd_bonds <- ''
    sessionlist$pdb_latest <- ''
    sessionlist$lig_id <- ''
    sessionlist$xtalroot <- ''
    sessionlist$rowname <- ''
    sessionlist$fumeta <- ''
    sessionlist$fv_warn <- ''
    sessionlist$busterReport <- ''
    sessionlist$isotype <- 'value'

    output$fv_warn <- renderPrint({sessionlist$fv_warn})

    observeEvent(input$aq_protein, {
        debugMessage(sID=sID, sprintf('(input$%s) Choosing Target: %s', 'aq_protein', input$aq_protein))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_protein', 'calling reactiviseData()'))
        r1 <- reactiviseData(inputData=inputData, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_protein', 'calling updateAQPTable()'))
        output$aqp <- updateAQPTable(r=r1, pl=100, input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_protein', 'calling prebuffer_proxy(aqpproxy)'))
        prebuffer_proxy(proxy=aqpproxy, input=input)
    })
    
    observeEvent(input$protein, {
        debugMessage(sID=sID, sprintf('(input$%s) Choosing Target %s', 'protein', input$protein))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'protein', 'calling reactiviseData()'))
        r1 <- reactiviseData(inputData=inputData, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'protein', 'calling updateMainTable()'))
        output$reviewtable <- updateMainTable(r1=r1, pl=100, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'protein', 'calling prebuffer_proxy(reviewtableproxy)'))
        prebuffer_proxy(proxy=reviewtableproxy, input=input)
    })

    observeEvent(input$out4, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out4', input$out4))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out4', 'calling reactiviseData()'))
        r1 <- reactiviseData(inputData=inputData, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out4', 'calling updateMainTable()'))
        output$reviewtable <- updateMainTable(r1=r1, pl=100, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out4', 'calling prebuffer_proxy(reviewtableproxy)'))
        prebuffer_proxy(proxy=reviewtableproxy, input=input)
    })

    observeEvent(input$out5, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out5', input$out5))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out5', 'calling reactiviseData()'))
        r1 <- reactiviseData(inputData=inputData, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out5', 'calling updateMainTable()'))
        output$reviewtable <- updateMainTable(r1=r1, pl=100, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out5', 'calling prebuffer_proxy(reviewtableproxy)'))
        prebuffer_proxy(proxy=reviewtableproxy, input=input)
    })

    observeEvent(input$out6, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out6', input$out6))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out6', 'calling reactiviseData()'))
        r1 <- reactiviseData(inputData=inputData, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out6', 'calling updateMainTable()'))
        output$reviewtable <- updateMainTable(r1=r1, pl=100, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'out6', 'calling prebuffer_proxy(reviewtableproxy)'))
        prebuffer_proxy(proxy=reviewtableproxy, input=input)
    })

    # Selector Stuff:
    debugMessage(sID=sID, sprintf('Fetch Review Table Data'))
    review_data <- trygetReviewData(configuration=configuration, target_list=target_list)
    mtzzz <- review_data[,c('crystal_name', 'mtz_latest')]
    updateSelectInput(session, 'protein', selected = '', choices=c('', sort(unique(as.character(review_data$target_name)))))
    updateSelectInput(session, 'fpe_target', selected = '', choices=c('', sort(unique(as.character(review_data$target_name)))))
    updateSelectInput(session, 'config_target', selected = '', choices=c('', sort(unique(as.character(review_data$target_name)))))
    inputData <- restartSessionKeepOptions(configuration=configuration, target_list=target_list)
    r1 <- reactiviseData(inputData=inputData, input=input)
    reviewtableproxy <- DT::dataTableProxy('reviewtable')
    aqpproxy <- DT::dataTableProxy('aqp')
    flexplotData <- flexPlotDataFun(r1=r1, input=input)
    output$flexplot1 <- updateFlexPlot(flexdata=flexplotData)

    interestingData <- reactive({
        di = inputData()$decision_str
        di <- sapply(di, function(x) ifelse(is.null(x), NA, x))

        if(!is.null(input$protein_to_summarize)){
            if(!input$protein_to_summarize == ''){
                idstocheck <- as.character(inputData()$target_name)==input$protein_to_summarize
                di = di[idstocheck]
            }
        }

        output = list(
            'total_ligands' = length(di),
            'total_reviewed' = sum(!is.na(di)),
            'to_review' = sum(is.na(di))
        )

        return(output)
    })
    debugMessage(sID=sID, sprintf('Finalised...'))
    debugMessage(sID=sID, sprintf('Fetch Fragview Data'))
    fvd <- getFragalysisViewData(configuration=configuration, target_list=target_list2)
    fragview_data <- reactivegetFragalysisViewData(configuration=configuration, target_list=target_list2)
    updateSelectInput(session, 'fragSelect', selected='', choices=fragfolders)

    fragview_input <- react_fv_data(data=fragview_data, input=input)
    fragview_table_data <- react_fv_data2(data=fragview_data, input=input)

    fragviewproxy <- DT::dataTableProxy('therow')
    output$as_message <- renderUI({HTML('Select Atom: Alt + Click <br/> Select Side Chain: Alt + Ctrl + Click <br/> Select Whole Residue: Shift + Ctrl + Alt + Click')})
    debugMessage(sID=sID, sprintf('Finalised...'))



    observeEvent(input$as_clear, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'as_clear', 'Resetting Atoms'))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'as_clear', 'Calling js as_resetclicked'))
        session$sendCustomMessage(type = 'as_resetclicked', list())
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'as_clear', 'Re-init atomstoquery table'))
        atomstoquery$data <- data.frame(name = character(),
            index = character(),
            comment = character(),
            stringsAsFactors = FALSE)

    })

    fv_values <- reactiveValues()
    fv_values$apofiles <- c()
    fv_values$molfiles <- c()
    fv_values$molfil <- c()

    observeEvent(input$fragSelect,{
        debugMessage(sID=sID, sprintf('(input$%s) Choosing fragview target: %s', 'fragSelect', input$fragSelect))
        #folderPath <- getAlignedStagingFolder()
        if(input$fragSelect == ''){
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'Nothing Selected. Render nothing...'))
            output$therow <- DT::renderDataTable({DT::datatable(data.frame('a' = '', 'Nothing Selected' = 'Please Select a Target', stringsAsFactors=FALSE))}, server=FALSE)
        } else {
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'Fetch apo Files'))
            fv_values$apofiles <- as.character(isolate(fragview_table_data()$apo_pdb))
            apo_existing <- sapply(fv_values$apofiles, file.exists)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'Fetch mol Files'))
            fv_values$molfiles <- as.character(isolate(fragview_table_data()$lig_mol_file))
            fv_values$apofiles <- fv_values$apofiles[apo_existing]
            fv_values$molfil <- gsub('.mol', '', basename(fv_values$molfiles))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'Update input$goto'))
            updateSelectInput(session, 'goto', choices = fv_values$molfil)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'calling react_fv_data()'))
            fragview_input <- react_fv_data(data=fragview_data, input=input) # Filter missing files here??
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'calling updateMainTable2()'))
            output$therow <- updateMainTable2(r1=fragview_input, pl=100, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'Calling prebuffer_proxy(fragviewproxy)'))
            prebuffer_proxy(proxy=fragviewproxy, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'Calling uploadApoPDB()'))
            tryAddPDB <- try(uploadApoPDB(filepath=fv_values$apofiles[1], repr='cartoon', focus=TRUE, session=session), silent=T)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'fragSelect', 'Calling uploadUnfocuessedMol(), multiple times...'))
            molout <- try(sapply(fv_values$molfiles, uploadUnfocussedMol, session=session), silent=T)
        }   
    })

    observeEvent(input$bfactor, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'bfactor', 'Toggle Bfactors in Review Tab'))
        if(input$bfactor){
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'bfactor', 'Toggle ON'))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'bfactor', 'Clear existing bfactors'))
            uploadBFactors(filepath=sessionlist$apo_file, clear=TRUE, session=session)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'bfactor', 'Turn off original mol visibility'))
            updateVisability(name='mol', bool=FALSE, session=session)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'bfactor', 'calling uploadBFactors()'))
            uploadBFactors(filepath=gsub('.mol', '.pdb', sessionlist$mol_file), clear=FALSE, session=session)
        } else {
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'bfactor', 'Toggle OFF'))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'bfactor', 'Clear existing bfactors'))
            clearWindowField(id='bfactors', session=session)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'bfactor', 'Turn on original mol visibility'))
            updateVisability(name='mol', bool=TRUE, session=session)  
        }
    })

    observeEvent(input$aq_bfactor,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_bfactor', 'Toggle Bfactors in Atom Quality zone'))
        updateTabsetPanel(session, 'ctab', selected = 'Controls')
        if(input$aq_bfactor){
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_bfactor', 'Toggle ON'))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_bfactor', 'Clear existing bfactors'))
            uploadBFactors(filepath=sessionlist$apo_file, clear=TRUE, session=session)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_bfactor', 'Turn off original mol visibility'))
            updateVisability(name='mol', bool=FALSE, session=session)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_bfactor', 'calling uploadBFactors()'))
            uploadBFactors(filepath=gsub('.mol', '.pdb', sessionlist$mol_file), clear=FALSE, session=session)
        } else {
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_bfactor', 'Toggle OFF'))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_bfactor', 'Clear existing bfactors'))
            clearWindowField(id='bfactors', session=session)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_bfactor', 'Turn off original on visibility'))
            updateVisability(name='mol', bool=TRUE, session=session) 
        }
    })

    observeEvent(input$gonext, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'gonext', ''))
        molfiles <- fv_values$molfiles
        molbase <- fv_values$molfil
        nmol <- length(molfiles)
        id <- which(molbase == input$goto)
        next_id <- id + 1
        if(next_id > nmol) next_id <- 1 # Overflow back to start of list
        # Cycle along to next ligand in molfil
        debugMessage(sID=sID, sprintf('Switching to: %s', molbase[next_id]))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'gonext', 'Update input$goto'))
        updateSelectInput(session, 'goto', selected = molbase[next_id], choices=molbase)
    })

    observeEvent(input$goback, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'goback', ''))
        molfiles <- fv_values$molfiles
        molbase <- fv_values$molfil
        nmol <- length(molfiles)
        id <- which(molbase == input$goto)
        next_id <- id - 1
        if(next_id < 1) next_id <- nmol # Underflow to end of list
        debugMessage(sID=sID, sprintf('Switching to: %s', molbase[next_id]))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'goback', 'Update input$goto'))
        updateSelectInput(session, 'goto', selected = molbase[next_id], choices=molbase)
    })

    observeEvent(input$goto, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'goto', ''))
        output$metastatus <- renderText({'STATUS: Pending...'})
        molfiles <- fv_values$molfiles
        molbase <- fv_values$molfil
        names(molfiles) <- molbase
        if(!input$goto == ''){
            folder <- dirname(molfiles[input$goto])
            mol_file <- molfiles[input$goto]
            smi_file <- gsub('.mol', '_smiles.txt', mol_file)
            smilestr <- system(sprintf('cat %s', smi_file), intern=T)

            choices <- unique(c('', as.character(isolate(fragview_table_data()[, 'site_Label']))))
            # Fill Form as seen
            updateTextInput(session, 'crysname', value = input$goto)
            updateTextInput(session, 'smiles', value = smilestr)
            tofilldata <- inputData()[inputData()$ligand_name == input$goto,]

            alt_name = as.character(isolate(fragview_table_data()[input$goto, 'alternate_name']))
            if(blankNAorNull(x=alt_name)) alt_name = tools::file_path_sans_ext(basename(tofilldata$cif))
            updateTextInput(session, 'alternate_name', value = alt_name)
            updateTextInput(session, 'new_smiles', value = as.character(isolate(fragview_table_data()[input$goto, 'new_smiles'])))

            site_name = as.character(isolate(fragview_table_data()[input$goto, 'site_Label']))
            if(blankNAorNull(x=site_name)){
                site_name = try(as.character(read.csv(gsub('.mol', '_meta.csv', mol_file), header=F, stringsAsFactors=F))[7], silent=T)
                if(inherits(site_name, 'try-error')) site_name = ''
            }
            updateSelectizeInput(session, 'site_name', selected = site_name, choices=c(site_name,choices))

            updateTextInput(session, 'pdb_entry', value = as.character(isolate(fragview_table_data()[input$goto, 'pdb_id'])))
            # Go to specific ligand do not edit go next loop
            debugMessage(sID=sID, sprintf('Selected: %s', input$goto))
            debugMessage(sID=sID, sprintf('trying to view: %s', molfiles[input$goto]))
            gogogo <- try(uploadMolAndFocus2(filepath=mol_file, session=session), silent=T)
            if(!file.exists(mol_file)){
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'goto', 'File is missing, check pipelines or soakdb?'))
                sessionlist$fv_warn <- 'WARNING: .mol File is MISSING SET TO IGNORE OR INVESTIGATE' 
	        } else {
	            sessionlist$fv_warn <- '.mol File found!'
	        }
        }
    })

    output$writeButton <- renderUI({
        if(is.null(input$desync)){
            actionButton('write', 'Write metadata to table')
        } else if(input$desync) {
            actionButton('write', 'Write metadata to table', style="background-color: #FF0000")
        } else {
            actionButton('write', 'Write metadata to table')
        }
    })

    observeEvent(input$updateTable,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'updateTable', ''))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'updateTable', 'calling reactivegetFragalysisViewData'))
        fragview_data <- reactivegetFragalysisViewData(configuration=configuration, target_list=target_list2)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'updateTable', 'calling react_fv_data'))
        fragview_input <- react_fv_data(data=fragview_data, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'updateTable', 'calling react_fv_data2 . Nice one tyler...'))
        fragview_table_data <- react_fv_data2(data=fragview_data, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'updateTable', 'calling updateMainTable2'))
        output$therow <- updateMainTable2(r1=fragview_input, pl=100, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'updateTable', 'calling prebuffer_proxy'))
        prebuffer_proxy(proxy=fragviewproxy, input=input)
    })

    observeEvent(input$write, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'write', 'Writing a Row'))
        rn <- rownames(isolate(fragview_table_data()))
        ids <- isolate(fragview_table_data()$id)
        names(ids) <- rn
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'write', 'calling updateOrCreateRow'))
        updateOrCreateRow(ligand_name_id=as.character(ids[input$goto]),
                          fragalysis_name=as.character(input$goto),
                          original_name=as.character(rsplit(input$goto, '_', 1)[1]),
                          site_label=as.character(input$site_name),
                          new_smiles=as.character(input$new_smiles),
                          alternate_name=as.character(input$alternate_name),
                          pdb_id=as.character(input$pdb_entry),
                          configuration=configuration)
        output$metastatus <- renderText({'STATUS: Written!'})
        if(!input$desync){
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'write', 'calling reactivegetFragalysisViewData'))
            fragview_data <- reactivegetFragalysisViewData(configuration=configuration, target_list=target_list2)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'write', 'react_fv_data'))
            fragview_input <- react_fv_data(data=fragview_data, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'write', 'react_fv_data2'))
            fragview_table_data <- react_fv_data2(data=fragview_data, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'write', 'updateMainTabe2'))
            output$therow <- updateMainTable2(r1=fragview_input, pl=100, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'write', 'prebuffer_proxy'))
            prebuffer_proxy(proxy=fragviewproxy, input=input)
        }
    })

    # On Table Rowclick # Potentially slow? Unneeded? # Go back to
    observeEvent(input$therow_rows_selected, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'therow_rows_selected', 'Selecting Row!'))
        molfiles <- fv_values$molfiles
        molbase <- fv_values$molfil
        choice = isolate(rownames(fragview_input())[input$therow_rows_selected])
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'therow_rows_selected', 'Updating input$goto'))
        updateSelectizeInput(session, 'goto', selected = choice, choices=molbase)
    })

    # Not sure I remember what this does to be honest. I think it's to maintain some form of sanity between tabs...
    observeEvent(input$protein,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'protein', 'Selecting protein'))
        updateSelectInput(session, 'fpe_target', selected=input$protein)
    })
    observeEvent(input$fpe_target,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'fpe_target', 'Selecting fragview target'))
        updateSelectInput(session, 'protein', selected=input$fpe_target)
    })

    output$progressBox1 <- renderInfoBox({
        infoBox('Total Reviewed', isolate(interestingData()$total_reviewed), icon = icon('thumbs-up', lib = 'glyphicon'), color = 'red')
    })
    output$approvalBox1 <- renderInfoBox({
        infoBox('Number of Ligands', isolate(interestingData()$total_ligands), icon = icon('thumbs-up', lib = 'glyphicon'), color = 'red')
    })

    # Form Handler
    fieldsAll <- c("name", 'ligand', "decision", "reason", "comments")
    formData <- reactive({
        data <- sapply(fieldsAll, function(x) paste0(input[[x]], collapse='; '))
        # Get Crystal ID
        # This bits are wrong?
        # data[2] should be: sessionlist$rowname
        xtalname <- data[2]
        data[2] <- sessionlist$rowname
        data <- c(data[1], possDec_int[data[3]] ,data[3:4], timestamp = epochTime(), review_data[data[2],'ligand_id'], review_data[data[2],'crystal_id'], data[5])
        list(data=data.frame(
            fedid=data[1],
            decision_int=as.numeric(data[2]),
            decision_str=data[3],
            reason=data[4],
            time_submitted=data[5],
            Ligand_name_id=as.integer(data[6]),
            crystal_id=as.integer(data[7]),
            comment=as.character(data[8]),
            stringsAsFactors=FALSE
        ), xtalname=xtalname)
    })

    observeEvent(input$ok, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ok', ''))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ok', 'calling restartSessionKeepOptions'))
        inputData <- restartSessionKeepOptions(configuration=configuration, target_list=target_list)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ok', 'calling reactiviseData'))
        r1 <- reactiviseData(inputData=inputData, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ok', 'calling updateMainTable'))
        output$reviewtable <- updateMainTable(r1=r1, pl=100, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ok', 'calling prebuffer_proxy'))
        prebuffer_proxy(proxy=reviewtableproxy, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ok', 'calling flexPlotDataFun'))
        flexplotData <- flexPlotDataFun(r1=r1, input=input)
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ok', 'calling updateFlexPlot'))
        output$flexplot1 <- updateFlexPlot(flexdata=flexplotData)
        sessionTime <- reactive({epochTime()})
        removeModal()
    })

    # Upon Main Page Submit
    observeEvent(input$submit, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'Submitting Review'))
        fData <- formData()[[1]]
        xtaln <- formData()[[2]]
        if(any(fData[1:7] %in% c('', ' '))) {
            showModal(modalDialog(title = "Please fill all fields in the form",
                "One or more fields have been left empty. Please provide your FedID, a decision and reason(s) before clicking submit.",
                easyClose=TRUE, footer = tagList(modalButton("Cancel"))
            ))
        } else {
             # Get ID and check ID ...
            xId <- fData[ ,'Ligand_name_id']
            if(sessionGreaterThanMostRecentResponse(id=xId, sessionTime=sessionTime(), configuration=configuration)){
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'Attempting to save Response (calling saveData)'))
                saveData(fData, xtaln, email_list=emailListperStructure, configuration=configuration)
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'Response Saved'))
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'calling resetForm'))
                inputData <- resetForm(configuration=configuration, session=session, r1=r1, possDec=possDec, target_list=target_list)
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'calling reactiviseData'))
                r1 <- reactiviseData(inputData=inputData, input=input)
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'calling updateMainTable'))
                output$reviewtable <- updateMainTable(r1=r1, pl=100, input=input)
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'calling prebuffer_proxy'))
                prebuffer_proxy(proxy=reviewtableproxy, input=input)
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'calling flexPlotDataFun'))
                flexplotData <- flexPlotDataFun(r1=r1, input=input)
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'calling updateFlexPlot'))
                output$flexplot1 <- updateFlexPlot(flexdata=flexplotData)
                sessionTime <- reactive({epochTime()})
            } else {
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit', 'calling displayModalWhoUpdated'))
                displayModalWhoUpdated(id=xId, configuration=configuration)
            }
        }
    })


    observeEvent(input$submit_atoms, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit_atoms', 'Submitting bad atoms'))
        if(any(as.character(atomstoquery$data$comment) %in% c('', ' '))){
            showModal(modalDialog(title = 'You have flagged some atoms',
                'Please annotate the selected atoms in the Atom Selection tab by double clicking on the comment cells. If you accidentally flagged an atom, try reloading the structure and resubmitting your review!',
                easyClose=TRUE))
        } else {
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit_atoms', 'Attempting to Write Atoms'))
            writeAtoms(ligand_id=as.character(sessionlist$ligand_id), configuration=configuration, atomstoquery=atomstoquery, sessionlist=sessionlist)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit_atoms', 'Success!')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit_atoms', 'calling restartSessionKeepOptions')))
            inputData <- restartSessionKeepOptions(configuration=configuration, target_list=target_list)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit_atoms', 'calling reactiviseData'))
            r1 <- reactiviseData(inputData=inputData, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit_atoms', 'calling updateAQPTable'))
            output$aqp <- updateAQPTable(r=r1, pl=100, input)
            aqpproxy <- DT::dataTableProxy('aqp')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'submit_atoms', 'calling prebuffer_proxy'))
            prebuffer_proxy(proxy=aqpproxy, input=input)
        }
    })

    atomstoquery <- reactiveValues()
    atomstoquery$data <- data.frame(name=character(),
                 index=character(),
                 comment=character(),
                 stringsAsFactors=FALSE)

    output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))}, server=FALSE)

    observeEvent(input$clickedAtoms, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'clickedAtoms', ''))
        newdat <- isolate(atomstoquery$data)
        new <- which(!as.character(input$clickNames) %in% as.character(newdat$name))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'clickedAtoms', new))
        for(i in new){
            newdat <- rbind(newdat, data.frame(name = input$clickNames[i], index = input$clickedAtoms[i], comment = '', stringsAsFactors=FALSE))
        }
        tokeep <- as.character(newdat$name) %in% as.character(input$clickNames)
        newdat <- newdat[tokeep,]
        atomstoquery$data <- newdat
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'clickedAtoms', 'Updating Table'))
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))), options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))}, server=FALSE)
    })

    observeEvent(input$atoms_cell_edit, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'atoms_cell_edit', 'Editting a Cell'))
        info = input$atoms_cell_edit
        str(info)
        i = info$row
        j = info$col
        v = info$value
        update <- isolate(atomstoquery$data)
        update[i, j] <- as.character(v)
        atomstoquery$data <- update
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))), options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))}, server=FALSE)
    })


    observeEvent(input$write_selected, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'write_selected', 'Writing to selected rows'))
        idx <- isolate(input$atoms_rows_selected)
        update <- isolate(atomstoquery$data)
        update[idx, 3] <- as.character(input$atom_text)
        atomstoquery$data <- update
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))), options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))}, server=FALSE)
    })

    observeEvent(input$write_all,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'write_all', 'Writing to all rows'))
        update <- isolate(atomstoquery$data)
        idx <- 1:nrow(update)
        update[idx, 3] <- as.character(input$atom_text)
        atomstoquery$data <- update
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))), options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))}, server=FALSE)
    })


    messageData <- data.frame(
        from = c('Frank', 'Frank'),
        message  = c('Hello', 'Frank Here')
    )

    output$messageMenu <- renderMenu({
        msgs <- apply(messageData, 1, function(row) messageItem(from = row[['from']], message = row[['message']]))
        dropdownMenu(type = 'messages', .list = msgs)
    })

    output$notifications <- renderMenu({
        dropdownMenu(
            type = 'notifications',
            notificationItem(text = '5 Users Todays', icon('users')),
            notificationItem(text = '12 Items Delivered', icon('truck'), status='success'),
            notificationItem(text = 'Server load at 86%', icon = icon('exclamation-triangle'), status = 'warning')
        )
    })

    output$tasks <- renderMenu({
        dropdownMenu(
            type = 'tasks', badgeStatus = 'success',
            taskItem(value = 10, color = 'red', '70X'),
            taskItem(value = 95, color = 'green', 'Mpro'),
            taskItem(value = 40, color = 'yellow', 'All')
        )
    })

    # Dynamically Rendered UI elements!
    output$flex1 <- renderUI({
        switch(input$tab,
            review = tagList(
                selectInput('protein', 'Select Protein', selected = '', choices=c('', sort(unique(as.character(review_data$target_name))))), #input$protein
                div(
                    id = 'form',
                    # Ligand/Xtal Select????
                    textInput('name', 'FedID', fedid), #input$name
                    selectInput('ligand', 'Ligand', selected='', choices = rownames(isolate(r1())), multiple=FALSE), #input$ligand
                    selectInput("decision", "Decision", choices = possDec), #input$decision
                    selectizeInput("reason", "Reason(s)", list(), multiple=TRUE, options=list(create=TRUE)), #input$reason
                    textInput('comments', 'Additional  Comments', value = "", width = NULL,placeholder = NULL), #input$comments
                    fluidRow(
                        column(6, actionButton('submit', 'Submit', class = 'btn-primary')), #input$submit
                        column(6, actionButton('clear', 'Clear', class = 'btn-primary')) #input$clear
                    )
                ),
                fluidRow(
                    column(4, checkboxInput('out4', 'Comp Chem Ready', value = TRUE)), #input$out4
                    column(4, checkboxInput('out5', 'Deposition Ready', value = FALSE)) #input$out5
                ),
                fluidRow(
                    column(4, checkboxInput('out6', 'Deposited', value = FALSE)) #input$out6
                )
            ),
            aqz = tagList(
                selectInput('aq_protein', 'Select Protein', selected = '', choices=fragfolders), #input$aq_protein
                selectInput('aq_ligand', 'Ligand', selected='', choices = list(), multiple=FALSE) #input$aq_ligand
            ),
            fragview = tagList(
                    selectInput('fragSelect', 'Project Select', selected = '', choices=fragfolders), #input$fragSelect
                    checkboxInput('desync', 'Turn off automatic Updates', value = FALSE), #input$desync
                    actionButton('goback', 'Prev Ligand'), #input$goback
                    actionButton('gonext', 'Next Ligand'), #input$gonext
                    selectInput('goto', 'Go to Ligand', choices=list()), #input$goto
                    textInput('crysname', 'Ligand Name', '' ), #input$crysname
                    textInput('smiles', 'Smiles String', ''), #input$smiles
                    textInput('new_smiles', 'New Smiles String', ''), #input$new_smiles
                    textInput('alternate_name', 'Alternate Fragment Name', ''), #input$alternate_name
                    selectizeInput('site_name', 'Site Label', list(), multiple=FALSE, options=list(create=TRUE)), #input$site_name
                    textInput('pdb_entry', 'PDB Entry', ''),#input$pd_entry
                    # textInput('fv_status', 'Status', ''), #input$fv_status, currently unused, will be important someday.
                    textOutput('metastatus'),
                    uiOutput('writeButton'),
                    actionButton('updateTable', 'Refresh Metadata Table') #input$updateTable
            ),
            help = tagList(),
            launchpad = tagList(),
            summary = tagList()
        )
    })

    observeEvent(input$decision,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'decision', 'Decision Selected'))
        possAns <- possRes[[input$decision]]
        updateSelectizeInput(session,'reason', choices=possAns)
    })

    # On session init, set control panel values to defaults.
    ngl_control_values <- reactiveValues()
    ngl_control_values$defaults <- loadDefaultParams()

    # Control Panel Listeners
    observeEvent(input$controls, ignoreNULL = TRUE, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'controls', 'Opening control panel'))
        showModal(
            controlPanelModal(
                values = isolate(ngl_control_values$defaults),
                title = 'NGL Viewer Controls'
            )
        )
    })

    observeEvent(input$tab, ignoreNULL=TRUE, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'Switching Tab'))
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', input$tab))
        if(input$tab == 'review'){
            output$nglShiny <- renderNglShiny(
                nglShiny(name = 'nglShiny', list(), width = NULL, height = NULL)
            )
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'Rebinding review table data'))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling restartSessionKeepOptions'))
            inputData <- restartSessionKeepOptions(configuration=configuration, target_list=target_list)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling reactiviseData'))
            r1 <- reactiviseData(inputData=inputData, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling updateMainTable'))
            output$reviewtable <- updateMainTable(r1=r1, pl=100, input=input)
            reviewtableproxy <- DT::dataTableProxy('reviewtable')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling prebuffer_proxy'))
            prebuffer_proxy(proxy=reviewtableproxy, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'Opening Control Panel'))
            showModal(
                controlPanelModal(
                    values = isolate(ngl_control_values$defaults),
                    title = 'As part of setup please confirm NGL Viewer Controls'
                )
            )
        }
        if(input$tab == 'fragview'){
            output$FragViewnglShiny <- renderNglShiny(
                nglShiny(name = 'nglShiny', list(), width=NULL, height=100)
            )
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'Rebinding fragview table data'))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling react_fv_data'))
            fragview_input <- react_fv_data(data=fragview_data, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling react_fv_data2'))
            fragview_table_data <- react_fv_data2(data=fragview_data, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'updateMainTable2'))
            output$therow <- updateMainTable2(r1=fragview_input, pl=100, input=input)
            fragviewproxy <- DT::dataTableProxy('therow')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling prebuffer_proxy'))
            prebuffer_proxy(proxy=fragviewproxy, input=input)
        }
        if(input$tab == 'aqz'){
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'Rebinding AQ table data'))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling restartSessionKeepOptions'))
            inputData <- restartSessionKeepOptions(configuration=configuration, target_list=target_list)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling reactiviseData'))
            r1 <- reactiviseData(inputData=inputData, input=input)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling updateAQPTable'))
            output$aqp <- updateAQPTable(r=r1, pl=100, input=input)
            aqpproxy <- DT::dataTableProxy('aqp')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'tab', 'calling prebuffer_proxy'))
            prebuffer_proxy(proxy=aqpproxy, input=input)
        }
    })


    observeEvent(input$aq_buster, ignoreNULL=TRUE,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_buster', 'Attempting to Open Buster Report'))
     	pdf_files = list.files(sessionlist$xtalroot, rec=T, pattern='report.pdf', full=T)
    	pdf_file = tail(pdf_files,1)
        to_file <- file.path(configuration$xcr_path, 'www', sprintf('report_%s.pdf', sID))
        copied <- file.copy(from=pdf_file, to=to_file, overwrite=TRUE)
        addResourcePath("www", file.path(configuration$xcr_path, 'www'))
    	output$buster_pdf <- renderUI({
      		tags$iframe(style="height:800px; width:100%", src=sprintf('www/report_%s.pdf', sID))
    	})
    	showModal(
    		customDraggableModalDialog(title=pdf_file,
                if(length(pdf_file)>0) uiOutput("buster_pdf")
                else 'Unable to find Buster Report', size='l', easyClose=FALSE)
   	    )       
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_buster', 'Closing Buster Report'))
    })

    observeEvent(input$buster, ignoreNULL=TRUE, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'buster', 'Attempting to Open Buster Report'))
    	pdf_files = list.files(sessionlist$xtalroot, rec=T, pattern='report.pdf', full=T)
    	pdf_file = tail(pdf_files,1)
        message(pdf_file)
        to_file <- file.path(configuration$xcr_path, 'www', sprintf('report_%s.pdf', sID))
        copied <- file.copy(from=pdf_file, to=to_file, overwrite=TRUE)
        message(copied)
        addResourcePath("www", file.path(configuration$xcr_path, 'www'))
    	output$buster_pdf <- renderUI({
      		tags$iframe(style="height:800px; width:100%", src=sprintf('www/report_%s.pdf', sID))
    	})
    	showModal(
    		customDraggableModalDialog(title=pdf_file,
                if(length(pdf_file)>0) uiOutput("buster_pdf")
                else 'Unable to find Buster Report', size='l', easyClose=FALSE)
   	    )
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'buster', 'Closing Buster Report'))
    })

    observeEvent(input$interactions, ignoreNULL=TRUE,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'interactions', 'Attempting to Open Interaction Graph'))
        pliphtml = gsub('.mol', '_plip.html', sessionlist$mol_file)
        to_file <- file.path(configuration$xcr_path, 'www', sprintf('plip_%s.html', sID))
        copied <- file.copy(from=pliphtml, to=to_file, overwrite=TRUE)
        addResourcePath("www", file.path(configuration$xcr_path, 'www'))
        output$plipview <- renderUI({
		    tags$iframe(style="height:800px; width:100%", src=sprintf('www/plip_%s.html', sID))
    	})
        showModal(
            customDraggableModalDialog(title = 'Protein Ligand Interactions',
                if(file.exists(pliphtml)) uiOutput("plipview")
                else 'Unable to find Protein-Ligand Interactions (plot not generated)', size='l', easyClose=FALSE)
        )
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'interactions', 'Closing Interaction Graph'))
    })

    observeEvent(input$updateParams, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'updateParams', 'Updating Parameters'))
        removeModal()
        for(i in names(ngl_control_values$defaults)){
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'updateParams', i))
            ngl_control_values$defaults[[i]] <- input[[i]]
        }
    })

    observeEvent(input$backgroundColor, { 
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'backgroundColor', 'Updating Background Colour...'))
        updateParam(which='backgroundColor', what=as.character(input$backgroundColor), session=session)
    })
    observeEvent(input$cameraType, { 
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'cameraType', 'Updating camera Type...'))
        updateParam(which='cameraType', what=as.character(input$cameraType), session=session)
    })
    observeEvent(input$mousePreset, { 
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'mousePreset', 'Updating Mouse Preset...'))
        updateParam(which='mousePreset', what=as.character(input$mousePreset), session=session)
    })
    observeEvent(input$clipDist, { 
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'clipDist', 'Updating clip Distance...'))
        updateParam(which='clipDist', what=as.character(input$clipDist), session=session)
    })
    observeEvent(input$fogging, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'fogging', 'Updating fogging...'))
        updateParam(which='fogNear', what=as.character(input$fogging[1]) , session=session)
        updateParam(which='fogFar' , what=as.character(input$fogging[2]) , session=session)
    })
    observeEvent(input$clipping, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'clipping', 'Updating clipping...'))
        updateParam(which='clipNear', what=as.character(input$clipping[1]) , session=session)
        updateParam(which='clipFar' , what=as.character(input$clipping[2]), session=session )
    })


    # NGL Shiny Stages... (all share the same work...) 
    output$nglShiny <- renderNglShiny(
        nglShiny(name = 'nglShiny', list(), width = NULL, height = NULL)
    )

    output$FragViewnglShiny <- renderNglShiny(
        nglShiny(name = 'nglShiny', list(), width=NULL, height=100)
    )

    output$AVnglShiny <- renderNglShiny(
        nglShiny(name = 'nglShiny', list(), width=NULL, height=100)
    )

    # Map Listeners
    observeEvent(input$eventMap,   {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'eventMap', 'Toggle Eventmap'))
        updateVisability(name='eventmap', bool=input$eventMap, session=session)
    })
    observeEvent(input$twofofcMap, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'twofofcMap', 'Toggle 2fofc'))
        updateVisability(name='twofofc' , bool=input$twofofcMap, session=session)
    })
    observeEvent(input$fofcMap,    {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'fofcMap', 'Toggle fofc map'))
        updateVisability(name='fofcpos', bool=input$fofcMap, session=session) 
        updateVisability(name='fofcneg', bool=input$fofcMap, session=session) 
    })

    observeEvent(input$aq_eventMap,   {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_eventMap', 'Toggle Eventmap'))
        updateVisability(name='eventmap', bool=input$aq_eventMap, session=session)
    })
    observeEvent(input$aq_twofofcMap, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_twofofcMap', 'Toggle 2fofc'))
        updateVisability(name='twofofc' , bool=input$aq_twofofcMap, session=session)  })
    observeEvent(input$aq_fofcMap,    {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_fofcMap', 'Toggle fofc'))
        updateVisability(name='fofcpos', bool=input$aq_fofcMap, session=session) 
        updateVisability(name='fofcneg', bool=input$aq_fofcMap, session=session)
    })



    observeEvent(input$isoEvent, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'isoEvent', 'Changing isolevel'))
        updateDensityISO(name='eventmap', isolevel=input$isoEvent, session=session)
    })
    observeEvent(input$iso2fofc, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'iso2fofc', 'Changing isolevel'))
        updateDensityISO(name='twofofc', isolevel=input$iso2fofc, session=session)
    })
    observeEvent(input$isofofc , {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'isofofc', 'Changing isolevel'))
        updateDensityISO(name='fofcpos', isolevel=input$isofofc, session=session) 
        updateDensityISO(name='fofcneg', isolevel=input$isofofc, session=session) 
    })

    observeEvent(input$aq_isoEvent, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_isoEvent', 'Changing isolevel'))
        updateDensityISO(name='eventmap', isolevel=input$aq_isoEvent, session=session)
    })
    observeEvent(input$aq_iso2fofc, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_iso2fofc', 'Changing isolevel'))
        updateDensityISO(name='twofofc', isolevel=input$aq_iso2fofc, session=session)
    })
    observeEvent(input$aq_isofofc , {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_isofofc', 'Changing isolevel'))
        updateDensityISO(name='fofcpos', isolevel=input$aq_isofofc, session=session)
        updateDensityISO(name='fofcneg', isolevel=input$aq_isofofc, session=session)
    })


    observeEvent(input$boxsize , {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'boxsize', 'Changing boxsize'))
        for(windowname in c('eventmap', 'twofofc', 'fofcpos', 'fofcneg')){
            debugMessage(sID=sID, sprintf('(input$%s) Changing box size for: %s', 'boxsize', windowname))
            updateDensityBoxSize(name=windowname, boxsize=input$boxsize, session=session)
        }
    })

    observeEvent(input$fitButton, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'fitButton', 'calling js focus_on_mol'))
        try(session$sendCustomMessage(type='focus_on_mol', list()), silent=TRUE)
    })

    observeEvent(input$aq_fitButton, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_fitButton', 'calling js focus_on_mol'))
        try(session$sendCustomMessage(type='focus_on_mol', list()), silent=TRUE)
    })

    observeEvent(input$asuSwitch, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'asuSwitch', 'Swapping Assembly'))
        debugMessage(sID=sID, sprintf('(input$%s) Changing to %s', 'asuSwitch', input$asuSwitch ))
        try(session$sendCustomMessage('updateAssembly', list(isolate(input$asuSwitch))))
    })

    output$isoEventSlider <- renderUI({
            sliderInput("isoEvent", "",
                    min = 0, max = 3,
                    value = 1, step = 0.1) #input$isoEvent
    })

    output$iso2fofcSlider <- renderUI({
            sliderInput("iso2fofc", "",
                    min = 0, max = 3,
                    value = 1.5, step = 0.1) #input$iso2fofc
    })

    output$isofofcSlider <- renderUI({
            sliderInput("isofofc", "",
                min = 0, max = 3,
                value = 3, step = 0.1) #input$isofofc
    })

    observeEvent(input$aqp_rows_selected,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_rows_selected', 'Preparing to Render protein + mol'))
        rdat <- r1()[input$aqp_rows_selected,,drop=TRUE]
        sessionlist$rowname <- rownames(r1())[input$aqp_rows_selected]
        sessionlist$lig_name <- rdat$ligand_name
        sessionlist$ligand_id <- rdat$ligand_id
        #sessionlist$lig_id <- rdat[[1]][1]
        sessionlist$apo_file <- rdat$apo_pdb
        sessionlist$mol_file <- rdat$lig_mol_file
        sessionlist$event <- rdat$pandda_event
        sessionlist$twofofc_file <- rdat$two_fofc
        sessionlist$fofc_file <- rdat$fofc
        sessionlist$target_name <- rdat$target_name
        sessionlist$res <- rdat$res
        sessionlist$r_free <- rdat$r_free
        sessionlist$r_cryst <- rdat$rcryst
        sessionlist$ramachandran_outliers <- rdat$ramachandran_outliers
        sessionlist$rmsd_angles <- rdat$rmsd_angles
        sessionlist$rmsd_bonds <- rdat$rmsd_bonds
        sessionlist$pdb_latest <- rdat$pdb_latest
        sessionlist$xtalroot <- dirname(dirname(rdat$pdb_latest))
        updateSelectInput(session, 'aq_ligand', choices=isolate(sessionlist$lig_name), selected = isolate(sessionlist$lig_name))
    })

    observeEvent(input$reviewtable_rows_selected, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'reviewtable_rows_selected', 'Preparing to Render protein + mol'))
        rdat <- r1()[input$reviewtable_rows_selected,,drop=TRUE]
        sessionlist$rowname <- rownames(r1())[input$reviewtable_rows_selected]
        sessionlist$lig_name <- rdat$ligand_name
        #sessionlist$lig_id <- rdat[[1]][1]
        sessionlist$apo_file <- rdat$apo_pdb
        sessionlist$mol_file <- rdat$lig_mol_file
        sessionlist$event <- rdat$pandda_event
        sessionlist$twofofc_file <- rdat$two_fofc
        sessionlist$fofc_file <- rdat$fofc
        sessionlist$target_name <- rdat$target_name
        sessionlist$res <- rdat$res
        sessionlist$r_free <- rdat$r_free
        sessionlist$r_cryst <- rdat$rcryst
        sessionlist$ramachandran_outliers <- rdat$ramachandran_outliers
        sessionlist$rmsd_angles <- rdat$rmsd_angles
        sessionlist$rmsd_bonds <- rdat$rmsd_bonds
        sessionlist$pdb_latest <- rdat$pdb_latest
        sessionlist$xtalroot <- dirname(dirname(rdat$pdb_latest))
        updateSelectInput(session, 'ligand', choices=isolate(sessionlist$lig_name), selected = isolate(sessionlist$lig_name))
    })

    observeEvent(input$views, {
        debugMessage(sID=sID, sprintf('(input$%s) Switching view to... %s', 'views', input$views))
        waiter_show(id='nglShiny', html = waiting_screen, color = scales::alpha("black",.5))
        sessionlist$isotype='value'
        if(is.null(input$views)) updateRadioButtons(session, 'views', selected = 'aligned')
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'calling js setup'))
        session$sendCustomMessage(type = 'setup', message = list())
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'preparing stage...'))
        updateParam(which='mousePreset', what=as.character(input$mousePreset), session=session)
        updateParam(which='clipDist', what=as.character(input$clipDist), session=session)
        updateSelectInput(session, 'emap', choices = c('NotAMap.ccp4'), selected = c('NotAMap.ccp4'))
        updateSelectInput(session, 'asuSwitch', selected='AU', choices=c('AU', 'UNITCELL', 'SUPERCELL'))

        the_pdb_file <- isolate(sessionlist$apo_file)
        the_mol_file <- isolate(sessionlist$mol_file)
        the_emaps <- dir(dirname(isolate(sessionlist$apo_file)), pattern='event', full=TRUE)
        the_2fofc_map <- isolate(sessionlist$twofofc_file)
        the_fofc_map <- isolate(sessionlist$fofc_file)

        if(input$renderMisc & !isolate(sessionlist$apo_file) == ""){
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'Rending misc files and images'))
            #spfile <- tail(dir(isolate(sessionlist$xtalroot), pattern='A-1101.png', full.names=T, rec=T),1)
            #output$spiderPlot <- renderImage({
            #    if(length(spfile) == 1){
            #        list(src = spfile, contentType = 'image/png', width=200, height=200)
            #    } else {
            #        list(src = '', contentType = 'image/png', width=200, height=200)
            #    }
            #}, deleteFile=FALSE)
            ligfile <- tail(dir(sprintf('%s/compound', isolate(sessionlist$xtalroot)), pattern = '.png', full.names=T),1)
            renderedligfile <- gsub('.mol', '.png', isolate(sessionlist$mol_file))
            debugMessage(sID=sID, sprintf('rl_image: %s', renderedligfile))
            output$rlimage <- renderImage({
                if(length(renderedligfile) == 1){
                    list(src = renderedligfile, contentType = 'image/png', width=200, height=200)
                } else {
                    list(src = '', contentType = 'image/png', width=200,height=200)
                }
            }, deleteFile=FALSE)
            output$ligimage <- renderImage({
                if(length(ligfile) == 1){
                    list(src = ligfile,contentType = 'image/png',width=200,height=200)
                } else {
                    list(src = '',contentType = 'image/png',width=200,height=200)
                }
            }, deleteFile=FALSE)
            output$ligimage2 <- renderImage({
                if(length(ligfile) == 1){
                    list(src = ligfile,contentType = 'image/png',width=200,height=200)
                } else {
                    list(src = '',contentType = 'image/png',width=200,height=200)
                }
            }, deleteFile=FALSE)
        }
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'Trying to load files'))
        withProgress(message = sprintf('Loading %s Ligand', input$views), value = 0,{
            if(! isolate(sessionlist$apo_file) == ""){
                incProgress(.2, detail = 'Uploading Crystal + Ligand')
                switch(input$views,
                    ' ' = {
                        the_pdb_file <- ''
                        the_mol_file <- ''
                        the_emaps <- ''
                        the_2fofc_map <- ''
                        the_fofc_map <- ''
                    },
                    'aligned' = {
                        # Default Behaviour do not change anything!
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'Viewing mol in Aligned view...'))
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'pdb...'))
                        try(uploadApoPDB3(filepath=the_pdb_file, repr='line', focus=input$autocenter, molfile=isolate(sessionlist$mol_file),session=session), silent=F)
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'contacts...'))
                        try(addContacts(filepath=gsub('_apo', '_bound', the_pdb_file), session=session), silent=TRUE)
                        # Add stuff here:
                        debugMessage(sID=sID, sprintf('Render others?'))
                        clearWindowField(id='othermol', session=session)
                        glob = sprintf('%s*/*.mol', rsplit(dirname(the_mol_file), '_')[1])
                        other_mols = Sys.glob(glob)
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'other mols...'))
                        if(length(other_mols) > 1){
                            other_mols <- other_mols[!other_mols %in% the_mol_file]
                            for(i in other_mols){
                                debugMessage(sID=sID, sprintf('Rendering: %s', i))
                                uploadMolNoFocus(filepath=i, color='pink', session=session)
                            }
                        }
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'mol file...'))
                        try(uploadMolAndFocus3(filepath=the_mol_file, ext='mol', focus=input$autocenter, session=session), silent=F)
                        session$sendCustomMessage(type = 'restore_camera_pos', message = list())
                    },
                    'unaligned' = {
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'View mol in Unaligned view...'))
                        session$sendCustomMessage(type = 'save_camera_pos', message = list())
                        the_pdb_file <- gsub('pipeline_staging', 'pipeline_unaligned', the_pdb_file)
                        the_mol_file <- gsub('pipeline_staging', 'pipeline_unaligned', the_mol_file)
                        the_emaps <- dir(dirname(the_pdb_file), pattern='event', full=TRUE)
                        the_2fofc_map <- gsub('pipeline_staging', 'pipeline_unaligned', the_2fofc_map)
                        the_fofc_map <- gsub('pipeline_staging', 'pipeline_unaligned', the_fofc_map)
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'pdb......'))
                        try(uploadApoPDB3(filepath=the_pdb_file, repr='line', focus=TRUE, molfile=isolate(sessionlist$mol_file), session=session), silent=F)
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'contacts...'))
                        try(addContacts(filepath=gsub('_apo', '_bound', the_pdb_file), session=session), silent=TRUE)
                        # Add stuff here:
                        clearWindowField(id='othermol', session=session)
                        glob = sprintf('%s*/*.mol', rsplit(dirname(the_mol_file), '_')[1])
                        other_mols = Sys.glob(glob)
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'Other Mols...'))
                        if(length(other_mols) > 1){
                            other_mols <- other_mols[!other_mols %in% the_mol_file]
                            for(i in other_mols){
                                debugMessage(sID=sID, sprintf('Rendering: %s', i))
                                uploadMolNoFocus(filepath=i, color='pink', session=session)
                            }
                        } 
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'Mol!...'))
                        try(uploadMolAndFocus3(filepath=the_mol_file, ext='mol', focus=TRUE, session=session), silent=F)
                        
                    },
                    'crystallographic' = {
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'View mol in Crystallographic View...'))
                        session$sendCustomMessage(type = 'save_camera_pos', message = list())
                        splitted <-  rsplit(the_pdb_file, '/')
                        lig_pdb_file <- gsub('_apo', '', the_pdb_file)
                        the_folder <- dirname(gsub('aligned', 'crystallographic', splitted[1]))
                        the_xtal_name <- gsub('_[0-9][A-Z]_apo.pdb', '', splitted[2])
                        the_pdb_file <- sprintf('%s/%s.pdb', the_folder, the_xtal_name)
                        focus_point <- getResNum(pdb=lig_pdb_file)
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'pdb...'))
                        try(uploadPDB(filepath=the_pdb_file, input=input, focus_point=focus_point, session=session), silent=T)
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'contacts...'))
                        try(addContacts(filepath=the_pdb_file, session=session), silent=TRUE)
                        the_2fofc_map <- sprintf('%s/%s_2fofc.map', the_folder, the_xtal_name)
                        the_fofc_map <- sprintf('%s/%s_fofc.map', the_folder, the_xtal_name)
                        the_emaps <- dir(the_folder, pattern=sprintf('%s_event', the_xtal_name), full=TRUE)
                        sessionlist$isotype='sigma'
                    }
                )

                names(the_emaps) <- basename(the_emaps)
                sessionlist$current_emaps <- the_emaps
                print(the_emaps)
                render_bfactor <- switch(input$tab, 
                    'review' = isolate(input$bfactor),
                    'aqz' = isolate(input$aq_bfactor)
                )
                if(render_bfactor){
                    uploadBFactors(filepath=sessionlist$apo_file, clear=TRUE, session=session)
                    updateVisability(name='mol', bool=FALSE, session=session) 
                    uploadBFactors(filepath=gsub('.mol', '.pdb', sessionlist$the_mol_file), clear=FALSE, session=session)
                }

                incProgress(.2, detail = 'Uploading Event map')
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'loading Event Map'))
                updateSelectInput(session, 'emap', choices = names(isolate(sessionlist$current_emaps)), selected = names(isolate(sessionlist$current_emaps))[1])
                # Move this to a different part?
                message('Upload fofcs')
                incProgress(.2, detail = 'Uploading 2fofc map')
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'loading 2fofc map'))
                try(uploadVolumeDensity(filepath=the_2fofc_map,
                    color = 'blue', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$iso2fofc, visable=input$twofofcMap, windowname='twofofc', isotype=sessionlist$isotype, session=session), silent=T)
                incProgress(.1, detail = 'Uploading fofc map')
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'loading pos fofc map'))
                try(uploadVolumeDensity(filepath=the_fofc_map,
                    color = 'lightgreen', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcpos', isotype=sessionlist$isotype, session=session), silent=T)
                incProgress(.1, detail = 'Uploading fofc map')
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'views', 'loading neg fofc map'))
                try(uploadVolumeDensity(filepath=the_fofc_map,
                    color = 'tomato', negateiso = TRUE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcneg', isotype=sessionlist$isotype, session=session), silent=T)
            }
            setProgress(1)
        })
        waiter_hide()
    })

    observeEvent(input$aq_ligand, ignoreNULL=TRUE,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_ligand', 'Trying to load files'))
        waiter_show(id='nglShiny', html = waiting_screen, color = scales::alpha("black",.5))
        sessionlist$isotype <- 'value'
        atomstoquery$data <- data.frame(name=character(),
                 index=character(),
                 comment=character(),
                 stringsAsFactors=FALSE)
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data)}, server=FALSE)
        session$sendCustomMessage(type = 'setup', message = list())
        updateParam(which='mousePreset', what=as.character(input$mousePreset), session=session)
        the_pdb_file <- isolate(sessionlist$apo_file)
        the_mol_file <- isolate(sessionlist$mol_file)
        the_emaps <- dir(dirname(isolate(sessionlist$apo_file)), pattern='event', full=TRUE)
        the_2fofc_map <- isolate(sessionlist$twofofc_file)
        the_fofc_map <- isolate(sessionlist$fofc_file)

        withProgress(message = sprintf('Loading %s Ligand', 'Aligned'), value = 0,{
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_ligand', 'pdb file'))
            try(uploadApoPDB3(filepath=the_pdb_file, repr='line', focus=input$autocenter, molfile=isolate(sessionlist$mol_file), session=session), silent=F)   
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_ligand', 'contacts'))
            try(addContacts(filepath=gsub('_apo', '_bound', the_pdb_file), session=session), silent=TRUE)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_ligand', 'mol file'))
            try(uploadMolAndFocus(filepath=the_mol_file, ext='mol', focus=input$autocenter, session=session), silent=F)
            names(the_emaps) <- basename(the_emaps)
            sessionlist$current_emaps <- the_emaps
            # aq bfactor?
            if(input$tab == 'aqz' & input$aq_bfactor){
                uploadBFactors(filepath=sessionlist$apo_file, clear=TRUE, session=session)
                updateVisability(name='mol', bool=FALSE, session=session) 
                uploadBFactors(filepath=gsub('.mol', '.pdb', sessionlist$mol_file), clear=FALSE, session=session)
            }
            incProgress(.1, detail = 'Uploading Event map')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_ligand', 'event map'))
            try(uploadVolumeDensity(filepath=the_emaps[1], 
                color='orange', negateiso=FALSE, boxsize = input$boxsize, isolevel = input$aq_isoEvent, visable=input$aq_eventMap, windowname='eventmap', isotype=sessionlist$isotype,session=session), silent=T)
            incProgress(.1, detail = 'Uploading 2fofc map')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_ligand', '2fofc map'))
            try(uploadVolumeDensity(filepath=the_2fofc_map,
                color = 'blue', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$aq_iso2fofc, visable=input$aq_twofofcMap, windowname='twofofc', isotype=sessionlist$isotype, session=session), silent=T)
            incProgress(.1, detail = 'Uploading fofc map')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_ligand', 'pos fofc'))
            try(uploadVolumeDensity(filepath=the_fofc_map,
                color = 'lightgreen', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$aq_isofofc, visable=input$aq_fofcMap, windowname='fofcpos', isotype=sessionlist$isotype, session=session), silent=T)
            incProgress(.1, detail = 'Uploading fofc map')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'aq_ligand', 'new fofc'))
            try(uploadVolumeDensity(filepath=the_fofc_map,
                color = 'tomato', negateiso = TRUE, boxsize = input$boxsize, isolevel = input$aq_isofofc, visable=input$aq_fofcMap, windowname='fofcneg', isotype=sessionlist$isotype, session=session), silent=T)
        })

    })

    observeEvent(input$ligand, ignoreNULL = TRUE, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'Trying to load files'))
        waiter_show(id='nglShiny', html = waiting_screen, color = scales::alpha("black",.5))
        previous = isolate(input$views)
        if(previous == 'aligned'){
            session$sendCustomMessage(type = 'setup', message = list())
            updateParam(which='mousePreset', what=as.character(input$mousePreset), session=session)
            the_pdb_file <- isolate(sessionlist$apo_file)
            the_mol_file <- isolate(sessionlist$mol_file)
            the_emaps <- dir(dirname(isolate(sessionlist$apo_file)), pattern='event', full=TRUE)
            the_2fofc_map <- isolate(sessionlist$twofofc_file)
            the_fofc_map <- isolate(sessionlist$fofc_file)

            if(input$renderMisc & !isolate(sessionlist$apo_file) == ""){
                debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'Render Misc Files'))
                #spfile <- tail(dir(isolate(sessionlist$xtalroot), pattern='A-1101.png', full.names=T, rec=T),1)
                #output$spiderPlot <- renderImage({
                #    if(length(spfile) == 1){
                #        list(src = spfile, contentType = 'image/png', width=200, height=200)
                #    } else {
                #        list(src = '', contentType = 'image/png', width=200, height=200)
                #    }
                #}, deleteFile=FALSE)
                ligfile <- tail(dir(sprintf('%s/compound', isolate(sessionlist$xtalroot)), pattern = '.png', full.names=T),1)
                renderedligfile <- gsub('.mol', '.png', isolate(sessionlist$mol_file))
                debugMessage(sID=sID, sprintf('rl_image: %s', renderedligfile))
                output$rlimage <- renderImage({
                    if(length(renderedligfile) == 1){
                        list(src = renderedligfile, contentType = 'image/png', width=200, height=200)
                    } else {
                        list(src = '', contentType = 'image/png', width=200,height=200)
                    }
                }, deleteFile=FALSE)
                output$ligimage <- renderImage({
                    if(length(ligfile) == 1){
                        list(src = ligfile,contentType = 'image/png',width=200,height=200)
                    } else {
                        list(src = '',contentType = 'image/png',width=200,height=200)
                    }
                }, deleteFile=FALSE)
                output$ligimage2 <- renderImage({
                    if(length(renderedligfile) == 1){
                        list(src = renderedligfile,contentType = 'image/png',width=200,height=200)
                    } else {
                        list(src = '',contentType = 'image/png',width=200,height=200)
                    }
                }, deleteFile=FALSE)
            }
            withProgress(message = sprintf('Loading %s Ligand', input$views), value = 0,{
                if(! isolate(sessionlist$apo_file) == ""){
                    incProgress(.2, detail = 'Uploading Crystal + Ligand')
                    debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'pdb file'))
                    try(uploadApoPDB3(filepath=the_pdb_file, repr='line', focus=input$autocenter, molfile=isolate(sessionlist$mol_file), session=session), silent=F)
                    debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'contacts'))
                    try(addContacts(filepath=gsub('_apo', '_bound', the_pdb_file), session=session), silent=TRUE)
                    
                    # Add stuff here:
                    clearWindowField(id='othermol', session=session)
                    glob = sprintf('%s*/*.mol', rsplit(dirname(the_mol_file), '_')[1])
                    other_mols = Sys.glob(glob)
                    debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'Other Mols'))
                    if(length(other_mols) > 1){
                        other_mols <- other_mols[!other_mols %in% the_mol_file]
                        for(i in other_mols){
                            debugMessage(sID=sID, sprintf('Rendering: %s', i))
                            uploadMolNoFocus(filepath=i, color='pink', session=session)
                        }
                    }
                    debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'The OG mol file'))
                    try(uploadMolAndFocus3(filepath=the_mol_file, ext='mol', focus=input$autocenter,session=session), silent=F)
                    names(the_emaps) <- basename(the_emaps)
                    sessionlist$current_emaps <- the_emaps
                    if(input$tab == 'review' & input$bfactor){
                        uploadBFactors(filepath=sessionlist$apo_file, clear=TRUE, session=session)
                        updateVisability(name='mol', bool=FALSE, session=session) 
                        uploadBFactors(filepath=gsub('.mol', '.pdb', sessionlist$mol_file), clear=FALSE, session=session)
                    }
                    incProgress(.2, detail = 'Uploading Event map')
                    debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'Event Map'))
                    updateSelectInput(session, 'emap', choices = names(isolate(sessionlist$current_emaps)), selected = names(isolate(sessionlist$current_emaps))[1])
                    # Move this to a different part?
                    message('Upload fofcs')
                    incProgress(.2, detail = 'Uploading 2fofc map')
                    debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', '2fofc map'))
                    try(uploadVolumeDensity(filepath=the_2fofc_map,
                        color = 'blue', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$iso2fofc, visable=input$twofofcMap, windowname='twofofc', isotype=sessionlist$isotype,session=session), silent=T)
                    incProgress(.1, detail = 'Uploading fofc map')
                    debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'fofc map'))
                    try(uploadVolumeDensity(filepath=the_fofc_map,
                        color = 'lightgreen', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcpos', isotype=sessionlist$isotype, session=session), silent=T)
                    incProgress(.1, detail = 'Uploading fofc map')
                    try(uploadVolumeDensity(filepath=the_fofc_map,
                        color = 'tomato', negateiso = TRUE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcneg', isotype=sessionlist$isotype, session=session), silent=T)
                    sites_df <- rel_df <- data.frame()
                    if(file.exists(gsub('.mol', '_sites.json', the_mol_file))){
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'Loading Sites json'))
                        js <- jsonlite::read_json(gsub('.mol', '_sites.json', the_mol_file))
                        sites_df <- data.frame(names(js), sapply(js, function(x) paste(unlist(x[['site_id']]), collapse=';')), row.names=NULL)
                        colnames(sites_df) <- c('Site Type', 'Site Indicies')
                    }
                    if(file.exists(gsub('.mol', '_relationships.json', the_mol_file))){
                        debugMessage(sID=sID, sprintf('(input$%s) %s', 'ligand', 'Loading relationships json'))
                        js <- jsonlite::read_json(gsub('.mol', '_relationships.json', the_mol_file))
                        rel_df <- data.frame(names(js), t(sapply(js, unlist)),row.names=NULL)
                        colnames(rel_df) <- c('Related Ligand', 'Similarity', 'COM_Distance', 'Nearest Atom Distance', 'Relationship')
                    }
                    output$site_table <- DT::renderDataTable({DT::datatable(sites_df)}, server=FALSE)
                    output$relationship_table <- DT::renderDataTable({DT::datatable(rel_df)}, server=FALSE)
                }
                setProgress(1)
                residues <- get_residues(the_pdb_file)
                updateSelectInput(session, 'gotores', choices=residues)
                updateSelectizeInput(session, 'highlight_res', choices=residues, selected=input$highlight_res)
                if(!is.null(input$highlight_res)){
                    if(!input$highlight_res == ''){
                    pos <- paste(sapply(strsplit(input$highlight_res, '_'), '[', 2), collapse=', ')
                    try(session$sendCustomMessage(type='highlight_residues', list(pos)))
                }}
            })
        } else {
            # There is a problem with observeEvents not rendering stale references therefore we have to manually the loading if the event state does not change.
            updateRadioButtons(session, 'views', selected = 'aligned')
        }
        waiter_hide()
    })

    observeEvent(input$gotores, ignoreNULL=TRUE, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'gotores', 'Going to Residue...'))
        if(!input$gotores == ''){
            pos <- strsplit(input$gotores, '_')[[1]][2]
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'gotores', pos))
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'gotores', 'calling js go_to_res'))
            try(session$sendCustomMessage(type='go_to_residue', list(pos)))
        }
    })

    observeEvent(input$highlight_res, ignoreNULL=TRUE,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'highlight_res', 'Highlighting Residues'))
        if(!input$highlight_res == ''){
            pos <- paste(sapply(strsplit(input$highlight_res, '_'), '[', 2), collapse=', ')
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'highlight_res', 'calling js highlist_residues'))
            try(session$sendCustomMessage(type='highlight_residues', list(pos)))
        }
    })
    observeEvent(input$emap, ignoreNULL = TRUE, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'emap', 'Uploading Specific Eventmap'))
        sel <- isolate(sessionlist$current_emaps)[input$emap]
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'emap', 'callig UploadVolumeDensity'))
        try(uploadVolumeDensity(filepath=sel,
            color = 'orange', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$isoEvent, visable=input$eventMap, windowname='eventmap', isotype=sessionlist$isotype, session=session), silent=T)
    })

    output$plotElement <- renderUI({
        plotOutput('plottoRender', width = "100%",click = NULL,dblclick = NULL,hover = NULL,brush = NULL,inline = FALSE)
    })

    output$flexPlotElement <- renderUI({
        plotlyOutput('flexplot1')
    })

    flexplotData <- reactive({
        rowidx = as.character(r1()$target_name) == input$fpe_target
        ligand_name = rownames(r1())[rowidx]

        col = reactive({
            d <- event_data("plotly_click")
            vec = as.character(r1()[rowidx, 'decision_str'])
            vec[vec == 'NULL'] <- 'Needs Review'
            vec[is.na(vec)] <- 'Needs Review'
            vec[vec=='NA'] <- 'Needs Review'
            print(vec)
            list(
                col=vec,
                button1 = input$submit,
                button2 = input$ok
            )
        })
        list('Target' = input$fpe_target,
             'xlab' = input$fpex,
             'ylab' = input$fpey,
             'rownames' = rownames(r1())[rowidx],
             'data' = data.frame(
                x=jitter(as.numeric(r1()[rowidx,input$fpex])),
                y=jitter(as.numeric(r1()[rowidx,input$fpey])),
                ligand_name=ligand_name,
                status=col()$col
            )
        )
    })

    observeEvent(event_data("plotly_click"),{
        d <- event_data("plotly_click")
        choice=d$customdata
        if(!is.na(choice)){
            rdat <- r1()[choice,,drop=TRUE]
            sessionlist$lig_name <- rdat$ligand_name
            #sessionlist$lig_id <- rdat[[1]][1]
            sessionlist$apo_file <- rdat$apo_pdb
            sessionlist$mol_file <- rdat$lig_mol_file
            sessionlist$event <- rdat$pandda_event
            sessionlist$twofofc_file <- rdat$two_fofc
            sessionlist$fofc_file <- rdat$fofc
            sessionlist$target_name <- rdat$target_name
            sessionlist$res <- rdat$res
            sessionlist$r_free <- rdat$r_free
            sessionlist$r_cryst <- rdat$rcryst
            sessionlist$ramachandran_outliers <- rdat$ramachandran_outliers
            sessionlist$rmsd_angles <- rdat$rmsd_angles
            sessionlist$rmsd_bonds <- rdat$rmsd_bonds
            sessionlist$pdb_latest <- rdat$pdb_latest
            sessionlist$xtalroot <- dirname(dirname(rdat$pdb_latest))
            updateSelectInput(session, 'ligand', choices=isolate(sessionlist$lig_name), selected = isolate(sessionlist$lig_name))
        }
    })

    plotData <- reactive({
        lapply(c('res', 'r_free', 'rcryst', 'ramachandran_outliers', 'rmsd_angles', 'rmsd_bonds'),
            xformplot, alldata=review_data, extradata=isolate(sessionlist)
        )
    })

    output$plottoRender <- renderPlot({
        hmapbar(data=plotData(), title = isolate(sessionlist$lig_name), target_name=(isolate(sessionlist$target_name)))
    })


    ## Launch Pad Stuff:
    output$launchpad_stuff <- renderUI({
        fluidPage(
            selectInput('lp_selection','Select Target', selected = '', choices=fragfolders),
	        checkboxInput('lp_copymaps', 'Copy MapFiles?', value=TRUE),
            textInput('lp_proposal', 'Proposal Number', value = "", placeholder = 'OPEN or number e.g. 12345'),
            textInput('lp_email', 'Email Address', value = ""),
            #downloadButton("downloadMeta", "Download Metadata"),
            actionButton('lp_launcher', "Upload Data"),
            fluidRow(
                column(12,
                    DT::dataTableOutput('lp_meta')
                )
            )
        )
    })

    observeEvent(input$lp_selection, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'lp_selection', 'Selecting Target to upload'))
        if(!isolate(input$lp_selection) == ''){
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'lp_selection', 'calling createUniqueMetaData'))
        sessionlist$fumeta <- createUniqueMetaData(configuration=configuration, target = isolate(input$lp_selection))
        output$lp_meta <- DT::renderDataTable({
            DT::datatable(
                sessionlist$fumeta, callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
                extensions = "Buttons",
                options = list(
                    dom = 'Bfrtip', buttons = c('pageLength', 'copy', 'csv', 'excel'),
                    pageLength = -1,
                    lengthMenu = list(c(10, 25, 100, -1), c('10', '25', '100','All'))
                )
            )
        }, server=FALSE)
        } else {
            DT::renderDataTable({DT::datatable(data.frame('a'='', 'No Target Selected'= 'Select a Target', stringsAsFactors=FALSE))}, server=FALSE)
        }
    })

    # Upload to fragalysis, this need fixing...
    observeEvent(input$lp_launcher, {
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'lp_launcher', 'Attempting to Upload to Fragalysis'))
        # Check if things are sound...
        lpprop <- isolate(input$lp_proposal)
        if((grepl('^[0-9]{5}$', lpprop) | lpprop == "OPEN")){
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'lp_launcher', 'Attaching atom qualities to mol files'))
            rewriteMols(target=isolate(input$lp_selection), configuration=configuration)
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'lp_launcher', 'Creating FragUploadFolder '))
            sessionlist$fullpath_frag <- createFragUploadFolder(meta=sessionlist$fumeta, target=isolate(input$lp_selection), copymaps=input$lp_copymaps, mtz=mtzzz, configuration=configuration)
            # Upload to stuff???
            debugMessage(sID=sID, sprintf('(input$%s) %s', 'lp_launcher', 'UploadingFragFolder'))
            task <- uploadFragFolder(filepath = sessionlist$fullpath_frag, target = isolate(input$lp_selection), proposal = isolate(input$lp_proposal), email = isolate(input$lp_email), configuration=configuration)
            # Clean up... providing fullpath_frag
            #if(grep('zip', basename(sessionlist$fullpath_frag)) & dirname(dirname(sessionlist$fullpath_frag)) == '/dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/'){
            #    try(system('rm -rf %s', dirname(sessionlist$fullpath_frag)), silent=T)
            #}
            showModal(modalDialog(
                    title = 'Upload appears to have been successful', 
                    sprintf('Please be patient - you should recieve an email when the upload on fragalysis has succeeded. 
                    You can track the progress of the upload here: https://fragalysis.diamond.ac.uk/viewer/upload_task/%s/
                    Please do not resubmit! Things are happening behind the scenes!', task), 
                    footer = modalButton("Dismiss"),
                    easyClose = FALSE
            ))
            debugMessage(sID=sID, sprintf('(input$%s) https://fragalysis.diamond.ac.uk/viewer/upload_task/%s/', 'lp_launcher', task))

        } else {
            showModal(modalDialog(
                title = 'Please fill out the form!', 
                'Clicking upload data will trigger an upload event, to do this we require a suitable proposal. 
                If you want to release the data to the public please input OPEN in the Proposal number box, otherwise enter the 5 digit number from your proposal. 
                E.g. if your proposal is lb12345 you should enter 12345', 
                footer = modalButton("Dismiss"), 
                easyClose = FALSE
            ))
        }
    })

    # Button behaviour to download metadata csv file (non-functional, but placeholder...)
    output$downloadMeta <- downloadHandler(
        filename = function() {return('metadata.csv')},
        content = function(file){
            write.csv(sessionlist$fumeta, file)
        }
    )

    # Button Behaviour for download button (non-functional, but placeholder...)
    output$downloadFragData <- downloadHandler(
        filename = function() {
            message('Downloading')
            print(basename(sessionlist$fullpath_frag))
            return(basename(sessionlist$fullpath_frag))
        },
        content = function(file) {
            print(sessionlist$fullpath_frag)
            print(file)
            file.copy(sessionlist$fullpath_frag, file)
        }
    )

    observeEvent(input$config_target, ignoreNULL = TRUE,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'config_target', 'Fetching pipeline params'))
        # Fetch the params from XCDB
        resp <- fetchPipelineOptions(configuration=configuration, target = isolate(input$config_target))
        updateCheckboxInput(session, 'monomeric', value = as.logical(resp$pl_monomeric))
        updateCheckboxInput(session, 'reduce', value =  as.logical(resp$pl_reduce_reference_frame))
        updateCheckboxInput(session, 'covalent', value = as.logical(resp$pl_reduce_reference_frame))
        updateCheckboxInput(session, 'active', value = as.logical(resp$pl_active))
    })

    observeEvent(input$config_change, ignoreNULL = TRUE,{
        debugMessage(sID=sID, sprintf('(input$%s) %s', 'config_change', 'Updating pipeline parameters...'))
    })

    # Stop Timeout.
    autoInvalidate <- reactiveTimer(10000)
    observe({
        autoInvalidate()
        cat("")
    })
}
