# List of Functions used in XCR...
# Many of the functions that involve uploading an element to the NGL stage
# These will use session$sendCustomMessage, all of which could be greatly optimised, alas I have run out of time.
# All sendCustomMessages are described in: https://github.com/xchem/nglShiny/blob/master/inst/htmlwidgets/nglShiny.js

#' Read a text file and return it as a single line separated by `\n` characters
#'
#' @param file A string corresponding to a filepath. 
#' @return Returns the contents of a file as a single string. 
readTxtToOneLine <- function(file){
    message('(readTxtToOneLine) Reading: ', file)
    return(paste0(readLines(file), collapse='\n'))
}

#' Connect to XCDB using provided configuration
#'
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @return Returns a connection object that can be used to make queries. 
xcdbConnect <- function(configuration){
    message('(xcdbConnect): Connecting to xcdb...')
    con <- dbConnect(
        RMariaDB::MariaDB(), 
        dbname = configuration$db, 
        host=configuration$host, 
        port=configuration$port, 
        user=configuration$user, 
        password=configuration$password
    )
    return(con)
}

#' Fetch pipeline parameters for a specific target
#'
#' pipeline parameters are stored in the target table and prefixed with pl_
#'
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @param target Name of the target to pull from xcdb
#' @return A single row data.frame that contains the contents of the below query.
fetchPipelineOptions <- function(configuration, target){
    con <- xcdbConnect(configuration=configuration)
    on.exit(con)
    message('(fetchPipelineOptions): Fetching Target info')
    res <- dbGetQuery(con, sprintf('SELECT * FROM target WHERE target_name = "%s"', target))
    return(res)
}

updatePipelineOptions <- function(configuration, input){
    con <- xcdbConnect(configuration=configuration)
    message('(updatePipelineOptions): Updating Target info')
    dbExecute(con, sprintf('UPDATE target SET %s WHERE target_name="%s"',
        sprintf("pl_monomeric=\'%s\', pl_reduce_reference_frame=\'%s\', pl_covalent_attachments=\'%s\', pl_active=\'%s\'", 
                as.numeric(input$monomeric), as.numeric(input$reduce), as.numeric(input$covalent), as.numeric(input$active), 
        input$config_target)
    )
    dbDisconnect(con) 
}

#' Fetch targets an authenticated fedid is able to see on ISPyB and XCDB.
#'
#' @param fedid String, corresponding to a fedid provided by shinyproxy/active directory when you log in.
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @return Vector of Targets the user is able to see. 
fetchTargets <- function(fedid, configuration){
    con <- xcdbConnect(configuration=configuration)
    on.exit(con)
    message('(fetchTargets): Fetching Targets for ', fedid)
    # Sort this out with INNER joins?
    pid <- dbGetQuery(con, sprintf('SELECT personId from ispyb.Person WHERE login = "%s"', fedid))
    sessions <- dbGetQuery(con, sprintf('SELECT sessionId from ispyb.Session_has_Person WHERE personId = "%s"', pid))
    sessionstr <- paste(sessions[,1], collapse=',')
    visits <- dbGetQuery(con, sprintf('SELECT proposalId,visit_number from ispyb.BLSession WHERE sessionId IN (%s)', sessionstr))
    propid_str <- paste(unique(visits[,1]), collapse=',')
    props <- dbGetQuery(con, sprintf('SELECT proposalId,proposalCode,proposalNumber from ispyb.Proposal WHERE proposalId IN (%s)', propid_str))
    rownames(props) <- as.character(props[,1])
    vis <- cbind(props[as.character(visits[,1]),-1], vn = visits[,2])
    propvis <- unique(sprintf('%s%s-%s', vis[,1], vis[,2], vis[,3]))
    xcrvis <- dbGetQuery(con,sprintf('SELECT id,visit from soakdb_files'))
    xcrvisit_to_see <- paste(xcrvis$id[xcrvis$visit %in% propvis], collapse=',')
    target_ids <- paste(unique(dbGetQuery(con, sprintf('SELECT target_id FROM crystal WHERE visit_id IN (%s)', xcrvisit_to_see))[,1]), collapse=',')
    target_list <- dbGetQuery(con, sprintf("SELECT target_name from target WHERE id IN (%s)", target_ids))[,1]
    return(target_list)
}

#' Get the IDs of the Atoms from a PDB file for a specific ligand and chain.
#'
#' @param pdb String, Filepath corresponding to pdb file.
#' @param lignum The residue number of the ligand
#' @param chain The chain where the ligand resides.
#' 
#' @return Updates proxy with new data...
getAtomIDs <- function(pdb, lignum, chain){
    message('(getAtomIDs) reading pdb...', pdb)
    pdb <- bio3d::read.pdb(pdb)
	ligtable <- pdb$atom[pdb$atom$resid == 'LIG',]
	if(chain == ''){
		chaintab <- split(ligtable, ligtable$chain)
		restab <- split(chaintab[[chain]], chaintab[[chain]]$resno)[[as.numeric(lignum)+1]]
		output <- rownames(restab)
	} else {
		restab <- split(ligtable, ligtable$resno)[[as.numeric(lignum)+1]]
		output <- rownames(restab)
	}
	return(sprintf('@%s',paste(output, collapse=',')))
}

#' Get Residue Number for a specific atom
#'
#' Well this doesn't seem right... I genuinely forgot what this does huh...
#'
#' @param pdb A string corresponding to file path of a pdb file
#' @return Returns a string of pdb atoms in a `@ 1,2,3,4,` format aparently...
getResNum <- function(pdb){
    message('(getResNum): reading pdb...', pdb)
    pdb <- bio3d::read.pdb(pdb)
    return(sprintf('@%s',paste(rownames(pdb$atom), collapse=',')))
}

#' Write data to the bottom of a mol_file
#'
#' @param mol_file String, corresponding to the filepath of a mol file
#' @param id_str A string to add to the file, usually integers seperated by semicolons
#' @param comment_str A string to add to the file, usually a couple of words seperated by semicolons
#' @param name_str A string to add to the file, usually a ATOM names e.g. [LYS]124.C1 seperated by semicolons
#' @return Writes the _str arguments to the mol_file.
write_to_mol_file <- function(mol_file, id_str, comment_str, name_str){
    message('(write_to_mol_file) reading file...', mol_file)
    lines <- readLines(mol_file)
    # Check if the fields exist... 
    # Add the necessary lines IF needed.
    if(any(lines == '> <BADATOMS>')){
        badid_line <- which(lines == '> <BADATOMS>') + 1
        lines[badid_line] <- id_str
    } else {
        lines <- c(lines, '> <BADATOMS>', id_str)
    }
    if(any(lines == '> <BADCOMMENTS>')){
        badcomment_line <- which(lines == '> <BADCOMMENTS>') + 1
        lines[badcomment_line] <- comment_str
    } else {
        lines <- c(lines, '> <BADCOMMENTS>', comment_str)
    }
    if(any(lines == '> <BADATOMNAMES>')){
        badname_line <- which(lines == '> <BADATOMNAMES>') + 1
        lines[badname_line] <- name_str
    } else {
        lines <- c(lines, '> <BADATOMNAMES>', name_str)
    }
    message('(write_to_mol_file) write to file...', mol_file)
    cat(paste(lines, collapse='\n'), file = mol_file)
}

#' rewriteMol files for a specific target
#'
#' When this function is called, all the mol_files and bad atom qualities for a specific target
#' will be checked and have their bad atom qualities added incase they have been lost during an 
#' update.
#'
#' @param target String, the Name of the target to edit the mol files.
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @return Writes the mol files for an entire target in XCDB.
rewriteMols <- function(target, configuration){ # This need refactoring later......
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))
    message('(rewriteMols) Fetching mol file list...')
    fligands <- dbGetQuery(con, 'SELECT * from fragalysis_ligand')
    fligand_id <- as.character(fligands$id)
    rownames(fligands) <- as.character(fligand_id)
    fligands <- fligands[grep(target, fligands$ligand_name), ]
    ligands <- dbGetQuery(con, sprintf('SELECT * from ligand WHERE fragalysis_ligand_id IN (%s)', paste(rownames(fligands), collapse=',')))
    ligand_id <- ligands$id
    rownames(ligands) <- as.character(ligand_id)
    atoms <- dbGetQuery(con, sprintf('SELECT * from bad_atoms WHERE ligand_id IN (%s)', paste(rownames(ligands), collapse=',')))
    rownames(atoms) <- as.character(atoms$ligand_id)
    ats <- na.omit(atoms[rownames(ligands),])
    towrite <- ats[!ats$atomid == '',]
    message('(rewriteMols) Adding atoms to files...')
    for(j in 1:nrow(towrite)){
        atomstowrite <- towrite[j,]
        lig <- as.character(atomstowrite[, 'ligand_id'])
        mol_file = fligands[as.character(ligands[lig,'fragalysis_ligand_id']),]$lig_mol_file
        id_str <- atomstowrite[1,'atomid']
        comment_str <- atomstowrite[1,'comment']
        name_str <- atomstowrite[1,'atomname']
        id_str2 <- paste0(id_str, collapse=';')
        comment_str2 <- paste0(comment_str, collapse=';')
        name_str2 <- paste0(name_str, collapse=';')
        try(write_to_mol_file(mol_file=mol_file, id_str=id_str2, comment_str = comment_str2, name_str=name_str2), silent=T)
    }
}

#' Split a string starting from the right!
#'
#' Functionally identical to python str.rsplit() method.
#'
#' @param string String to split,
#' @param split_by String of character(s) to split input string by
#' @param n Integer, How many splits do you want to perform from the right.
#' @return A vector of length N+1 split by desired string.
rsplit <- function(string, split_by, n=1){
    spl <- strsplit(string, split_by)[[1]]
    c(paste(unlist(head(spl, length(spl)-n)), collapse=split_by), unlist(tail(spl, n)))
}

#' Get Residue Numbers for a pdb file
#'
#' @param pdb_file String corresponding to a file-path to a pdb file.
#' @return Vector of RES Numbers for all residues in a pdb file RES_1, RES_2 etc.
get_residues <- function(pdb_file){
    message('(get_residues): reading pdb...', pdb_file)
    struc <- try(bio3d::read.pdb(pdb_file), silent=T)
    if(inherits(struc, 'try-error')) return('')
    return(c('', unique(paste(struc$atom$resid, struc$atom$resno, sep='_'))))
}

#' Try Catch on Getting the Review data
#'
#' This function calls the getReviewData within a try to prevent error-states.
#'
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @param target_list List of Targets
#' @return data.Frame of Review data.
trygetReviewData <- function(configuration, target_list){
   out <- try(getReviewData(configuration=configuration, target_list=target_list), silent=T)
   return(out)
}

#' Get Data from Review Responses Table in XCDB
#'
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @param target_list List of Targets to get
#' @return data.frame of review data.
getReviewData <- function(configuration, target_list){ 
    message('(getReviewData): getting Review Data')
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))

    targs <- dbGetQuery(con, 'SELECT * from target')
    target_ids <- paste(targs[targs[,2] %in% target_list, 1], collapse=',')
    
    ligand_data <- dbGetQuery(con, sprintf("SELECT * from ligand WHERE target_id IN (%s)", target_ids))
    lig_crys_ids <- paste(ligand_data[,'crystal_id'], collapse=',')

    ligand_crystal_data <- dbGetQuery(con, sprintf("SELECT id, crystal_name, target_id FROM crystal WHERE id IN (%s)", lig_crys_ids))
    rownames(ligand_crystal_data) <- as.character(ligand_crystal_data$id)
    lig_fl_ids <- paste(ligand_data[,'fragalysis_ligand_id'], collapse=',')

    ligand_fl_data <- dbGetQuery(con, sprintf("SELECT * from fragalysis_ligand WHERE id IN (%s)", lig_fl_ids))
    rownames(ligand_fl_data) <-  as.character(ligand_fl_data$id)

    ligand_refinement_data <- dbGetQuery(con, sprintf("SELECT id, crystal_name_id, r_free, rcryst, ramachandran_outliers, res, rmsd_angles, rmsd_bonds, lig_confidence_string, spacegroup, outcome, cif, pdb_latest, mtz_latest FROM refinement WHERE crystal_name_id IN (%s)", lig_crys_ids))
    rownames(ligand_refinement_data) <- as.character(ligand_refinement_data$crystal_name_id)

    ligand_target_data <- dbGetQuery(con, sprintf("SELECT * FROM target WHERE id IN (%s)", paste(ligand_crystal_data[,'target_id'], collapse=',')))
    rownames(ligand_target_data) <- as.character(ligand_target_data$id)

    atom_target_data <- dbGetQuery(con, sprintf("SELECT * FROM bad_atoms WHERE ligand_id IN (%s)", paste(as.character(ligand_data$id), collapse=',')))
    rownames(atom_target_data) <- as.character(atom_target_data$ligand_id)

    # This needs tweaking... Is this going to be slow??
    q1 <- dbGetQuery(con, sprintf("SELECT * FROM crystal_compound_pairs WHERE crystal_id IN (%s)", paste(as.character(ligand_data$crystal_id), collapse=',')))
    q2 <- dbGetQuery(con, sprintf("SELECT * FROM compound WHERE id IN (%s)", paste(as.character(q1$compound_id), collapse=',')))
    

    # Can this be vectorised???
    ccp <- t(sapply(as.character(ligand_data$crystal_id), function(x){
        ccpids <- q1[x == q1$crystal_id, ]
        cids <- q2[q2$id %in% ccpids$compound_id, ]
        ccpids[ccpids$product_smiles == 'None' | is.na(ccpids$product_smiles),'product_smiles'] <- ''
        return(
            c(
                'SMILES'= paste(cids[,'smiles'], collapse=';'), 
                'product_smiles'=paste(ccpids[,'product_smiles'], collapse=';'), 
                'compound_code'= paste(cids[,'compound_string'], collapse=';')
            )
        )
    }))
    rownames(ccp) <- as.character(ligand_data$crystal_id)

    ligand_response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
    mostrecent <- as.data.frame(t(sapply(split(ligand_response_data, ligand_response_data$ligand_name_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
    if(nrow(mostrecent)>1){
        rownames(mostrecent) <- as.character(mostrecent$ligand_name_id)
    } else {
        mostrecent <- ligand_response_data
        rownames(mostrecent) <- as.character(mostrecent$ligand_name_id)
    }
    output <- cbind(
        'ligand_name' = ligand_fl_data[as.character(ligand_data$fragalysis_ligand_id),2],
        'decision_str' = sapply(mostrecent[as.character(ligand_data$id),4], function(x) ifelse(is.null(x), NA, as.character(x))),
        'Review Date' = sapply(mostrecent[as.character(ligand_data$id),6], function(x) ifelse(is.null(x), NA, format(as_datetime(x), '%d/%m/%Y'))),
        'Review Comment' = sapply(mostrecent[as.character(ligand_data$id),'comment'], function(x) ifelse(is.null(x), NA, as.character(x))),
        'Current State' = sapply(ligand_refinement_data[as.character(ligand_data$crystal_id), 11], function(x) ifelse(x==4,'Comp Chem Ready', ifelse(x==5, 'Deposition Ready', 'Deposited'))),
        ligand_refinement_data[as.character(ligand_data$crystal_id),c(9,10,3,4,5,6,7,8,11,12,13,14)],
        ligand_fl_data[as.character(ligand_data$fragalysis_ligand_id),3:13],
        'target_name' = as.character(ligand_target_data[as.character(ligand_crystal_data[as.character(ligand_data$crystal_id),]$target_id),2]),
        ccp,
        'crystal_name' = ligand_crystal_data[as.character(ligand_data$crystal_id),2],
        'decision_int' = sapply(mostrecent[as.character(ligand_data$id),5], function(x) ifelse(is.null(x), NA, as.character(x))),
        'time_submitted' = sapply(mostrecent[as.character(ligand_data$id),6], function(x) ifelse(is.null(x), NA, format(as_datetime(x), '%Y%m%d%H%M%S'))),
        'ligand_id' = ligand_data$id,
        'crystal_id' = ligand_crystal_data[as.character(ligand_data$crystal_id),1],
        'out_of_date' = as.character(mapply(
            function(x,y){
                ifelse(is.na(x), FALSE, y > x)
            },
            as.numeric(sapply(mostrecent[as.character(ligand_data$id),6], function(x) ifelse(is.null(x), NA, format(as_datetime(x), '%Y%m%d%H%M%S')))),
            as.numeric(ligand_fl_data[as.character(ligand_data$fragalysis_ligand_id),13])
            )),
        'bad_atom_index' = atom_target_data[as.character(ligand_data$id), 'atomid'],
        'bad_atom_name' = atom_target_data[as.character(ligand_data$id), 'atomname'],
        'bad_atom_comment' = atom_target_data[as.character(ligand_data$id), 'comment']
    )
    rownames(output) <- make.names(as.character(output$ligand_name), unique=TRUE)
    message('(getReviewData): ... completed...')
    return(output)
}

#' Get the Fragview data from XCDB
#'
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @param target_list List of Targets
#' @return data.frame of Fragalysis View data from XCDB
getFragalysisViewData <- function(configuration, target_list){
    message('(getFragalysisViewData): getting Fragview Data')
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))

    targs <- dbGetQuery(con, "SELECT * from fragalysis_target") 
    target_ids <- paste(targs$id[targs$target %in% target_list], collapse=',')
    targets <- targs[targs$target %in% target_list,]
    rownames(targets) <- as.character(targets$id)

    fvdat <- dbGetQuery(con, sprintf("SELECT * from fragalysis_ligand WHERE fragalysis_target_id IN (%s)", target_ids))
    md_ids <- paste(fvdat$id, collapse=',')

    md <- dbGetQuery(con, sprintf("SELECT * FROM meta_data WHERE ligand_name_id in (%s)", md_ids))
    rownames(md) <- md$ligand_name_id

    output <- cbind(fvdat, targetname=targets[as.character(fvdat$fragalysis_target_id), 'target'], md[as.character(fvdat$id),])

    rns <- gsub('[.]', '-', make.names(as.character(output$ligand_name), unique=TRUE))
    numbers <- grepl('^[0-9]', output$ligand_name)
    rns[numbers] <-  gsub('^X{1}', '', rns[numbers])
    rownames(output) <- rns
    message('(getFragalysisViewData): ...completed...')
    # Test for file existence and filter table???
    return(output)
}

#' Create a metadata data.frame
#'
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @param target The target to generate the metadata for.
#' @return data.frame consisting of the targets metadata.
createUniqueMetaData <- function(configuration, target){
    message('(createUniqueMetaData): Creating Metadata')
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))

    targs <- dbGetQuery(con, "SELECT * from fragalysis_target") 
    target_ids <- paste(targs$id[targs$target %in% target], collapse=',')

    fvdat <- dbGetQuery(con, sprintf("SELECT * from fragalysis_ligand WHERE fragalysis_target_id IN (%s)", target_ids))
    md_ids <- paste(fvdat$id, collapse=',')

    md <- dbGetQuery(con, sprintf("SELECT * FROM meta_data WHERE ligand_name_id in (%s)", md_ids))
    prot <- target
    md <- md[!(md$site_Label == 'IGNORE' | is.na(md$site_Label)),] # Might be a problem for later...
    meta <- md[grep(paste0(prot, '-'), md$fragalysis_name),c(6,7,3,3,4,2,5)] #Ensure valid targets get included.
    colnames(meta) <- c('crystal_name', 'RealCrystalName', 'smiles', 'new_smiles', 'alternate_name', 'site_name', 'pdb_entry')

    # Dodgy Hack for DLS
    rootf <- file.path(configuration$staging, prot, 'aligned/')
    #rootf = sprintf('/dls/science/groups/i04-1/fragprep/pipeline_staging/%s/aligned/', prot)
    for(i in meta$crystal_name){
        id = meta$crystal_name == i
        smifi = sprintf('%s%s/%s_smiles.txt', rootf, i,i)
        if(!file.exists(smifi)){
            mfile <-  sprintf('%s%s/%s_meta.csv', rootf, i,i)
            meta$smiles[id] <- strsplit(readLines(mfile), ',')[[1]][4]
        } else {
            meta$smiles[id] <- readLines(smifi)[1]
        }
    }
    message('(createUniqueMetaData): ...Completed')
    return(meta)
}


#' Update or Crate a row in the metadata XCDB table
#'
#' @param ligand_name_id Fragalysis Ligand ID
#' @param fragalysis_name Name of the Fragalysis Ligand to be rendered in fragalysis
#' @param original_name Name of the crystal from which the ligand was generated from
#' @param site_label Name of the site
#' @param new_smiles New smiles string to fix things downstream if needed.
#' @param alternate_name Additional name to render in fragalysis
#' @param pdb_id 4 Digit PDB code
#' @param status Currently an unused field to dictate how data should be replace in fragalysis... 
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @return Either updates an exsiting row OR creates a new row in the metadata table in XCDB
updateOrCreateRow <- function(
    ligand_name_id, 
    fragalysis_name, 
    original_name, 
    site_label='', 
    new_smiles='', 
    alternate_name='', 
    pdb_id='', 
    status='', 
    configuration
    ){
    df = data.frame(site_Label=as.character(site_label),
                    new_smiles=as.character(new_smiles),
                    alternate_name=as.character(alternate_name),
                    pdb_id=as.character(pdb_id),
                    fragalysis_name=as.character(fragalysis_name),
                    original_name=as.character(original_name),
                    ligand_name_id=as.character(ligand_name_id),
                    status=as.character(status)
                    )
    con <- xcdbConnect(configuration=configuration)
    id = dbGetQuery(con, sprintf("SELECT id from meta_data WHERE ligand_name_id=%s", ligand_name_id))[1,1]
    dbDisconnect(con)
    if(is.na(id)){ # Append 
        message('(updateOrCreateRow): Creating Row')
        dbAppendWrapper(configuration = configuration, table_name = 'meta_data', data = df)
    } else { # Update
        message('(updateOrCreateRow): Updating Row')
        con <- xcdbConnect(configuration=configuration)
        dbExecute(con, sprintf("UPDATE meta_data SET %s WHERE ligand_name_id=%s",
            sprintf("site_Label=\'%s\', new_smiles=\'%s\', alternate_name=\'%s\', pdb_id=\'%s\', fragalysis_name=\'%s\', original_name=\'%s\', status=\'%s\'", site_label, new_smiles, alternate_name, pdb_id, fragalysis_name, original_name, status),
            ligand_name_id)
        )
        dbDisconnect(con)
    }
}

#' Debug Message Wrapper
#'
#' @param sID ID for the user who is active in the active session.
#' @param text The message to output to std::out
#' @return Prints a message to std::out
debugMessage <- function(sID, text){
    message(sprintf('sid: %s | %s | %s', sID, text, Sys.time()))
}

#' Open the control panel model
#'
#' @param values List of values to update the inputs.
#' @param title Title to be described at the top of the modal.
#' @return Generates the modal UI element to control some behaviour.
controlPanelModal <- function(values, title){
    # Function that opens up a modal dialog to contain accessory ngl controls that are not needed to be accessed immediately.
    customDraggableModalDialog(
        title=title,
        numericInput('boxsize', 'Box Size', value = values$boxsize, min = 0, max = 100, width = '100px'),
        numericInput('clipDist', 'Clipping Distance', value = values$clipDist, min = 0, max = 100, width = '100px'),
        sliderInput('fogging', 'Fogging:', min = 0, max = 100, value = values$fogging),
        sliderInput('clipping', 'Clipping:', min = 0, max = 100, value = values$clipping),
        selectInput('backgroundColor', 'Background Colour', selected = values$backgroundColor, choices = c('black', 'white')),
        selectInput('cameraType', 'Camera Type', selected = values$cameraType, choices = c('orthographic', 'perspective')),
        selectInput('mousePreset', 'Mouse Preset', selected = values$mousePreset, choices = c('coot', 'default', 'pymol')),
        easyClose = FALSE,
        footer = tagList(actionButton('updateParams', 'Update Controls'))
    )
}

#' Change the colour of a specific rank in a raster object.
#'
#' @param raster A raster object
#' @param rank The percntile rank to change the colour of
#' @param colour The colour to change
#' @return Updated Raster object.
changeRasterRank <- function(raster, rank, colour){
    raster[floor(rank*100)] <-  colour
    return(raster)
}

hmapbar <- function(data, title, target_name){

    colfunc <- colorRampPalette(c("red", "white", "blue"))

    globalcolour = '#000000'
    experimentcolour = '#00FF00'
    par(xpd=TRUE)
    plot.new()
    text(x=0.5, y=.99, labels=title)
    text(x=0.05, y = .95, labels = 'Metric')
    text(x=0.5, y = .95, labels = 'Percentile Ranks')
    text(x=0.95, y = .95, labels = 'Value')

    text(x= 0.05, y = .8, labels = 'Resolution')
    text(x= 0.95, y = .8, labels = data[[1]][3])
    # Reverse the Resoluion Ranks...
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), 1-data[[1]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, 1-data[[1]][2], globalcolour)
    rasterImage(legend_image,.15,.75,.9,.85)

    text(x= 0.05, y = .65, labels = 'R Free')
    text(x= 0.95, y = .65, labels = data[[2]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), 1-data[[2]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, 1-data[[2]][2], globalcolour)
    rasterImage(legend_image,.15,.6,.9,.7)

    text(x= 0.05, y = .5, labels = 'R Cryst')
    text(x= 0.95, y = .5, labels = data[[3]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), 1-data[[3]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, 1-data[[3]][2], globalcolour)
    rasterImage(legend_image,.15,.45,.9,.55)

    text(x= 0.05, y = .35, labels = 'Ram Outliers')
    text(x= 0.95, y = .35, labels = data[[4]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), 1-data[[4]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, 1-data[[4]][2], globalcolour)
    rasterImage(legend_image,.15,.3,.9,.4)

    text(x= 0.06, y = .2, labels = 'RMSD Angles')
    text(x= 0.95, y = .2, labels = data[[5]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), 1-data[[5]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, 1-data[[5]][2], globalcolour)
    rasterImage(legend_image,.15,0.15,.9,.25)

    text(x= 0.06, y = .05, labels = 'RMSD Bonds')
    text(x= 0.95, y = .05, labels = data[[6]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), 1-data[[6]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, 1-data[[6]][2], globalcolour)
    rasterImage(legend_image,.15,0,.9,.1)

    text(x=0.15, y=-0.05, 'Worst')
    text(x=0.9, y=-0.05, 'Best')

    legend(x=0.10, y=-0.05, legend=c(sprintf('Percentile Relative to other %s crystals marked as comp. chem ready (or above)', target_name), 'Percentile Relative to crystals marked as comp. chem ready (or above) from all XChem Experiments '), fill=c(experimentcolour , globalcolour), bty='n')
}


#' Create and Populate a folder of data on dls disk to be uploaded to Fragalysis
#'
#' @param meta data.frame of metadata contraining the manifest to upload
#' @param target The Name of the Target/Protein to upload
#' @param copymaps Logical, Copy map files from staging file or root crystal directories (mtz files)
#' @param mtz Matrix of mtz files.
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @return Returns nothing, but there should be a folder in configuration$prep_path
createFragUploadFolder <- function(meta, target, copymaps=FALSE, mtz, configuration){
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Creating Fragalysis Folder", value = 0)
    message('(createFragUploadFolder): Creating Folder')
    prot = target
    protsuffix <- paste(prot, format(Sys.time(), "%Y%m%d_%H%M"), sep='_')
    base_root <- sprintf('%s/%s/',configuration$staging, prot)
    extrafiles1 <- file.path(base_root, 'extra_files')
    rootf <- sprintf('%s/%s',configuration$prep_path, protsuffix)
    basef <- sprintf('%s/%s/%s', configuration$prep_path,protsuffix, prot)
    align_dir <- sprintf('%s/%s/%s/aligned',configuration$prep_path, protsuffix, prot)
    crys_dir <- sprintf('%s/%s/%s/crystallographic',configuration$prep_path, protsuffix, prot)

    message('(createFragUploadFolder): Making Directories')
    dir.create(rootf)
    dir.create(basef)
    dir.create(align_dir)
    dir.create(crys_dir)

    # Maybe use file.copy(...,rec=T)?
    # Make Investigation.cif file here...
    try(system(sprintf('cp -r %s %s', extrafiles, basef)), silent = TRUE)

    # aligned data copy
    progress$set(message = "Copying aligned Files", value = 0)
    increment = (1/nrow(meta))/2.2
    message('(createFragUploadFolder): Copying Files')
    for(frag in meta$crystal_name){
        progress$inc(increment, detail=frag)
        cf <- sprintf('%saligned/%s', base_root, frag)
        nf <- paste(align_dir, frag, sep='/')
        files <- list.files(cf)
	    files2 <- files
        if(copymaps){
		    # Copy event_0 and save it as event.cpp4
            event_0 <- grepl('event_0', files)
            event_19 <- grepl('event_[1-9]', files)
		    if(sum(event_0) == 1) files2[event_0] <- gsub('event_0.ccp4', 'event.ccp4', files[event_0])
		    if(sum(event_19) > 0){
			    files[event_19] <- NA
			    files2[event_19] <- NA
			    files <- na.omit(files)
			    files2 <- na.omit(files2)
		    }
        } else {
		    files <- files[!grepl("(.map$|.ccp4$)", files)]
		    files2 <- files
	    }   
	    system(sprintf('mkdir %s', nf))
        file.copy(file.path(cf,files), file.path(nf, files2))
    }

    # crystallographic copy
    progress$set(message = "Copying Input Files", value = 0.5)
    for(rcn in meta$RealCrystalName){
        progress$inc(increment, detail=rcn)
        cf <- sprintf('%scrystallographic/', base_root)
        #nf <- paste(crys_dir, frag, sep='/')
        nf <- crys_dir
        files <- list.files(cf, pattern=rcn)
	    # Remove maps from crystallographic...
        files <- files[!grepl("(.map$|.ccp4$)", files)]
        print(files)
        file.copy(file.path(cf,files), file.path(nf, files))
        if(copymaps){
		    # Replace uncut maps with mtz files...
		    # IF RCN is in mtz
		    mtz_file = mtz[,1] == rcn
            if(sum(mtz_file)>0){
            	id <- which(mtz_file)[1]
                mtz_file <- mtz[id,2]
                print(mtz_file)
		        file.copy(mtz_file, file.path(nf, sprintf('%s_refine.mtz', rcn)))
		        rootmtz <- dir(dirname(dirname(mtz_file)), pattern='.mtz')
		        event_mtz <- rootmtz[grep('event',rootmtz)]
			    if(length(event_mtz)>0){
				    newnames <- basename(event_mtz)
			        ends <- sapply(strsplit(newnames, '_'), tail, 1)
			        print(file.path(dirname(dirname(mtz_file)), event_mtz))
				    file.copy(file.path(dirname(dirname(mtz_file)), event_mtz), file.path(nf, sprintf('%s_event_%s', rcn, ends)))
			    }
		    }
	    }
    }
    message('(createFragUploadFolder): Copied!')
    # Write and copy stuff
    file.copy(sprintf('%sreference.pdb', base_root), sprintf("%s/reference.pdb", basef))
    write.csv(meta, sprintf("%s/metadata.csv", basef), quote = FALSE)
    write.csv(meta, sprintf("%s/metadata.csv", align_dir), quote = FALSE)

    progress$set(message = "Zipping File!", value = .9)

    # Zip File 
    zipf <- sprintf('%s.zip', prot)
    #zipcommand <- sprintf("(cd %s && zip -r %s . && touch done.log)", rootf, prot)
    #system(zipcommand, wait=FALSE)
    currentwd <- getwd()
    setwd(rootf)
    message('(createFragUploadFolder): Zipping folder')
    zip(zipfile=zipf, files=prot)
    message('(createFragUploadFolder): Zipped!')
    cat('', file='done.log')
    setwd(currentwd)
    system(sprintf('chmod 775 %s', rootf))
    system(sprintf('chmod 775 %s/%s', rootf, prot))
    system(sprintf('chmod 775 %s/%s', rootf, zipf))
    full_path_zipf <- sprintf('%s/%s', rootf, zipf)
    message('(createFragUploadFolder): finished without error')
    return(full_path_zipf)
}


#' Upload Frag Folder to fragalysis using external script
#'
#' @param filepath The file path of the zip file to upload
#' @param target Name of the Target
#' @param proposal The proposal number or OPEN
#' @param email Email eddress for whom to email when process finishes
#' @param configuration List containing various config options including xcdb credentials.
#' Described in configuration.toml
#' @return This function should upload stuff to fragalysis.diamond.ac.uk
uploadFragFolder <-  function(filepath, target, proposal, email, configuration){
    message('(uploadFragFolder): Uploading FragFolder')
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Uploading Fragalysis Folder", value = 0)
    external_script <- file.path(configuration$script_path, 'upload_2_fragalysis.sh')
    message('(uploadFragFolder): Running External Command')
    command <- sprintf('%s %s %s %s %s %s', external_script, basename(filepath), target, proposal, email, filepath)
    message('(uploadFragFolder) :', command)
    task <- tail(system(command, intern=T),1)
    return(task)
}

#' Convert values to percentile ranks
#'
#' @param x Values to rank
#' @return Percentile Ranks
perc.rank <- function(x) trunc(rank(x))/length(x)

#' Not sure what does this does now...
perc.rank2 <- function(x, xo=NULL)  length(x[x <= xo])/length(x)

xformplot <- function(x, alldata, extradata){
    y <- x
    if(x=='rcryst') y <- 'r_cryst'
    currentvalue <- as.numeric(extradata[[y]])
    idx <-  which(as.character(alldata$target_name) == extradata$target_name)
    return(
        c(
            # experiment, global, value
            perc.rank2(x=as.numeric(alldata[idx, x]), xo=currentvalue),
            perc.rank2(x=as.numeric(alldata[,x]), xo=currentvalue),
            currentvalue
        )
    )
}

#' Update the boxsize of a window object in nglShiny
#'
#' @param name Name of the window object to configure
#' @param boxsize Size of the box,
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Updates the boxsize of the desired object   
updateDensityBoxSize <- function(name, boxsize, session) session$sendCustomMessage('updateVolumeDensityBoxSize', list(name, boxsize))

#' Update the isolevel of a window object in nglShiny
#'
#' @param name Name of the window object to configure
#' @param isolevel The isolevel to set to.
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Updates the iso level of the desired object    
updateDensityISO <- function(name, isolevel, session) session$sendCustomMessage('updateVolumeDensityISO', list(name, isolevel))
        

#' Update the visibilty of a window object in nglShiny
#'
#' @param name Name of the window object to configure
#' @param bool Logical value indicating visibility
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Updates the visibility of the window object    
updateVisability <- function(name, bool, session){
    session$sendCustomMessage(
        type = 'updateVolumeDensityVisability',
        list(
            as.character(name),
            tcl(bool)
        )
    )
}

#' Upload a volume density
#'
#' @param filepath File path for the map file (.ccp4) or (.map)
#' @param color Colour of the map
#' @param negateiso Whether or not to invert the negative ISO (only useful to do red fofc map)
#' @param boxsize Initial Size of the box
#' @param visable Logical, Whether or not to see the density 
#' @param windowname Name of the window object in js to add the density to
#' @param isotype Either (sigma or value)
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Read a volume density file, converts it to base64 and then it get transferred for NGL to be rendered.
uploadVolumeDensity <- function(filepath, color, negateiso = FALSE, boxsize, isolevel, visable, windowname, isotype=c('sigma', 'value'), session){
    isotype <- match.arg(isotype)
    message('(uploadVolumeDensity): Reading map file')
    volume_bin <- readBin(filepath, what='raw', file.info(filepath)$size)
    volume_b64 <- base64encode(volume_bin, size=NA, endian=.Platform$endian)
    message('(uploadVolumeDensity): Uploading map file')
    session$sendCustomMessage(
        type = 'addVolumeDensity',
        message = list(
            as.character(volume_b64), #0
            as.character(isolevel),#1
            as.character(color),#2
            tcl(negateiso),#3
            as.character(getExt(filepath)),#4
            as.character(boxsize),#5
            tcl(visable),#6
            as.character(windowname),#7
            as.character(isotype)#8
        )
    )
}

#' Get file Extension
#'
#' @param x The file name
#' @return the extension of the file.
getExt <- function(x) sapply(strsplit(x, '[.]'), tail, 1)

#' upload a mol file without zooming in on it.
#'
#' @param filepath File path of mol file to upload
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Uploads a mol file to the ngl stage.
uploadUnfocussedMol <- function(filepath, session){
    message('(uploadUnfocussedMol) Uploading ...', filepath)
    if(file.exists(filepath)){
        choice <- readTxtToOneLine(file=filepath)
        session$sendCustomMessage(
            type='addMol',
            list(choice)
        )
    } else {
        message('(uploadUnfocussedMol) ',filepath, 'Not a File')
    }
}

#' Upload a mol file, and zoom in on it!
#'
#' @param filepath File path of mol file to upload
#' @param ext File Extension
#' @param focus Logical, indicate whether or not to focus on the mol file when it is loaded
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Uploads a mol file to the ngl stage.
uploadMolAndFocus <- function(filepath, ext, focus, session){
    message('(uploadMolAndFocus) Uploading ...', filepath)
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'addMolandfocus',
        list(choice, ext, tcl(focus))
    )
}

#' Upload bfactors to NGL Stage
#'
#' @param filepath File of pdb file containing bfactors.
#' @param clear Logical, if true, removes existing bfactor render
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Adds Bfactors to NGL stage
uploadBFactors <- function(filepath, clear=TRUE, session){
    if(clear) clearWindowField(id='bfactor', session=session)
    choice <- readTxtToOneLine(file=filepath)
    message('(uploadBFactors) Uploading bfactor')
    session$sendCustomMessage(
        type = 'setBFactor',
        message = list(
            choice
        )
    )
}

#' Upload contact map to NGL Stage
#'
#' @param filepath File of pdb file containing contacts to render
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Adds contact information from a pdb file to NGL Stage
addContacts <- function(filepath, session){
    choice <- readTxtToOneLine(file=filepath)
    message('(addContacts) Uploading contacts')
    session$sendCustomMessage(
        type = 'addContacts',
        message = list(
            choice
        )
    )
}

#' Read Mol File plus bad atom quality and render in a particular fashion in NGL stage
#'
#' @param filepath File of mol file to upload
#' @param ext Extension of the file
#' @param focus Logical, indicate whether or not uploaded mol is focussed on
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return This will take a particular type of mol file with BADATOM comments mentioned and use that information
#' to render spikyballs and stripy bonds!
uploadMolAndFocus3 <- function(filepath, ext, focus, session){
    choice <- readTxtToOneLine(file=filepath)
    ml <- readLines(filepath)
    n <- as.numeric(strsplit(trimws(ml[4]),' ')[[1]][1])
    ban <- unlist(strsplit(ml[grep('> <BADATOMNAMES>', ml) + 1], ';'))
    bid <- unlist(strsplit(ml[grep('> <BADATOMS>', ml) + 1], ';'))
    bic <- unlist(strsplit(ml[grep('> <BADCOMMENTS>', ml) + 1], ';'))
    if(all(!is.na(bid))){
        # check ban for non-HETs
        mol_atoms <- ifelse(is.na(ban), TRUE, grepl('HET', ban))
        # Check if any oob?
        mid <- bid[as.numeric(bid) < n]
        mic <- bic[as.numeric(bid) < n]
    } else {
        mid <- ''
        mic <- ''
    }
    message('(uploadMolAndFocus3) Uploading mol file with bad atoms')
    session$sendCustomMessage(
        type = 'setBadMolandfocus',
        list(choice,ext, tcl(focus), paste(mid, collapse=';'), paste(mic,collapse=';'))
    )
}

#' Read pdb and mol file containing bad atoms and render then in a particular style in ngl.
#'
#' @param filepath Filepath of pdb file
#' @param repr ngl structure representation to use
#' @param focus Logical to indicicate whether or not to focus on uploaded component
#' @param molfile filepath to mol file containing bad atoms
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return This will take a particular type of mol file with BADATOM comments mentioned and use that information
#' to render a pdb file with spikyballs and stripy bonds!
uploadApoPDB3 <- function(filepath, repr, focus, molfile, session){
    choice <- readTxtToOneLine(file=filepath)
    # Read mol-file for badatomids and comments...
    ml <- readLines(molfile)
    n <- as.numeric(strsplit(trimws(ml[4]),' ')[[1]][1])
    ban <- unlist(strsplit(ml[grep('> <BADATOMNAMES>', ml) + 1], ';'))
    bid <- as.numeric(unlist(strsplit(ml[grep('> <BADATOMS>', ml) + 1], ';')))
    bic <- unlist(strsplit(ml[grep('> <BADCOMMENTS>', ml) + 1], ';'))
    if(all(!is.na(bid))){
        # check ban for non-HETs
        protein_atoms <- ifelse(is.na(ban), FALSE, !grepl('HET', ban))
        # Check if any oob?
        pid <- bid[protein_atoms | as.numeric(bid) > n]
        pic <- bic[protein_atoms | as.numeric(bid) > n]
    } else {
        pid <- ''
        pic <- ''
    }
    message('(uploadApoPDB3) Uploading pdb file with bad atoms')
    session$sendCustomMessage(
        type = 'setBadAtomsPDB',
        message = list(
            choice,
            repr,
            tcl(focus),
            paste(pid, collapse=';'),
            paste(pic, collapse=';')
        )
    )
}

#' convert TRUE/FALSE to "true"/"false"
#'
#' @param x Logical Value
#' @return returns the logical value into "true" or "false"
tcl <- function(x) tolower(as.character(as.logical(x))) # or you can do: ifelse(as.logical(x), 'true', 'false')

#' Upload a pdb file with a particular flavour
#'
#' @param filepath filepath to pdb file
#' @param repr nglshiny repr to render pdb file in
#' @param focus boolean to indicate whether or not to focus on newly uploaded object
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Uploads a pdb file to the ngl stage
uploadApoPDB <- function(filepath, repr, focus, session){
    choice <- readTxtToOneLine(file=filepath)
    message('(uploadApoPDB) Uploading pdb file')
    session$sendCustomMessage(
        type = 'setapoPDB',
        message = list(
            choice,
            repr,
            tcl(focus)
        )
    )
}

#' Remove a component that is stored inside the javascript window object
#'
#' @param objectname Window object name to remove
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Removes component from ngl stage
removeNamedComponent <- function(objectname, session) session$sendCustomMessage(type='removeNamedComponent', list(objectname))

#' Upload a pdb file for the nth time...
#'
#' @param filepath Filepath to pdb file
#' @param input Shiny input R6 object contains various information about user input
#' @param focus_point NGL selection language to specify where to focus
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Upload a pdb file...
uploadPDB <- function(filepath, input, focus_point = 'LIG', session){
    choice <- readTxtToOneLine(file=filepath)
    message('(uploadPDB) Uploading pdb file')
    session$sendCustomMessage(
        type = 'setPDB2', # See xchem/NGLShiny for details on setPDB2
        message = list(
            choice,
            isolate(input$clipDist),
            isolate(input$clipping)[1],
            isolate(input$clipping)[2],
            isolate(input$fogging)[1],
            isolate(input$fogging)[2],
            focus_point
        )
    )
}

#' Update a parameter in ngl
#'
#' @param which The option to change 
#' @param what The value to change the option to
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Updates the parameter in ngl
updateParam <- function(which, what, session) session$sendCustomMessage('updateaparam', list(which, what))

#' Pull current parameters from users input
#'
#' @return Returns a list of the users input.
getCurrentParams <- function(input){
    list(
        fogging = input$fogging,
        clipping = input$clipping,
        boxsize = input$boxsize,
        clipDist = input$clipDist,
        backgroundColor = input$backgroundColor,
        cameraType = input$cameraType,
        mousePreset = input$mousePreset
    )
}

#' Create a list with some default parameters
#'
#' @return Returns a list with the following values
loadDefaultParams <- function(){
    list(
        fogging = c(50,62),
        clipping = c(42,100),
        boxsize = 5,
        clipDist = 10,
        backgroundColor = 'black',
        cameraType = 'orthographic',
        mousePreset = 'coot'
    )
}

#' write bad atoms to xcdb
#'
#' @param ligand_id Pthe ligand id from the in the bad_atoms table
#' @param configuration The aforementioned configuration file generated from configuration.toml
#' @param atomstoquery bad atoms data.frame produced by XCR
#' @param sessionlist XCR session information (reactiveValues object)
#' @return appends or updates bad atoms to xcdb
writeAtoms <- function(ligand_id, configuration, atomstoquery, sessionlist){
    newdat <- data.frame(
        ligand_id=ligand_id,
        atomid=paste(as.numeric(atomstoquery$data[,'index']), collapse=';'),
        comment=paste(atomstoquery$data[,'comment'], collapse=';'),
        atomname=paste(atomstoquery$data[,'name'], collapse=';')
    )
    con <- xcdbConnect(configuration=configuration)
    id = dbGetQuery(con, sprintf("SELECT id from bad_atoms WHERE ligand_id=%s", ligand_id))[1,1]
    dbDisconnect(con)
    if(is.na(id)){
        message('(writeAtoms) Appending bad atoms')
        dbAppendWrapper(configuration = configuration, table_name = 'bad_atoms', data = newdat)
    } else {
        message('(writeAtoms) Updating ligand with atom qualities')
        con <- xcdbConnect(configuration=configuration)
        dbExecute(con, sprintf(
            "UPDATE bad_atoms SET %s WHERE ligand_id=%s",
            sprintf("atomid=\'%s\', comment=\'%s\', atomname=\'%s\'", newdat$atomid, newdat$comment, newdat$atomname),
            ligand_id)
        )
        dbDisconnect(con)        
    }
    # Write atom data to .mol file???
    message('(writeAtoms) Updating ligand with atom qualities')
    badnamestr = paste(atomstoquery$data[,'name'], collapse=';')
    badidsstr = paste(atomstoquery$data[,'index'], collapse=';')
    badcommentstr = paste(atomstoquery$data[,'comment'], collapse=';')
    if(!badidsstr == ''){
        lines <- readLines(isolate(sessionlist$mol_file))
        if(any(lines == '> <BADATOMS>')){
            badid_line <- which(lines == '> <BADATOMS>') + 1
            lines[badid_line] <- badidsstr
        } else {
            lines <- c(lines, '> <BADATOMS>', badidsstr)
	    }
        if(any(lines == '> <BADCOMMENTS>')){
            badcomment_line <- which(lines == '> <BADCOMMENTS>') + 1
            lines[badcomment_line] <- badcommentstr
	    } else {
		    lines <- c(lines, '> <BADCOMMENTS>', badcommentstr)
	    }
        if(any(lines == '> <BADATOMNAMES>')){
		    an_line <- which(lines == '> <BADATOMNAMES>') + 1
            lines[an_line] <- badnamestr
        } else {
            lines <- c(lines, '> <BADATOMNAMES>', badnamestr)
        }
        cat(paste(lines, collapse='\n'), file = isolate(sessionlist$mol_file))
        message('(writeAtoms) ...Finished')
    }
}

#' Check XCDB and then create a modal if someone has reviewed a ligand before you
#'
#' @param id ID of ligand
#' @param configuration The aforementioned configuration file generated from configuration.toml
#' @return Makes a modal for the user to interact with.
displayModalWhoUpdated <- function(id, configuration){
    message('(displayModalWhoUpdated)')
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))
    response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
    mostrecent <- as.data.frame(t(sapply(split(response_data, response_data$ligand_name_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
    rownames(mostrecent) <- as.character(mostrecent$ligand_name_id)
    user <- mostrecent[as.character(id), 'fedid'][[1]]
    showModal(
        modalDialog(
            title = "Someone has recently reviewed this ligand!",
            sprintf("A User (%s) has recently reviewed this structure. Restarting the session to update their response. If you disagree with the current response, please submit another response or select another crystal.", user),
            footer = tagList(actionButton("ok", "Okay"))
        )
    )
}

#' Check if users session is later than db entry for a specific ligand
#'
#' @param id id of the ligand
#' @param sessionTime Time when session was last refreshed
#' @param configuration The aforementioned configuration file generated from configuration.toml
#' @return Returns a TRUE or FALSE depending if sessionTime is greater than a timestamp in db.
sessionGreaterThanMostRecentResponse <- function(id, sessionTime, configuration){
    message('(sessionGreaterThanMostRecentResponse)')
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))
    response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
    if(nrow(response_data) > 0){
        mostrecent <- as.data.frame(t(sapply(split(response_data, response_data$ligand_name_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
        rownames(mostrecent) <- as.character(mostrecent$ligand_name_id)
        t0 <- mostrecent[as.character(id), 'time_submitted'][[1]]
        output <- ifelse(is.null(t0), TRUE, sessionTime > t0)
    } else {
        output <- TRUE
    }
    return(TRUE) # LOL. Why does this function exist anymore.
}


#' Upload mol file but dont focus on it
#'
#' @param filepath file path to mol file
#' @param color color of ligand
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Adds molecule to ngl
uploadMolNoFocus <- function(filepath, color, session){
    message('(uploadMolNoFocus)')
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'fv_addMolandfocus_withcolor',
        list(choice, color)
    )
}

#' Clears window.id from javascript
#'
#' @param id Window object id to clear
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Returns nothing, clears up stuff.
clearWindowField <- function(id, session){
    message('(clearWindowField)')
    session$sendCustomMessage(
        type = 'clear_window_field',
        list(id)
    )
}

#' Check if string is blank, NA or NULL
#'
#' @param x String to Check
#' @return Returns TRUE is '', NA or NULL
blankNAorNull <- function(x){
    ifelse(is.null(x), TRUE, is.na(x) | x == '')
}

#' UUpload mol file but do focus on it
#'
#' @param filepath Filepath of mol file
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @return Uploads mol file to nfl
uploadMolAndFocus2 <- function(filepath, session){
    message('(uploadMolNoFocus2)')
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'fv_addMolandfocus',
        list(choice)
    )
}

#' Get fragalysis annotation data from xcdb and make it reactive
#'
#' @param configuration The aforementioned configuration file generated from configuration.toml
#' @param target_list Vector of Targets to fetch
#' @return Updates proxy with new data...
reactivegetFragalysisViewData <- function(configuration, target_list){
    message('(reactivegetFragalysisViewData)')
    reactive({getFragalysisViewData(configuration=configuration, target_list=target_list)})
}

#' Get fragalysis annotation data from xcdb and make it reactive with a slight variation
#'
#' @param data ReactiveData object ideally from reactivegetFragalysisViewData
#' @param input users shiny input object
#' @return Creates a version of the reactive datsa that may or may not be subsetted.
react_fv_data2 <- function(data, input){
    message('(react_fv_data2)')
    reactive({
        if(!is.null(input$fragSelect)){
            if(!input$fragSelect == '' | input$fragSelect == 'Select'){
                toshow <- data()$targetname == input$fragSelect
                dat <- data()[toshow, ]
                rn <- rownames(dat)
                rownames(dat) <- gsub('-[0-9]$', '', rn)
                dat
            } else {
                data()[NULL,]
            }
        }
    })
}

#' Get fragalysis annotation data from xcdb and make it reactive with a slight variation
#'
#' @param data ReactiveData object ideally from reactivegetFragalysisViewData
#' @param input users shiny input object
#' @return Creates a version of the reactive data that may or may not be subsetted. But with specific column names!
react_fv_data <- function(data, input){
    message('(react_fv_data2)')
    reactive({
        if(!is.null(input$fragSelect)){
            if(!input$fragSelect == '' | input$fragSelect == 'Select'){
                toshow <- data()$targetname == input$fragSelect
                dat <- data()[toshow, c('original_name', 'fragalysis_name', 'alternate_name','site_Label', 'new_smiles', 'pdb_id')]
                rn <- rownames(dat)
                rownames(dat) <- gsub('-[0-9]$', '', rn)
                dat
            } else {
                data()[NULL,]
            }
        }
    })
}

#' Hopefull wont break...
#'
#' @param flexdata Don't worry about this...
#' @return Dont worry about it
updateFlexPlot <- function(flexdata){
    message('(updateFlexPlot)')
    renderPlotly({
        plot_ly(flexdata()$data, x=~x, y=~y, text=~ligand_name, color=~status, customdata = ~ligand_name, size=20) %>% config(scrollZoom = TRUE)
    })
}

#' Converts data into a plotly format
#'
#' @param r1 ...
#' @param input ....
#' @return Technically makes the data for a plotly plot
flexPlotDataFun <- function(r1, input){
    message('(flexPlotDataFun)')
    reactive({
        rowidx = as.character(r1()$target_name) == input$fpe_target
        ligand_name = rownames(r1())[rowidx]
        col = reactive({
            d <- event_data("plotly_click")
            vec = as.character(r1()[rowidx, 'decision_str'])
            vec[vec == 'NULL'] <- 'Needs Review'
            vec[is.na(vec)] <- 'Needs Review'
            vec[vec == 'NA'] <- 'Needs Review'
            list(
                col=vec
            )
        })
        list(
            'Target' = input$fpe_target,
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
}

#' Update the main table function v2
#'
#' @param r1 Reactive Data
#' @param pl Pagelength
#' @param input input R6 object
#' @return Returns a datatable with particular formatting
updateMainTable2 <- function(r1, pl=100, input){
    message('(updateMainTable2)')
    if (is.null(isolate(input$therow_state))) {
        dtt <- DT::datatable(
            r1(),
            selection = 'single',
            callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
            options = list(
                stateSave = TRUE, 
                pageLength=pl, 
                columnDefs=list(list(orderable=TRUE, targets=c(0,1,2,3,4,5,6)))
            )
        )
    } else {
        dtt <- DT::datatable(
            r1(),
            selection = 'single',
            callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
            options = list(
                stateSave = TRUE,
                order = isolate(input$therow_state$order),
                paging = TRUE,
                pageLength = isolate(input$therow_state$length), 
                columnDefs=list(list(orderable=TRUE, targets=c(0,1,2,3,4,5,6)))
            )   
        )
    }
    DT::renderDataTable({
        dtt %>% DT::formatStyle(columns = 1:ncol(r1()),"white-space"="nowrap")
    }, server=FALSE)
}

#' updateMainTable version 1
#'
#' @param r1 Reactive data
#' @param pl pagelength
#' @param input input R6 object
#' @return Returns a datatable with slightly different styling compared to updateMainTable2
updateMainTable <- function(r1, pl=100, input){
    message('(updateMainTable)')
    if (is.null(isolate(input$reviewtable_state))) {
        dtt <- DT::datatable(
            r1(),
            selection = 'single',
            callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
            extensions = "Buttons",
            options = list(
                dom = 'Bfrtip', buttons = c('pageLength', 'copy', 'csv', 'excel'),
                stateSave = TRUE, 
                pageLength=pl, 
                pageLength = isolate(input$reviewtable_state$length),
                lengthMenu = list(c(10, 25, 100, -1), c('10', '25', '100','All')),
                columnDefs=list(list(orderable=TRUE))
            )
        )
    } else {
        dtt <- DT::datatable(
            r1(),
            selection = 'single',
            callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
            extensions = "Buttons",
            options = list(
                dom = 'Bfrtip', buttons = c('pageLength', 'copy', 'csv', 'excel'),
                stateSave = TRUE,
                order = isolate(input$reviewtable_state$order),
                paging = TRUE,
                pageLength = isolate(input$reviewtable_state$length),
                lengthMenu = list(c(10, 25, 100, -1), c('10', '25', '100','All')),
                columnDefs=list(list(orderable=TRUE))
            )   
        )
    }
    DT::renderDataTable({
        dtt %>% DT::formatStyle(
            'decision_str',
            target = 'row',
            backgroundColor = DT::styleEqual(
                c('Release', 'More Refinement', 'More Experiments', 'Reject'),
                c('#648FFF', '#FFB000',         '#FE6100',          '#DC267F')
            )
        ) %>% DT::formatStyle(
            'out_of_date',
            target = 'row',
            backgroundColor = DT::styleEqual(
                c('true', TRUE, 'TRUE'), c('#FFFFFF', '#FFFFFF', '#FFFFFF')
            )
        ) %>% DT::formatStyle(columns = 1:ncol(r1()),"white-space"="nowrap")
    }, server=FALSE
    )
}

#' Prebuffer a data table proxy
#'
#' @param proxy Proxy of datatable.
#' @param input input R6 object
#' @return Takes users state of data-table, refreshes the data and returns the state.
prebuffer_proxy <- function(proxy, input) { 
    message('(prebuffer_proxy)')
    updateSearch(proxy, keywords = list(global = input$aqp_state$search$search, columns = NULL)) # see input$therow_state$columns if needed
    selectPage(proxy, page = input$aqp_state$start/input$aqp_state$length+1)
}

#' Similar to updateMainTable functions but for the atom table!
#'
#' @param r Reactive data.frame of Ligand table.
#' @param pl Initial Pagelength
#' @param input input R6 object
#' @return data-table to be rendered
updateAQPTable <- function(r, pl=100, input=input){
    message('(updateAQPTable)')
    if(is.null(isolate(input$aqp_state))){
        dtt <- DT::datatable(
            r()[,c(1,39,40,41)],
            selection = 'single',
            options = list(
                stateSave=TRUE, 
                pageLength=pl, 
                columnDefs=list(list(orderable=TRUE, targets=c(0,1,2,3)))
            )
        )
    } else {
        dtt <- DT::datatable(
            r()[,c(1,39,40,41)],
            selection = 'single',
            options = list(
                stateSave = TRUE,
                order = isolate(input$aqp_state$order),
                paging = TRUE,
                pageLength = isolate(input$aqp_state$length), 
                columnDefs=list(list(orderable=TRUE, targets=c(0,1,2,3)))
            )
        )
    }
    DT::renderDataTable({
        dtt %>% DT::formatStyle(
            'bad_atom_index',
            target = 'row',
            backgroundColor = DT::styleEqual(
                c(''),
                c('#648FFF')
            )
        )
    }, server=FALSE)
}

#' Create view of data based on inputs selected...
#'
#' @param inputData Reactive data.frame of review data.
#' @param input R6 Object, specifying the state of UI elements.
#' @return A reactive data obtain, that responds to user interactions.
reactiviseData <- function(inputData, input){
    message('(ReactiviseData)')
    reactive({
        if(input$tab == 'review'){
            rowidx <- rep(FALSE, nrow(inputData())) # This broken?
            outcome <- as.numeric(as.character(inputData()$outcome))
            #rowidx <- rep(FALSE, length(outcome))
            current_state <- as.character(inputData()[['Current State']])
            review <- inputData()$decision_str
            if(any(c(is.null(input$protein), is.null(input$out4), is.null(input$out5), is.null(input$out6)))){
                inputData()[,]
            } else if(input$protein == '') {
                if(input$out4) rowidx[outcome == 4 | current_state == 'Comp Chem Ready'] <- TRUE
                if(input$out5) rowidx[outcome == 5 | current_state == 'Deposition Ready'] <- TRUE
                if(input$out6) rowidx[outcome == 6 | current_state == 'Deposited'] <- TRUE
                inputData()[rowidx,]
            } else {
                if(input$out4) rowidx[outcome == 4 | current_state == 'Comp Chem Ready'] <- TRUE
                if(input$out5) rowidx[outcome == 5 | current_state == 'Deposition Ready'] <- TRUE
                if(input$out6) rowidx[outcome == 6 | current_state == 'Deposited'] <- TRUE
                inputData()[rowidx & grepl(input$protein, as.character(inputData()$target_name)),]
            }
        } else if(input$tab == 'aqz'){
            if(is.null(input$aq_protein)){
                inputData()[,]
            } else {
                inputData()[grepl(input$aq_protein, as.character(inputData()$target_name)),]
            }
        }
    })
}

#' Restart R Session but keep users settings
#'
#' @param configuration List containing various config options including credentials.
#' @param target_list List targets they should see
#' @return reactive data of review table
restartSessionKeepOptions <- function(configuration, target_list){
    message('(restartSessionKeepOptions)')
    dbdat <- trygetReviewData(configuration = configuration, target_list=target_list)
    inputData <- reactive({dbdat})
    return(inputData)
}

#' Reset Review Form
#'
#' @param configuration List containing various config options including credentials.
#' @param session A shiny session object - which is generated when a user connects to XCR.
#' @param r1 Previous data... oh no...
#' @param possDec List of possible decisions...
#' @param target_list List targets they should see
#' @return Resets the review form to accept new responses.
resetForm <- function(configuration, session, r1, possDec, target_list){
    message('(resetForm)')
    updateSelectizeInput(session, "ligand", selected = '', choices = sort(rownames( r1() )))
    updateSelectInput(session, 'decision', selected ='', choices = possDec)
    updateSelectInput(session, 'reason', selected='', choices='')
    updateTextInput(session, 'comments', value = "")
    return(restartSessionKeepOptions(configuration=configuration, target_list=target_list))
}


#' Send an email using STFC SMP server
#'
#' @param structure String, Name of the structure
#' @param user String, The fedid of the user sending the email
#' @param decision String, the decision provided
#' @param reason String, The reason provided
#' @param comments String, Any additional Comments
#' @param email_list List object, containing email addresses as described in ./config.R
#' @return Does not return anything, will send an email to whoever is listed in email_list for the target the structure belongs to.
sendEmail <- function(structure, user, decision, reason, comments, email_list){
    message('(sendEmail)')
    target <- gsub('-[a-zA-Z]*[0-9]+_[0-9]*[A-Z]*', '', structure)

    # If a email list has not been specified for the target, use the defaults.
    if(is.null(email_list[[target]])){
        to <- sort(unique(email_list[['defaultUsers']]))
    } else {
        to <- sort(unique(email_list[[target]]))
    }

    sendmailR::sendmail(
        from = '<XChemStructureReview@diamond.ac.uk>',
        to = to,
        subject = sprintf('%s has been labelled as %s', structure, decision),
        msg = sprintf(
            '%s has been labelled as %s by %s for the following reason(s): %s.
With these additional comments:
%s
-------------------------------
If you wish to review this change please go to xchemreview.diamond.ac.uk while
connected to the diamond VPN or via NX.
This email was automatically sent by The XChem Review app
If you believe you have been sent this message in error, please email tyler.gorrie-stone@diamond.ac.uk',
            structure, decision, user, reason, comments, structure, target
        ),
        control = list(
            smtpServer = 'exchsmtp.stfc.ac.uk',
            smtpPort = 25
        )
    )
    modalDialog(
        title = 'Submission Sent', 
        'Your review has been sucessfully recorded, please select another structure!', 
        footer = modalButton("Dismiss"),
        size = 's', 
        easyClose = TRUE, 
        fade = TRUE
    )
}

#' Connect to, and then append to a db table, safely.
#'
#' @param configuration List containing configuration
#' @param table_name Name of the table to append to
#' @param data single row dataframe to append to the table.
#' @return Does not return anything.
dbAppendWrapper <- function(configuration, table_name, data){
    message('(dbAppendWrapper)')
    con <- dbConnect(RMariaDB::MariaDB(), dbname = configuration$db, host=configuration$host, 
        port=configuration$port, user=configuration$user, password=configuration$password)
    on.exit(dbDisconnect(con))
    dbAppendTable(conn = con, name = table_name, value = data, row.names=NULL)
}

#' Save a Review to XCDB...
#'
#' @param data Single row data.frame containing the row to append to the review table.
#' @param xtaln Name of the ligand.
#' @param email_list List containing emails of people who may need to be emails.
#' @param configuration List containing various config options
#' @return Does not return anything. Will append response to review table for the ligand and then send an email to whoever is listed in email_list for the target the structure belongs to.
saveData <- function(data, xtaln, email_list, configuration) {
    message('(saveData)')
    dbAppendWrapper(configuration = configuration, table_name = 'review_responses', data = data)
    sendEmail(xtaln, data[,'fedid'], data[,'decision_str'], data[,'reason'], data[,'comment'], email_list)
}

#' Get the epochTime
#'
#' @return An Integer representing the time
epochTime <- function() as.integer(Sys.time())

#' Get the time in a human readable format
#'
#' @return A string with the time YYYY-MM-DD-HH-MM-SS
humanTime <- function() format(Sys.time(), "%Y%m%d%H%M%OS")

# Don't ask about this kthx. It's yanked from SO
customDraggableModalDialog <- function(..., title = NULL,
                                 footer = shiny::modalButton("Dismiss"),
                                 size = c("m", "s", "l"),
                                 easyClose = FALSE, fade = FALSE) {
  size <- match.arg(size)
  cls <- if (fade) { "modal fade" } else { "modal" }
  shiny::div(
    id = "shiny-modal",
    class = cls,
    # tabindex = "-1", This line should be commented out or removed
    #`data-backdrop` = if (!easyClose) { "static" } ,
    #`data-keyboard` = if (!easyClose) { "false" } ,
    shiny::div(
      class = "modal-dialog",
      class = switch(size, s = "modal-sm", m = NULL, l = "modal-lg"),
      jqui_draggable(shiny::div(
        class = "modal-content",
        if (!is.null(title)) {
          shiny::div(
            class = "modal-header",
            shiny::tags$h4(class = "modal-title",  title)
          )
        },
        shiny::div(class = "modal-body", ...),
        if (!is.null(footer)) {
          shiny::div(class = "modal-footer", footer)
        }
      ))
    ),
    shiny::tags$script("$('#shiny-modal').modal().focus();"),
    shiny::tags$style(HTML("
    .modal-backdrop{
      display: none;
    }
    .modal {
      pointer-events: none;
    }
    .modal-content {
      pointer-events: all;
    }"))
  )
}
