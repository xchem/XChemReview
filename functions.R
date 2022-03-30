# List of Functions used in XCR...

#' Read a text file and return it as a single line separated by `\n` characters
#'
#' @param file A string corresponding to a filepath. 
#' @return Returns the contents of a file as a single string. 
readTxtToOneLine <- function(file){
    return(paste0(readLines(file), collapse='\n'))
}

#' Connect to XCDB using provided configuration
#'
#' @param configuration List containing various config options including xcdb credentials.
#' @return Returns a connection object that can be used to make queries. 
xcdbConnect <- function(configuration){
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

#' Fetch targets an authenticated fedid is able to see on ISPyB and XCDB.
#'
#' @param fedid String, corresponding to a fedid provided by shinyproxy/active directory when you log in.
#' @param configuration List containing various config options including xcdb credentials.
#' @return Vector of Targets the user is able to see. 
fetchTargets <- function(fedid, configuration){
    con <- xcdbConnect(configuration=configuration)
    on.exit(con)
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

#' Update datatable proxy.
#'
#' @param pdb a
#' @return Updates proxy with new data...
getResNum <- function(pdb){
    pdb <- bio3d::read.pdb(pdb)
    #return(as.character(pdb$atom$resno[1]))
    return(sprintf('@%s',paste(rownames(pdb$atom), collapse=',')))
}

#' Update datatable proxy.
#'
#' @param mol_file a
#' @param id_str a
#' @param comment_str a
#' @param name_str a
#' @return Updates proxy with new data...
write_to_mol_file <- function(mol_file, id_str, comment_str, name_str){
    lines <- readLines(mol_file)
    if(any(lines == '> <BADATOMS>')){
        # Don't update
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
    # And bad atom names...
    print(mol_file)
    print(lines)
    cat(paste(lines, collapse='\n'), file = mol_file)
}

#' Update datatable proxy.
#'
#' @param target a
#' @param configuration a
#' @return Updates proxy with new data...
rewriteMols <- function(target, configuration){ # This need refactoring later......
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))
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

#' Update datatable proxy.
#'
#' @param string a
#' @param split_by a
#' @return Updates proxy with new data...
rsplit <- function(string, split_by, n=1){
    spl <- strsplit(string, split_by)[[1]]
    c(paste(unlist(head(spl, length(spl)-n)), collapse=split_by), unlist(tail(spl, n)))
}

#' Update datatable proxy.
#'
#' @param pdb_file a
#' @return Updates proxy with new data...
get_residues <- function(pdb_file){
    struc <- try(bio3d::read.pdb(pdb_file), silent=T)
    if(inherits(struc, 'try-error')) return('')
    return(c('', unique(paste(struc$atom$resid, struc$atom$resno, sep='_'))))
}

#' Update datatable proxy.
#'
#' @param configuration a
#' @param target_list a
#' @return Updates proxy with new data...
trygetReviewData <- function(configuration, target_list){
   out <- try(getReviewData(configuration=configuration, target_list=target_list), silent=T)
   return(out)
}

#' Update datatable proxy.
#'
#' @param configuration a
#' @param target_list a
#' @return Updates proxy with new data...
getReviewData <- function(configuration, target_list){ 
    # Fetch Atom Q at the end??
    # Get data from target_list only first...
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
    return(output)
}

#' Update datatable proxy.
#'
#' @param configuration a
#' @param target_list a
#' @return Updates proxy with new data...
getFragalysisViewData <- function(configuration, target_list){

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
    return(output)
}

#' Update datatable proxy.
#'
#' @param configuration a
#' @param target a
#' @return Updates proxy with new data...
createUniqueMetaData <- function(configuration, target){
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))

    targs <- dbGetQuery(con, "SELECT * from fragalysis_target") 
    target_ids <- paste(targs$id[targs$target %in% target], collapse=',')

    fvdat <- dbGetQuery(con, sprintf("SELECT * from fragalysis_ligand WHERE fragalysis_target_id IN (%s)", target_ids))
    md_ids <- paste(fvdat$id, collapse=',')

    md <- dbGetQuery(con, sprintf("SELECT * FROM meta_data WHERE ligand_name_id in (%s)", md_ids))
    prot <- target
    md <- md[!(md$site_Label == 'IGNORE' | is.na(md$site_Label)),] # Might be a problem for later...
    meta <- md[grep(paste0(prot, '-'), md$fragalysis_name),c(6,7,3,3,4,2,5)]
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
    return(meta)
}


#' Update datatable proxy.
#'
#' @param ligand_name_id a
#' @param fragalysis_name a
#' @param original_name a
#' @param site_label a
#' @param new_smiles a
#' @param alternate_name a
#' @param pdb_id a
#' @param status a
#' @param configuration a
#' @return Updates proxy with new data...
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
        dbAppendWrapper(configuration = configuration, table_name = 'meta_data', data = df)
    } else { # Update
        con <- xcdbConnect(configuration=configuration)
        dbExecute(con, sprintf("UPDATE meta_data SET %s WHERE ligand_name_id=%s",
            sprintf("site_Label=\'%s\', new_smiles=\'%s\', alternate_name=\'%s\', pdb_id=\'%s\', fragalysis_name=\'%s\', original_name=\'%s\', status=\'%s\'", site_label, new_smiles, alternate_name, pdb_id, fragalysis_name, original_name, status),
            ligand_name_id)
        )
        dbDisconnect(con)
    }
}

#' Update datatable proxy.
#'
#' @param sID a
#' @param text a
#' @return Updates proxy with new data...
debugMessage <- function(sID, text){
    message(sprintf('sid: %s | %s | %s', sID, text, Sys.time()))
}

#' Update datatable proxy.
#'
#' @param values a
#' @param title a
#' @return Updates proxy with new data...
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

#' Update datatable proxy.
#'
#' @param raster a
#' @param rank a
#' @param colour a
#' @return Updates proxy with new data...
changeRasterRank <- function(raster, rank, colour){
    raster[floor(rank*100)] <-  colour
    return(raster)
}

#' Update datatable proxy.
#'
#' @param data a
#' @param title a
#' @param target_name a
#' @return Updates proxy with new data...
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


#' Update datatable proxy.
#'
#' @param meta a
#' @param target a
#' @param copymaps a
#' @param mtz a
#' @param configuration a
#' @return Updates proxy with new data...
createFragUploadFolder <- function(meta, target, copymaps=FALSE, mtz, configuration){
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Creating Fragalysis Folder", value = 0)
    prot = target
    protsuffix <- paste(prot, format(Sys.time(), "%Y%m%d_%H%M"), sep='_')
    base_root <- sprintf('%s/%s/',configuration$staging, prot)
    extrafiles1 <- file.path(base_root, 'extra_files')
    rootf <- sprintf('%s/%s',configuration$prep_path, protsuffix)
    basef <- sprintf('%s/%s/%s', configuration$prep_path,protsuffix, prot)
    align_dir <- sprintf('%s/%s/%s/aligned',configuration$prep_path, protsuffix, prot)
    crys_dir <- sprintf('%s/%s/%s/crystallographic',configuration$prep_path, protsuffix, prot)

    dir.create(rootf)
    dir.create(basef)
    dir.create(align_dir)
    dir.create(crys_dir)

    # Maybe use file.copy(...,rec=T)?
    try(system(sprintf('cp -r %s %s', extrafiles, basef)), silent = TRUE)

    # aligned data copy
    progress$set(message = "Copying aligned Files", value = 0)
    increment = (1/nrow(meta))/2.2
    print(meta)
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
    zip(zipfile=zipf, files=prot)
    cat('', file='done.log')
    setwd(currentwd)
    system(sprintf('chmod 775 %s', rootf))
    system(sprintf('chmod 775 %s/%s', rootf, prot))
    system(sprintf('chmod 775 %s/%s', rootf, zipf))
    full_path_zipf <- sprintf('%s/%s', rootf, zipf)
    return(full_path_zipf)
}


#' Update datatable proxy.
#'
#' @param filepath a
#' @param target a
#' @param proposal a
#' @param email a
#' @param configuration a
#' @return Updates proxy with new data...
uploadFragFolder <-  function(filepath, target, proposal, email, configuration){
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Uploading Fragalysis Folder", value = 0)
    external_script <- file.path(configuration$script_path, 'upload_2_fragalysis.sh')
    command <- sprintf('%s %s %s %s %s %s', external_script, basename(filepath), target, proposal, email, filepath)
    task <- tail(system(command, intern=T),1)
    return(task)
}

#' Update datatable proxy.
#'
#' @param x a
#' @return Updates proxy with new data...
perc.rank <- function(x) trunc(rank(x))/length(x)

#' Update datatable proxy.
#'
#' @param x a
#' @param xo a
#' @return Updates proxy with new data...
perc.rank2 <- function(x, xo=NULL)  length(x[x <= xo])/length(x)

#' Update datatable proxy.
#'
#' @param x a
#' @param alldata a
#' @param extradata a
#' @return Updates proxy with new data...   
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


#' Update datatable proxy.
#'
#' @param name a
#' @param boxsize a
#' @param session a
#' @return Updates proxy with new data...   
updateDensityBoxSize <- function(name, boxsize, session) session$sendCustomMessage('updateVolumeDensityBoxSize', list(name, boxsize))

#' Update datatable proxy.
#'
#' @param name a
#' @param isolevel a
#' @param session a
#' @return Updates proxy with new data...      
updateDensityISO <- function(name, isolevel, session) session$sendCustomMessage('updateVolumeDensityISO', list(name, isolevel))
        

#' Update datatable proxy.
#'
#' @param name a
#' @param bool a
#' @param session a
#' @return Updates proxy with new data...        
updateVisability <- function(name, bool, session){
    session$sendCustomMessage(
        type = 'updateVolumeDensityVisability',
        list(
            as.character(name),
            tcl(bool)
        )
    )
}

#' Update datatable proxy.
#'
#' @param filepath a
#' @param color a
#' @param negateiso a
#' @param boxsize a
#' @param visable a
#' @param windowname a
#' @param isotype a
#' @param session a
#' @return Updates proxy with new data...
uploadVolumeDensity <- function(filepath, color, negateiso = FALSE, boxsize, isolevel, visable, windowname, isotype='sigma', session){
    volume_bin <- readBin(filepath, what='raw', file.info(filepath)$size)
    volume_b64 <- base64encode(volume_bin, size=NA, endian=.Platform$endian)
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

#' Update datatable proxy.
#'
#' @param x Proxy of datatable.
#' @return Updates proxy with new data...
getExt <- function(x) sapply(strsplit(x, '[.]'), tail, 1)

#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data...
uploadUnfocussedMol <- function(filepath, session){
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type='addMol',
        list(choice)
    )
}

#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param ext Proxy of datatable.
#' @param focus Proxy of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data...
uploadMolAndFocus <- function(filepath, ext, focus, session){
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'addMolandfocus',
        list(choice, ext, tcl(focus))
    )
}

#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param clear Proxy of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data...
uploadBFactors <- function(filepath, clear=TRUE, session){
    if(clear) clearWindowField(id='bfactor', session=session)
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'setBFactor',
        message = list(
            choice
        )
    )
}


#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data..
addContacts <- function(filepath, session){
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'addContacts',
        message = list(
            choice
        )
    )
}

#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param ext Proxy of datatable.
#' @param focus of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data..
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
    session$sendCustomMessage(
        type = 'setBadMolandfocus',
        list(choice,ext, tcl(focus), paste(mid, collapse=';'), paste(mic,collapse=';'))
    )
}

#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param repr Proxy of datatable.
#' @param focus of datatable.
#' @param molfile of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data..
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

#' Update datatable proxy.
#'
#' @param x Proxy of datatable.
#' @return Updates proxy with new data.
tcl <- function(x) tolower(as.character(as.logical(x)))

#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param repr Proxy of datatable.
#' @param focus of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data..
uploadApoPDB <- function(filepath, repr, focus, session){
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'setapoPDB',
        message = list(
            choice,
            repr,
            tcl(focus)
        )
    )
}

#' Update datatable proxy.
#'
#' @param objectname Proxy of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data..
removeNamedComponent <- function(objectname, session) session$sendCustomMessage(type='removeNamedComponent', list(objectname))

#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param input Proxy of datatable.
#' @param focus_point of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data..
uploadPDB <- function(filepath, input, focus_point = 'LIG', session){
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'setPDB2', # See TJGorrie/NGLShiny for details on setPDB2
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

#' Update datatable proxy.
#'
#' @param which Proxy of datatable.
#' @param what Proxy of datatable.
#' @param session of datatable.
#' @return Updates proxy with new data..
updateParam <- function(which, what, session) session$sendCustomMessage('updateaparam', list(which, what))

#' Update datatable proxy.
#'
#' @return Updates proxy with new data..
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

#' Update datatable proxy.
#'
#' @return Updates proxy with new data..
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

#' Update datatable proxy.
#'
#' @param id Proxy of datatable.
#' @param configuration Proxy of datatable.
#' @param atomstoquery Proxy of datatable.
#' @return Updates proxy with new data..
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
        dbAppendWrapper(configuration = configuration, table_name = 'bad_atoms', data = newdat)
    } else {
        con <- xcdbConnect(configuration=configuration)
        dbExecute(con, sprintf(
            "UPDATE bad_atoms SET %s WHERE ligand_id=%s",
            sprintf("atomid=\'%s\', comment=\'%s\', atomname=\'%s\'", newdat$atomid, newdat$comment, newdat$atomname),
            ligand_id)
        )
        dbDisconnect(con)        
    }
    # Write atom data to .mol file???
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
    }
}

#' Update datatable proxy.
#'
#' @param id Proxy of datatable.
#' @param configuration Proxy of datatable.
#' @return Updates proxy with new data..
displayModalWhoUpdated <- function(id, configuration){
    con <- xcdbConnect(configuration=configuration)
    on.exit(dbDisconnect(con))
    response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses"))
    mostrecent <- as.data.frame(t(sapply(split(response_data, response_data$ligand_name_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
    rownames(mostrecent) <- as.character(mostrecent$ligand_name_id)
    user <- mostrecent[as.character(id), 'fedid'][[1]]
    showModal(
        modalDialog(
            title = "Someone has recently reviewed this crystal",
            sprintf("A User (%s) has recently reviewed this structure. Restarting the session to update their response. If you disagree with the current response, please submit another response or select another crystal.", user),
            footer = tagList(actionButton("ok", "Okay"))
        )
    )
}

#' Update datatable proxy.
#'
#' @param id Proxy of datatable.
#' @param sessionTime Proxy of datatable.
#' @param configuration Proxy of datatable.
#' @return Updates proxy with new data...
sessionGreaterThanMostRecentResponse <- function(id, sessionTime, configuration){
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
    return(TRUE) # Just force it... The app updates fairly rapidly now...
}


#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param color Proxy of datatable.
#' @param session Proxy of datatable.
#' @return Updates proxy with new data...
uploadMolNoFocus <- function(filepath, color, session){
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'fv_addMolandfocus_withcolor',
        list(choice, color)
    )
}

#' Update datatable proxy.
#'
#' @param id Proxy of datatable.
#' @param session Proxy of datatable.
#' @return Updates proxy with new data...
clearWindowField <- function(id, session){
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

#' Update datatable proxy.
#'
#' @param filepath Proxy of datatable.
#' @param session Proxy of datatable.
#' @return Updates proxy with new data...
uploadMolAndFocus2 <- function(filepath, session){
    # Used once...
    choice <- readTxtToOneLine(file=filepath)
    session$sendCustomMessage(
        type = 'fv_addMolandfocus',
        list(choice)
    )
}

#' Update datatable proxy.
#'
#' @param configuration Proxy of datatable.
#' @param target_list Proxy of datatable.
#' @return Updates proxy with new data...
reactivegetFragalysisViewData <- function(configuration, target_list){
    reactive({getFragalysisViewData(configuration=configuration, target_list=target_list)})
}

#' Update datatable proxy.
#'
#' @param flexdata Proxy of datatable.
#' @param input Proxy of datatable.
#' @return Updates proxy with new data...
react_fv_data2 <- function(data, input){
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

#' Update datatable proxy.
#'
#' @param flexdata Proxy of datatable.
#' @param input Proxy of datatable.
#' @return Updates proxy with new data...
react_fv_data <- function(data, input){
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

#' Update datatable proxy.
#'
#' @param flexdata Proxy of datatable.
#' @return Updates proxy with new data...
updateFlexPlot <- function(flexdata){
    renderPlotly({
        plot_ly(flexdata()$data, x=~x, y=~y, text=~ligand_name, color=~status, customdata = ~ligand_name, size=20) %>% config(scrollZoom = TRUE)
    })
}

#' Update datatable proxy.
#'
#' @param r1 Proxy of datatable.
#' @param input input R6 object
#' @return Updates proxy with new data...
flexPlotDataFun <- function(r1, input){
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

#' Update datatable proxy.
#'
#' @param r1 Proxy of datatable.
#' @param pl Proxy of datatable.
#' @param input input R6 object
#' @return Updates proxy with new data...
updateMainTable2 <- function(r1, pl=100, input){
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

#' Update datatable proxy.
#'
#' @param r1 Proxy of datatable.
#' @param pl Proxy of datatable.
#' @param input input R6 object
#' @return Updates proxy with new data...
updateMainTable <- function(r1, pl=100, input){
    if (is.null(isolate(input$reviewtable_state))) {
        dtt <- DT::datatable(
            r1(),
            selection = 'single',
            callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
            options = list(
                stateSave = TRUE, 
                pageLength=pl, 
                columnDefs=list(list(orderable=TRUE))
            )
        )
    } else {
        dtt <- DT::datatable(
            r1(),
            selection = 'single',
            callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
            options = list(
                stateSave = TRUE,
                order = isolate(input$reviewtable_state$order),
                paging = TRUE,
                pageLength = isolate(input$reviewtable_state$length), 
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
    }, server=FALSE)
}

#' Update datatable proxy.
#'
#' @param proxy Proxy of datatable.
#' @param input input R6 object
#' @return Updates proxy with new data...
prebuffer_proxy <- function(proxy, input) { 
    updateSearch(proxy, keywords = list(global = input$aqp_state$search$search, columns = NULL)) # see input$therow_state$columns if needed
    selectPage(proxy, page = input$aqp_state$start/input$aqp_state$length+1)
}

#' Create view of data based on inputs selected...
#'
#' @param r Reactive data.frame of Ligand table.
#' @param pl Initial Pagelength
#' @param input input R6 object
#' @return A reactive data obtain, that responds to user interactions.
updateAQPTable <- function(r, pl=100, input=input){
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
    reactive({
        if(input$tab == 'review'){
            rowidx <- rep(FALSE, nrow(inputData()))
            outcome <- as.numeric(as.character(inputData()$outcome))
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

#' Fetch review data to 'update' the application
#'
#' @param configuration List containing various config options including credentials.
#' @return Returns a data.frame for the data to be reviewed.
restartSessionKeepOptions <- function(configuration, target_list){
    dbdat <- trygetReviewData(configuration = configuration, target_list=target_list)
    inputData <- reactive({dbdat})
    return(inputData)
}

#' Reset Review Form
#'
#' @param configuration List containing various config options including credentials.
#' @param session Session R6 object. Probably not needed but included anyway... Can be obtained from higher scope...
#' @param r1 Previous data... oh no...
#' @param possDec List of possible decisions...
#' @return Resets the review form to accept new responses.
resetForm <- function(configuration, session, r1, possDec){
    updateSelectizeInput(session, "ligand", selected = '', choices = sort(rownames( r1() )))
    updateSelectInput(session, 'decision', selected ='', choices = possDec)
    updateSelectInput(session, 'reason', selected='', choices='')
    updateTextInput(session, 'comments', value = "")
    return(restartSessionKeepOptions(configuration=configuration))
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

# This should be functions.R but I am not sure that the code will behave correctly if I 
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