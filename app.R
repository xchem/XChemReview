local = FALSE
# Generic Shiny Libraries
library(httr)
library(shiny)
library(grDevices)
library(lubridate)
library(shinydashboard)
library(shinyjqui)
library(shinyWidgets)
library(DT)
# DB Lib
library(DBI)
library(sendmailR)
# Read Binaries to encode to b64
library(caTools)
# Read Pdb files
library(bio3d)
# Render structures in ngl window
if(local) {
    library(nglShiny)
    source('./my_config.R')
} else {
    source('/dls/science/groups/i04-1/software/xchemreview/config.R')
    #install.packages('/dls/science/groups/i04-1/software/nglshiny', repos=NULL, type='source', lib='/dls/science/groups/i04-1/software/xchemreview/xcrlib')
    library(nglShiny, lib.loc='/dls/science/groups/i04-1/software/xchemreview/xcrlib')
}

#fragfolders <- c('', 'Mpro', 'PlPro', 'PHIPA', 'NSP16','PGN_RS02895PGA', 'XX02KALRNA', 'CD44MMA')
#target_list <- c('PGN_RS02895PGA','Mpro', 'PlPro', 'PHIPA', 'NSP16', 'XX02KALRNA', 'CD44MMA')
target_list <- sort(c(
	#'CD44MMA'
	'Mpro'
	#'PlPro',
	#'NSP16',
	#'PGN_RS02895PGA',
	#'PHIPA'
))
fragfolders <- c('', target_list)
# Plotting Libs
library(ggplot2)
library(plotly)

sessionInfo()
#fragfolders <- c('', 'Mpro', 'PlPro', 'PHIPA', 'NSP16','PGN_RS02895PGA', 'XX02KALRNA')

# Update existing mol_files with latest atom comments
con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
atoms <- dbGetQuery(con, 'SELECT * from "BadAtoms"')
reviews <- dbGetQuery(con, 'SELECT * from "review_responses_new"')
fligands <- dbGetQuery(con, 'SELECT * from "FragalysisLigand"')
fligand_id <- fligands$id
rownames(fligands) <- as.character(fligand_id)
ligands <- dbGetQuery(con, 'SELECT * from "ligand"')
dbDisconnect(con)
ligand_id <- ligands$id
rownames(ligands) <- as.character(ligand_id)
latest_review <-  t(sapply(split(reviews, reviews$Ligand_name_id), function(x) x[which.max(x$time_submitted),]))



write_to_mol_file <- function(mol_file, id_str, comment_str){
    lines <- readLines(mol_file)
    if(any(lines == '> <BADATOMS>')){
        # Don't update
        #badid_line <- which(lines == '> <BADATOMS>') + 1
        #lines[badid_line] <- id_str
    } else {
        lines <- c(lines, '> <BADATOMS>', id_str)
    }
    if(any(lines == '> <BADCOMMENTS>')){
        #badcomment_line <- which(lines == '> <BADCOMMENTS>') + 1
        #        lines[badcomment_line] <- badcommentstr
    } else {
        lines <- c(lines, '> <BADCOMMENTS>', comment_str)
    }
    # And bad atom names...
    cat(paste(lines, collapse='\n'), file = mol_file)
}

for(i in 1:nrow(latest_review)){
    lig <- latest_review[i,'Ligand_name_id']
    rev <- latest_review[i,'id']
    mol_file <- fligands[as.character(ligands[as.character(lig),'fragalysis_ligand_id']),]$lig_mol_file
    ids_to_grab <- which(atoms$Review_id==rev)
    if(length(ids_to_grab)>0){
	print(mol_file)
        id_str <- atoms[ids_to_grab, 'atomid']
        comment_str <- atoms[ids_to_grab, 'comment']
        id_str2 <- paste0(id_str, collapse=';')
        comment_str2 <- paste0(comment_str, collapse=';')
        write_to_mol_file(mol_file=mol_file, id_str=id_str2, comment_str = comment_str2)
    }
}


# I can't believe this doesn't exist in R!
# Functional equivalent to string.rsplit('_', 1)
rsplit <- function(string, split_by, n=1){
    spl <- strsplit(string, split_by)[[1]]
    c(paste(unlist(head(spl, length(spl)-n)), collapse=split_by), unlist(tail(spl, n)))
}

get_residues <- function(pdb_file){
    struc <- try(bio3d::read.pdb(pdb_file), silent=T)
    if(inherits(struc, 'try-error')) return('')
    return(c('', unique(paste(struc$atom$resid, struc$atom$resno, sep='_'))))
}

getReviewData <- function(db, host_db, db_port, db_user, db_password, target_list){
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    ligand_data <- dbGetQuery(con, "SELECT * from ligand")
    lig_crys_ids <- paste(ligand_data[,'crystal_id'], collapse=',')
    ligand_crystal_data <- dbGetQuery(con, sprintf("SELECT id, crystal_name, compound_id, target_id FROM crystal WHERE id IN (%s)", lig_crys_ids))
    rownames(ligand_crystal_data) <- as.character(ligand_crystal_data$id)
    lig_fl_ids <- paste(ligand_data[,'fragalysis_ligand_id'], collapse=',')
    ligand_fl_data <- dbGetQuery(con, sprintf("SELECT * from \"FragalysisLigand\" WHERE id IN (%s)", lig_fl_ids))
    rownames(ligand_fl_data) <-  as.character(ligand_fl_data$id)
    ligand_refinement_data <- dbGetQuery(con, sprintf("SELECT id, crystal_name_id, r_free, rcryst, ramachandran_outliers, res, rmsd_angles, rmsd_bonds, lig_confidence_string, spacegroup, outcome, cif, pdb_latest, mtz_latest FROM refinement WHERE crystal_name_id IN (%s)", lig_crys_ids))
    rownames(ligand_refinement_data) <- as.character(ligand_refinement_data$crystal_name_id)

    ligand_target_data <- dbGetQuery(con, sprintf("SELECT * FROM target WHERE id IN (%s)", paste(ligand_crystal_data[,'target_id'], collapse=',')))
    rownames(ligand_target_data) <- as.character(ligand_target_data$id)
    ligand_compound_data <- dbGetQuery(con, sprintf("SELECT * FROM compounds WHERE id IN (%s)", paste(ligand_crystal_data[,'compound_id'], collapse=',')))
    rownames(ligand_compound_data) <- as.character(ligand_compound_data$id)
    ligand_response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses_new"))
    mostrecent <- as.data.frame(t(sapply(split(ligand_response_data, ligand_response_data$Ligand_name_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
    if(nrow(mostrecent)>1){
        rownames(mostrecent) <- as.character(mostrecent$Ligand_name_id)
    } else {
        mostrecent <- ligand_response_data
        rownames(mostrecent) <- as.character(mostrecent$Ligand_name_id)
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
        'SMILES' = ligand_compound_data[as.character(ligand_crystal_data[as.character(ligand_data$crystal_id),]$compound_id),2],
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
            ))
    )
    rownames(output) <- make.names(as.character(output$ligand_name), unique=TRUE)
    output <- output[output$target_name %in% target_list, ] # Add to list as more targets needed?
    dbDisconnect(con)

    return(output)
}

getFragalysisViewData <- function(db, host_db, db_port, db_user, db_password, target_list){
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    # Get All ligands that are reviewed as release or above. Mandatory...
    ligand_response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses_new"))
    mostrecent <- as.data.frame(t(sapply(split(ligand_response_data, ligand_response_data$Ligand_name_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
    to_release_ids <- unlist(mostrecent$Ligand_name_id[mostrecent$decision_int==1])
    liganded_ligands <- dbGetQuery(con, "SELECT fragalysis_ligand_id, id from ligand")
    ind <- as.character(liganded_ligands[,2]) %in% as.character(to_release_ids)
    fvdat <- dbGetQuery(con, "SELECT * from \"FragalysisLigand\"")
    md <- dbGetQuery(con, "SELECT * FROM \"MetaData\"")
    rownames(md) <- md$Ligand_name_id
    # Show all even if not reviewed...
    annotatable_fv_dat <- fvdat#fvdat[!fvdat$id %in% liganded_ligands[!ind, 1], ]
    targets <- dbGetQuery(con, sprintf("SELECT * from \"FragalysisTarget\" WHERE id IN (%s)", paste(unique(annotatable_fv_dat$fragalysis_target_id), collapse=',')))
    rownames(targets) <- as.character(targets$id)
    output <- cbind(annotatable_fv_dat, targetname=targets[as.character(annotatable_fv_dat$fragalysis_target_id), 'target'], md[as.character(annotatable_fv_dat$id),])
    dbDisconnect(con)
    rns <- gsub('[.]', '-', make.names(as.character(output$ligand_name), unique=TRUE))
    numbers <- grepl('^[0-9]', output$ligand_name)
    rns[numbers] <-  gsub('^X{1}', '', rns[numbers])
    rownames(output) <- rns
    output <- output[output$targetname %in% target_list, ]
    return(output)
}

createUniqueMetaData <- function(db, host_db, db_port, db_user, db_password, target){
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    fvdat <- dbGetQuery(con, "SELECT * from \"FragalysisLigand\"")
    md <- dbGetQuery(con, "SELECT * FROM \"MetaData\"")
    prot <- target
    md <- md[!(md$Site_Label == 'IGNORE' | is.na(md$Site_label)),]
    meta <- md[grep(paste0(prot, '-'), md$fragalysis_name),c(6,7,3,3,4,2,5)]
    colnames(meta) <- c('crystal_name', 'RealCrystalName', 'smiles', 'new_smiles', 'alternate_name', 'site_name', 'pdb_entry')

    # Dodgy Hack for DLS
    rootf = sprintf('/dls/science/groups/i04-1/fragprep/staging_test/%s/aligned/', prot)
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

createFragUploadFolder <- function(meta, target, copymaps=FALSE, mtz){
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Creating Fragalysis Folder", value = 0)
    prot = target
    base_root <- sprintf('/dls/science/groups/i04-1/fragprep/staging_test/%s/', prot)
    protsuffix <- paste(prot, format(Sys.time(), "%Y%m%d_%H%M"), sep='_')

    rootf <- sprintf('/dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/%s', protsuffix)
    basef <- sprintf('/dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/%s/%s', protsuffix, prot)
    align_dir <- sprintf('/dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/%s/%s/aligned', protsuffix, prot)
    crys_dir <- sprintf('/dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/%s/%s/crystallographic', protsuffix, prot)

    system(sprintf('mkdir /dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/%s', protsuffix))
    system(sprintf('mkdir /dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/%s/%s', protsuffix, prot))
    system(sprintf('mkdir /dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/%s/%s/aligned', protsuffix, prot))
    system(sprintf('mkdir /dls/science/groups/i04-1/fragprep/FragalysisUploadFolders/%s/%s/crystallographic', protsuffix, prot))

    # aligned data copy
    progress$set(message = "Copying aligned Files", value = 0)
    increment = (1/nrow(meta))/2.2
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
        file.copy(file.path(cf,files), file.path(nf, files))
        if(copymaps){
		# Replace uncut maps with mtz files...
		# IF RCN is in mtz
		mtz_file = mtz[,1] == rcn
                if(sum(mtz_file)>0){
                	id <- which(mtz_file)[1]
                        mtz_file <- mtz[id,2]
			file.copy(mtz_file, file.path(nf, sprintf('%s_refine.mtz', rcn)))
			# Find event.mtz
			print(mtz_file)
			print(dirname(mtz_file))
			print(dirname(dirname(mtz_file)))
			rootmtz <- dir(dirname(dirname(mtz_file)), pattern='.mtz')
			print(rootmtz)
			event_mtz <- rootmtz[grep('event',rootmtz)]
			print(event_mtz)
			if(length(event_mtz)>0){
				newnames <- basename(event_mtz)
				ends <- sapply(strsplit(newnames, '_'), tail, 1)
				print(ends)
				file.copy(file.path(dirname(dirname(mtz_file)), event_mtz), file.path(nf, sprintf('%s_event_%s', rcn, ends)))
			}
		}
	}
    }

    write.csv(meta, sprintf("%s/metadata.csv", basef), quote = FALSE)
    write.csv(meta, sprintf("%s/metadata.csv", align_dir), quote = FALSE)

    progress$set(message = "Zipping File!", value = .9)

    # Zip File
    zipf <- sprintf('%s.zip', prot)
    zipcommand <- sprintf("(cd %s && zip -r %s .)", rootf, prot)
    system(zipcommand)
    system(sprintf('chmod 775 %s', rootf))
    system(sprintf('chmod 775 %s/%s', rootf, prot))
    system(sprintf('chmod 775 %s/%s', rootf, zipf))

    full_path_zipf <- sprintf('%s/%s', rootf, zipf)

    return(full_path_zipf)
}


updateOrCreateRow <- function(ligand_name_id, fragalysis_name, original_name, site_label='', new_smiles='', alternate_name='', pdb_id='',
    dbname, host, port, user, password){
    df = data.frame(Site_Label=as.character(site_label),
                    new_smiles=as.character(new_smiles),
                    alternate_name=as.character(alternate_name),
                    pdb_id=as.character(pdb_id),
                    fragalysis_name=as.character(fragalysis_name),
                    original_name=as.character(original_name),
                    Ligand_name_id=as.character(ligand_name_id)
                    )
    print(df)
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    id = dbGetQuery(con, sprintf("SELECT id from \"MetaData\" WHERE \"Ligand_name_id\"=%s", df$Ligand_name_id))[1,1]
    if(is.na(id)){
        message('Creating MetaRow')
        dbAppendTable(con, "MetaData", value = df, row.names=NULL)
    } else {
        message('Updating MetaRow!')
        dbExecute(con, sprintf("UPDATE \"MetaData\" SET %s WHERE \"Ligand_name_id\"=%s",
            sprintf("\"Site_Label\"=\'%s\', new_smiles=\'%s\', alternate_name=\'%s\', pdb_id=\'%s\', fragalysis_name=\'%s\', original_name=\'%s\'", site_label, new_smiles, alternate_name, pdb_id, fragalysis_name, original_name),
            ligand_name_id)
        )
    }
    dbDisconnect(con)
}


getReviewRow <- function(data, db, host_db, db_port, db_user, db_password){
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    ids <- with(data, dbGetQuery(con, sprintf("SELECT \"Ligand_name_id\", id FROM review_responses_new WHERE fedid='%s' AND time_submitted='%s'", fedid, time_submitted)))[1,1:2]
    dbDisconnect(con)
    return(ids)
}

debugMessage <- function(sID, text){
    message(sprintf('sid: %s | %s | %s', sID, text, Sys.time()))
}

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

colfunc <- colorRampPalette(c("red", "white", "blue"))

changeRasterRank <- function(raster, rank, colour){
    raster[floor(rank*100)] <-  colour
    return(raster)
}

hmapbar <- function(data, title, target_name){
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


header <- dashboardHeader(
    title = 'XChemReview',
    tags$li(class='dropdown', actionButton('controls', 'Additional NGL Controls', class = 'btn-primary', icon = icon('cog', lib = 'glyphicon')))
)

sidebar <- dashboardSidebar(
    sidebarMenu(
        id = 'tab',
        menuItem('Summary', tabName = 'summary', icon=icon('th'), badgeLabel = 'new', badgeColor = 'green'),
        menuItem('FragView', tabName = 'fragview', icon = icon('dashboard')),
        menuItem('Review', tabName = 'review', icon = icon('dashboard')),
        menuItem('LaunchPad', tabName = 'launchpad', icon = icon('th'), badgeLabel = 'new', badgeColor = 'green'),
        menuItem('Help', tabName = 'help', icon = icon('th')),
        hr(),
        # Flexible Sidebar options depending on which menuitem is selected.
        uiOutput('flex1')
    )
)

# Copied from SO
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

body <- dashboardBody(
	tags$head(tags$script("$(function() {$.fn.dataTableExt.errMode = 'throw';});")),
    tags$head(shiny::tags$style(HTML("
    .modal-backdrop{
      display: none;
    }
    .modal {
      pointer-events: none;
    }
    .modal-content {
      pointer-events: all;
    }"))),
    tabItems(
        # First Tab
        tabItem(
            tabName = 'review',
            fluidRow( 
                nglShinyOutput('nglShiny', height = '500px'),
                jqui_draggable(
                    tabBox(
                        tabPanel(
                            title = 'NGL Controls',
                            fluidRow(
                                column(6, actionButton(
                                    "fitButton",
                                    "Center on Ligand"
                                )),
                                column(6,checkboxInput('autocenter', 'Automatically Center on load', value=TRUE))
                            ),
                            fluidRow(
                                shinyWidgets::chooseSliderSkin("Flat", color='#112446'),
                                column(6,
                                    fluidRow(
                                        column(2, checkboxInput('eventMap', 'Show Event Map', value = TRUE)),
                                        column(10, sliderInput("isoEvent", "", min = 0, max = 3, value = 1, step = 0.1))
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('twofofcMap', 'Show 2fofc Map', value = TRUE)),
                                        column(10, sliderInput("iso2fofc", "", min = 0, max = 3, value = 1.5, step = 0.1))
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('fofcMap', 'Show fofc Map', value = TRUE)),
                                        column(10, sliderInput("isofofc", "", min = 0, max = 3, value = 3, step = 0.1))
                                    ),
                            	    selectInput('gotores', 'Go to Residue:', choices = '', multiple=FALSE),
                            	    selectizeInput('highlight_res', 'Highlight Residues:', choices = '', multiple=TRUE)
                                ),
                                column(6,
                                    imageOutput('ligimage2', height='300px'),
                                    radioButtons('views', 'View Type', selected = 'aligned', inline = FALSE, width = NULL,
                                        choiceNames = c('Aligned (what will be in Fragalysis)', 'Unaligned (to check if the api alignment introduces problems)', 'Raw Input Files (What you should see in coot, maps may take long time to load)'),
                                        choiceValues = c('aligned', 'unaligned', 'crystallographic')
                                    ),
                                    selectInput('asuSwitch', 'Assembly Type (Only in Raw and Unalign)', selected='AU', choices=c('AU', 'UNITCELL', 'SUPERCELL'))
                                )
                            )
                        ),
                        tabPanel(
                            title = 'Ligand Information',
                            div(style='overflow-y:scroll;height:600px;',
                            fluidRow(
                                column(8,
                                    column(6, 
                                        helpText('Input Ligand'),
                                        imageOutput('ligimage')
                                    ),
                                    column(6, 
                                        helpText('Modelled Ligand'),
                                        imageOutput('rendered_ligimage')
                                    )
                                ),
                                column(4,
                                    div(style = "margin-top:-1em", checkboxInput('renderMisc', 'Render Image/Spider Plot', value = TRUE, width = NULL)),
                                    div(style = "margin-top:-1em", selectInput('emap', 'Select Eventmap', choices='', multiple=FALSE)),
                                    actionButton('buster', 'Buster Report'),
                                    #div(style = "margin-top:-1em", selectInput('scope', 'Scope', c('Experiment', 'Global'))),
                                    #div(style = "margin-top:-1em", selectInput('plotType', 'Statistic', c('res', 'r_free', 'rcryst', 'ramachandran_outliers', 'rmsd_angles', 'rmsd_bonds')))

                                )
                            ),
                            column(12,div(style = "margin-top:-15em",
                                fluidRow(
                                    uiOutput('plotElement')
                                )
                            ))
                            ),
                        ),
                        tabPanel(
                            title = 'Atom Selection (Alt + Left Click)',
                            textOutput('as_message'),
                            fluidRow(
                                column(3, actionButton('as_clear', label = 'Clear Atoms')),
                                column(3, fluidRow(
                                actionButton('write_all', 'Write to All Atoms?', value=FALSE),
                                actionButton('write_selected', label = 'Write to selected rows')
                                )),
                                column(3, selectizeInput('atom_text', 'Comment', choices=c('', 'Weak Density', 'No Density Evidence', 'Unexpected Atom', 'Multiple Conformations'), options=list(create=TRUE)))
                            ),
                            DT::dataTableOutput('atoms')
                        )
                    ), options = list(delay = '1000', cancel = '.selectize-control')
                ),
                jqui_draggable(
                    tabBox(
                        tabPanel(
                            title='Review Table',
                            div(
                                style='overflow-y:scroll;height:600px;',
                                DT::DTOutput('reviewtable') # Great Name!
                            )
                        ),
                        tabPanel(
                            title='Review Plots (Click points to load ligand)',
                            fluidRow(
                                column(4,
                                        selectInput('fpex', 'x', selected = 'res', choices=c('res', 'r_free', 'rcryst', 'ramachandran_outliers', 'rmsd_angles', 'rmsd_bonds'))
                                ),
                                column(4,
                                    selectInput('fpey', 'y', selected = 'r_free', choices=c('res', 'r_free', 'rcryst', 'ramachandran_outliers', 'rmsd_angles', 'rmsd_bonds'))
                                ),
                                column(4,
                                    selectInput('fpe_target', 'target', selected = '', choices=c('A', 'B', 'C'))
                                )
                            ),
                            verbatimTextOutput("info"),
                            uiOutput('flexPlotElement')
                        )
                    ), options = list(delay = '1000')
                )
            )
        ),
        tabItem(
            tabName = 'help',
            h2('Help documentation Goes Here')
        ),
        tabItem(
            tabName = 'summary',
            fluidRow(
                infoBoxOutput('progressBox1'),
                infoBoxOutput('progressBox2'),
                infoBoxOutput('progressBox3')
            ),
            fluidRow(
                infoBoxOutput('approvalBox1'),
                infoBoxOutput('approvalBox2'),
                infoBoxOutput('approvalBox3')
            )
        ),
        tabItem(
            tabName = 'fragview',
            fluidRow(
                nglShinyOutput('FragViewnglShiny', height = '500px'),
                jqui_draggable(tabBox(
                    div(style='overflow-y:scroll;height:600px;',DT::dataTableOutput('therow')), width=10)), textOutput('fv_warn')
            )
        ),
        tabItem(
            tabName = 'launchpad',
            uiOutput('launchpad_stuff')
        )
    )
)

ui <- dashboardPage(header, sidebar, body)

server <- function(input, output, session){
    sID <- sample(1:100000, 1)
    debug = TRUE
    if(debug) debugMessage(sID=sID, sprintf('Session init'))
    session$allowReconnect(FALSE)
    sessionDisconnect <- function(){
        # Clean-up:
        try(file.remove(sprintf('/dls/science/groups/i04-1/software/xchemreview/www/report_%s.pdf', sID)), silent=T)
        debugMessage(sID=sID, 'Disconnected')
    }
    session$onSessionEnded(sessionDisconnect)
    #observe({
    query <- parseQueryString(isolate(session$clientData$url_search))
    if(!is.null(query[['target']])){
        print(query[['target']])
        target_list <- sort(c(query[['target']]))
        fragfolders <- c('', target_list)
    }
    #})
    epochTime <- function() as.integer(Sys.time())
    humanTime <- function() format(Sys.time(), "%Y%m%d%H%M%OS")
    sessionTime <- reactive({epochTime()})
    sendEmail <- function(structure, user, decision, reason, comments){
        protein <- gsub('-[a-zA-Z]*[0-9]+_[0-9]*[A-Z]*', '', structure)
        if(is.null(emailListperStructure[[protein]])){
            emaillist <- defaultUsers
        } else {
            emaillist <- emailListperStructure[[protein]]
        }
        print(emaillist)
        if(debug) debugMessage(sID=sID, sprintf('Sending Email'))
        sendmailR::sendmail(
            from = '<XChemStructureReview@diamond.ac.uk>',
            to = sort(unique(emaillist)),
            subject = sprintf('%s has been labelled as %s', structure, decision),
            msg = sprintf(
'%s has been labelled as %s by %s for the following reason(s): %s.
With these additional comments:
%s
-------------------------------
If you wish to review this change please go to xchemreview.diamond.ac.uk while
connected to the diamond VPN or via NX.
Direct Link (must be connected to diamond VPN): https://xchemreview.diamond.ac.uk/?xtal=%s&protein=%s
This email was automatically sent by The XChem Review app
If you believe you have been sent this message in error, please email tyler.gorrie-stone@diamond.ac.uk',
            structure, decision, user, reason, comments, structure, protein),
            control = list(
                smtpServer = 'exchsmtp.stfc.ac.uk',
                smtpPort = 25
                )
        )
        modalDialog(title = 'Submission Sent', 'Your review has been sucessfully recorded, please select another structure!', footer = modalButton("Dismiss"),
  size = 's', easyClose = TRUE, fade = TRUE)
    }

    resetForm <- function(){
        if(debug) debugMessage(sID=sID, sprintf('Resetting Form'))
        updateSelectizeInput(session, "ligand", selected = '', choices = sort(rownames( r1() )))
        updateSelectInput(session, 'decision', selected ='', choices = possDec)
        updateSelectInput(session, 'reason', selected='', choices='')
        updateTextInput(session, 'comments', value = "")
        return(restartSessionKeepOptions())
    }


    possDec <- c("", "Release", "More Refinement", "More Experiments", "Reject")
    possAns <- possAns2 <- c('Select Decision')
    possRes <- list()
    possRes[['Release']] <- c('High Confidence', 'Clear Density, Unexpected Ligand', 'Correct Ligand, Weak Density', 'Low Confidence', 'No Ligand Present')
    possRes[["More Refinement"]] <- c('Check ligand conformation',
        'Check sidechain rotamers',
        'Check Rfactors',
        'Check that refinement converged',
        'Improve water model',
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
    possRes[["Reject"]] <- c('Density not convincing',
        'Too few interactions',
        'Binding site too noisy',
        'Not the ligand',
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

    output$fv_warn <- renderPrint({sessionlist$fv_warn})

    # Loading Data Gubbins:
    restartSessionKeepOptions <- function(){
        message('Updating Data')
        dbdat <- getReviewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password, target_list=target_list)
        print(dim(dbdat))
        inputData <- reactive({dbdat})
        return(inputData)
    }

    reactiviseData <- function(inputData, input){
        reactive({
            rowidx <- rep(FALSE, nrow(inputData()))
            outcome <- as.numeric(as.character(inputData()$outcome))
            if(any(c(is.null(input$protein), is.null(input$out4), is.null(input$out5), is.null(input$out6)))){
                inputData()[,]
            } else if(input$protein == '') {
                if(input$out4) rowidx[outcome==4] <- TRUE
                if(input$out5) rowidx[outcome==5] <- TRUE
                if(input$out6) rowidx[outcome==6] <- TRUE
                inputData()[rowidx,]
            } else {
                if(input$out4) rowidx[outcome==4] <- TRUE
                if(input$out5) rowidx[outcome==5] <- TRUE
                if(input$out6) rowidx[outcome==6] <- TRUE
                inputData()[rowidx & grepl(input$protein, as.character(inputData()$target_name)),]
            }
        })
    }

    updateMainTable <- function(r1, pl=25){
        DT::renderDataTable({
            DT::datatable(
                r1(), callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
                selection = 'single',
                options = list(
                    pageLength = pl,
                    columnDefs = list(list(width='100px', targets=c(4)))
                ), rownames= FALSE
            ) %>% DT::formatStyle(
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
        })
    }

    updateMainTable2 <- function(r1, pl=25){
        DT::renderDataTable({
            DT::datatable(
                r1(), callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
                selection = 'single',
                options = list(
                    pageLength = pl
                ), , rownames= TRUE
            ) %>% DT::formatStyle(columns = 1:ncol(r1()),"white-space"="nowrap")
        })
    }

    updateFlexPlot <- function(flexdata){
        renderPlotly({
            plot_ly(flexdata()$data, x=~x, y=~y, text=~ligand_name, color=~status, customdata = ~ligand_name, size=20) %>% config(scrollZoom = TRUE)
        })
    }


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
                print(vec)
                list(
                    col=vec
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
    }

    # Selector Stuff:
    review_data <- getReviewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password, target_list=target_list)
    mtzzz <- review_data[,c('crystal_name', 'mtz_latest')]
    updateSelectInput(session, 'protein', selected = '', choices=c('', sort(unique(as.character(review_data$target_name)))))
    updateSelectInput(session, 'fpe_target', selected = '', choices=c('', sort(unique(as.character(review_data$target_name)))))

    inputData <- restartSessionKeepOptions()
    r1 <- reactiviseData(inputData=inputData, input=input)
    output$reviewtable <- updateMainTable(r1=r1)
    reviewtableproxy <- DT::dataTableProxy('reviewtable')
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

    react_fv_data <- function(data, input){
        reactive({
            if(!is.null(input$fragSelect)){
                if(!input$fragSelect == '' | input$fragSelect == 'Select'){
                    toshow <- data()$targetname == input$fragSelect
                    dat <- data()[toshow, c('original_name', 'fragalysis_name', 'alternate_name','Site_Label', 'new_smiles', 'pdb_id')]
                    rn <- rownames(dat)
                    rownames(dat) <- gsub('-[0-9]$', '', rn)
                    dat
                } else {
                    data()[NULL,]
                }
            }
        })
    }

    reactivegetFragalysisViewData <- function(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password, target_list=target_list){
        reactive({getFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password, target_list=target_list)})
    }

    fvd <- getFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password, target_list=target_list)
    fragview_data <- reactivegetFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password, target_list=target_list)
    #fragfolders <- c('', sort(unique(fvd$targetname)))
    #fragfolders <- c('', 'Mpro', 'PlPro', 'PHIPA', 'NSP16','PGN_RS02895PGA', 'XX02KALRNA')
    updateSelectInput(session, 'fragSelect', selected='', choices=fragfolders)

    fragview_input <- react_fv_data(fragview_data, input)
    fragview_table_data <- react_fv_data2(fragview_data, input)

    output$therow <- updateMainTable2(fragview_input, pl=100)
    fragviewproxy <- DT::dataTableProxy('therow')
    output$as_message <- renderText({'Alt Click to select Atom'})

    observeEvent(input$as_clear, {
        session$sendCustomMessage(type = 'as_resetclicked', list())
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
        if(debug) debugMessage(sID=sID, sprintf('Selecting: %s', input$fragSelect))
        #folderPath <- getAlignedStagingFolder()
        if(input$fragSelect == ''){
            message('nothing selected')
        } else {
            fv_values$apofiles <- as.character(isolate(fragview_table_data()$apo_pdb))
            apo_existing <- sapply(fv_values$apofiles, file.exists)
            fv_values$molfiles <- as.character(isolate(fragview_table_data()$lig_mol_file))
            #molfiles_existing <- sapply(fv_values$molfiles, file.exists)
            fv_values$apofiles <- fv_values$apofiles[apo_existing]
            fv_values$molfil <- gsub('.mol', '', basename(fv_values$molfiles))
            updateSelectInput(session, 'goto', choices = fv_values$molfil)
            fragview_input <- react_fv_data(fragview_data, input) # Filter missing files here??
            fragviewproxy %>% replaceData(fragview_input(), rownames = TRUE, resetPaging = FALSE)
            #output$therow <- updateMainTable2(fragview_input, pl=100)
            tryAddPDB <- try(uploadApoPDB(filepath=fv_values$apofiles[1], repr='cartoon', focus=TRUE), silent=T)
            molout <- try(sapply(fv_values$molfiles, uploadUnfocussedMol), silent=T)
        }   
    })

    observeEvent(input$gonext, {
        molfiles <- fv_values$molfiles
        molbase <- fv_values$molfil
        nmol <- length(molfiles)
        id <- which(molbase == input$goto)
        next_id <- id + 1
        if(next_id > nmol) next_id <- 1 # Overflow back to start of list
        # Cycle along to next ligand in molfil
        if(debug) debugMessage(sID=sID, sprintf('Switching to: %s', molbase[next_id]))
        updateSelectInput(session, 'goto', selected = molbase[next_id], choices=molbase)
    })

    observeEvent(input$goback, {
        molfiles <- fv_values$molfiles
        molbase <- fv_values$molfil
        nmol <- length(molfiles)
        id <- which(molbase == input$goto)
        next_id <- id - 1
        if(next_id < 1) next_id <- nmol # Underflow to end of list
        if(debug) debugMessage(sID=sID, sprintf('Switching to: %s', molbase[next_id]))
        updateSelectInput(session, 'goto', selected = molbase[next_id], choices=molbase)
    })


    uploadMolAndFocus2 <- function(filepath){
        syscall <- sprintf('cat %s', filepath)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse = '\n')
        session$sendCustomMessage(
            type = 'fv_addMolandfocus',
            list(choice)
        )
    }

    observeEvent(input$goto, {
        output$metastatus <- renderText({'STATUS: Pending...'})
        molfiles <- fv_values$molfiles
        molbase <- fv_values$molfil
        names(molfiles) <- molbase
        if(!input$goto == ''){
            folder <- dirname(molfiles[input$goto])
            mol_file <- molfiles[input$goto]
            smi_file <- gsub('.mol', '_smiles.txt', mol_file)
            smilestr <- system(sprintf('cat %s', smi_file), intern=T)

            choices <- unique(c('', as.character(isolate(fragview_table_data()[, 'Site_Label']))))
            # Fill Form as seen
            updateTextInput(session, 'crysname', value = input$goto)
            updateTextInput(session, 'smiles', value = smilestr)
            updateTextInput(session, 'new_smiles', value = as.character(isolate(fragview_table_data()[input$goto, 'new_smiles'])))
            updateTextInput(session, 'alternate_name', value = as.character(isolate(fragview_table_data()[input$goto, 'alternate_name'])))
            updateSelectizeInput(session, 'site_name', selected = as.character(isolate(fragview_table_data()[input$goto, 'Site_Label'])), choices=choices)
            updateTextInput(session, 'pdb_entry', value = as.character(isolate(fragview_table_data()[input$goto, 'pdb_id'])))
            # Go to specific ligand do not edit go next loop
            if(debug) debugMessage(sID=sID, sprintf('Selected: %s', input$goto))
            if(debug) debugMessage(sID=sID, sprintf('trying to view: %s', molfiles[input$goto]))
            gogogo <- try(uploadMolAndFocus2(mol_file), silent=T)
            if(!file.exists(mol_file)){
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
        fragview_data <- reactivegetFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password, target_list=target_list)
        fragview_input <- react_fv_data(fragview_data, input)
        fragview_table_data <- react_fv_data2(fragview_data, input)
        fragviewproxy %>% replaceData(fragview_input(), rownames = TRUE, resetPaging = FALSE)
        #output$therow <- updateMainTable2(fragview_input, pl=100)
    })

    observeEvent(input$write, {
        if(debug) debugMessage(sID=sID, sprintf('Writing a row'))
        rn <- rownames(isolate(fragview_table_data()))
        ids <- isolate(fragview_table_data()$id)
        names(ids) <- rn
        updateOrCreateRow(ligand_name_id=as.character(ids[input$goto]),
                          fragalysis_name=as.character(input$goto),
                          original_name=as.character(rsplit(input$goto, '_', 1)[1]),
                          site_label=as.character(input$site_name),
                          new_smiles=as.character(input$new_smiles),
                          alternate_name=as.character(input$alternate_name),
                          pdb_id=as.character(input$pdb_entry),
                          dbname=db,
                          host=host_db,
                          port=db_post,
                          user=db_user,
                          password=db_password)
        output$metastatus <- renderText({'STATUS: Written!'})
        if(!input$desync){
            fragview_data <- reactivegetFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password)
            fragview_input <- react_fv_data(fragview_data, input)
            fragview_table_data <- react_fv_data2(fragview_data, input)
            fragviewproxy %>% replaceData(fragview_input(), rownames = TRUE, resetPaging = FALSE)
            #output$therow <- updateMainTable2(fragview_input, pl=100)
        }
    })

    # On Table Rowclick # Potentially slow? Unneeded? # Go back to
    observeEvent(input$therow_rows_selected, {
        if(debug) debugMessage(sID=sID, sprintf('Selecting Row'))
        molfiles <- fv_values$molfiles
        molbase <- fv_values$molfil
        choice = isolate(rownames(fragview_input())[input$therow_rows_selected])
        updateSelectizeInput(session, 'goto', selected = choice, choices=molbase)
    })

    observeEvent(input$protein,{
        updateSelectInput(session, 'fpe_target', selected=input$protein)
    })
    observeEvent(input$fpe_target,{
        updateSelectInput(session, 'protein', selected=input$fpe_target)
    })

    output$progressBox1 <- renderInfoBox({
        infoBox('Total Reviewed', isolate(interestingData()$total_reviewed), icon = icon('thumbs-up', lib = 'glyphicon'), color = 'red')
    })
    output$approvalBox1 <- renderInfoBox({
        infoBox('Number of Ligands', isolate(interestingData()$total_ligands), icon = icon('thumbs-up', lib = 'glyphicon'), color = 'red')
    })

    #total_fragalysis_ligands =
    #total_fragalysis_ligands_to_annotate =

    sessionGreaterThanMostRecentResponse <- function(id, sessionTime){
        con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
        response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses_new"))
        dbDisconnect(con)
        if(nrow(response_data) > 0){
            mostrecent <- as.data.frame(t(sapply(split(response_data, response_data$Ligand_name_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
            rownames(mostrecent) <- as.character(mostrecent$Ligand_name_id)
            t0 <- mostrecent[as.character(id), 'time_submitted'][[1]]
            output <- ifelse(is.null(t0), TRUE, sessionTime > t0)
        } else {
            output <- TRUE
        }
        return(TRUE) # Just force it... The app updates fairly rapidly now...
    }

    displayModalWhoUpdated <- function(id){
                con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
                response_data <- dbGetQuery(con, sprintf("SELECT * FROM review_responses_new"))
                dbDisconnect(con)

                mostrecent <- as.data.frame(t(sapply(split(response_data, response_data$Ligand_name_id), function(x) x[which.max(x$time_submitted),])), stringsAsFactors=F)
                rownames(mostrecent) <- as.character(mostrecent$Ligand_name_id)
                user <- mostrecent[as.character(id), 'fedid'][[1]]

                showModal(modalDialog(title = "Someone has recently reviewed this crystal",
                    sprintf("A User (%s) has recently reviewed this structure. Restarting the session to update their response. If you disagree with the current response, please submit another response or select another crystal.", user)
                    , footer = tagList(actionButton("ok", "Okay"))
                ))
    }

        # Save Responses.
    saveData <- function(data, xtaln, atoms) {
        con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
        dbAppendTable(con, 'review_responses_new', value = data, row.names=NULL)
        dbDisconnect(con)
        # Check for Atoms and tag them to review response?
        rr <- getReviewRow(data, db = db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password)
        for(atom in seq_len(nrow(atoms))){
            message(atom)
            newdat <- cbind(atomid=as.numeric(atoms[atom,2]), comment=atoms[atom,3], rr)
            str(newdat)
            colnames(newdat) <- c('atomid', 'comment', 'Ligand_id', 'Review_id')
            newdat$comment = as.character(newdat$comment)
            newdat$Ligand_id = as.numeric(newdat$Ligand_id)
            newdat$Review_id = as.numeric(newdat$Review_id)
            str(newdat)
            con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
            dbAppendTable(con, 'BadAtoms', value = newdat, row.names=NULL)
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
        sendEmail(xtaln, data[,'fedid'], data[,'decision_str'], data[,'reason'], data[,'comment'])
    }

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
        inputData <- restartSessionKeepOptions()
        r1 <- reactiviseData(inputData=inputData, input=input)
        reviewtableproxy %>% replaceData(r1(), rownames = FALSE, resetPaging = FALSE)
        #output$reviewtable <- updateMainTable(r1=r1)
        flexplotData <- flexPlotDataFun(r1=r1, input=input)
        output$flexplot1 <- updateFlexPlot(flexdata=flexplotData)
        sessionTime <- reactive({epochTime()})
        removeModal()
    })

    # Upon Main Page Submit
    observeEvent(input$submit, {
        fData <- formData()[[1]]
        xtaln <- formData()[[2]]
        if(debug) debugMessage(sID=sID, sprintf('Submitting Form'))
        if(debug) print(fData)
        if(any(fData[1:7] %in% c('', ' '))) {
            showModal(modalDialog(title = "Please fill all fields in the form",
                "One or more fields have been left empty. Please provide your FedID, a decision and reason(s) before clicking submit.",
                easyClose=TRUE, footer = tagList(modalButton("Cancel"))
            ))
        } else {
             # Get ID...
            xId <- fData[ ,'Ligand_name_id']
            # Check ID
            if(sessionGreaterThanMostRecentResponse(id=xId, sessionTime=sessionTime())){
                print(atomstoquery$data)
                if(any(as.character(atomstoquery$data$comment) %in% c('', ' '))){
                    showModal(modalDialog(title = 'You have flagged some atoms',
                        'Please annotate the selected atoms in the Atom Selection tab by double clicking on the comment cells. If you accidentally flagged an atom, try reloading the structure and resubmitting your review!',
                        easyClose=TRUE))
                } else {
                    saveData(fData, xtaln, atomstoquery$data)
                    message(sessionTime())
                    inputData <- resetForm()
                    r1 <- reactiviseData(inputData=inputData, input=input)
                    #output$reviewtable <- updateMainTable(r1=r1)
                    reviewtableproxy %>% replaceData(r1(), rownames = FALSE, resetPaging = FALSE)
                    flexplotData <- flexPlotDataFun(r1=r1, input=input)
                    output$flexplot1 <- updateFlexPlot(flexdata=flexplotData)
                    sessionTime <- reactive({epochTime()})
                }
            } else {
                displayModalWhoUpdated(id=xId)
            }
        }
    })

    atomstoquery <- reactiveValues()
    atomstoquery$data <- data.frame(name=character(),
                 index=character(),
                 comment=character(),
                 stringsAsFactors=FALSE)

    output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))})

    observeEvent(input$clickedAtoms, {
        newdat <- isolate(atomstoquery$data)
        # Check for 'new' rows:
        new <- which(!as.character(input$clickNames) %in% as.character(newdat$name))
        for(i in new){
            newdat <- rbind(newdat, data.frame(name = input$clickNames[i], index = input$clickedAtoms[i], comment = '', stringsAsFactors=FALSE))
        }
        tokeep <- as.character(newdat$name) %in% as.character(input$clickNames)
        newdat <- newdat[tokeep,]
        atomstoquery$data <- newdat
        print(atomstoquery$data)
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))), options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))})
    })

    observeEvent(input$atoms_cell_edit, {
        info = input$atoms_cell_edit
        str(info)
        i = info$row
        j = info$col
        v = info$value
        update <- isolate(atomstoquery$data)
        update[i, j] <- as.character(v)
        atomstoquery$data <- update
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))), options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))})
    })


    observeEvent(input$write_selected, {
        idx <- isolate(input$atoms_rows_selected)
        update <- isolate(atomstoquery$data)
        update[idx, 3] <- as.character(input$atom_text)
        atomstoquery$data <- update
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))), options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))})
    })

    observeEvent(input$write_all,{
        update <- isolate(atomstoquery$data)
        idx <- 1:nrow(update)
        update[idx, 3] <- as.character(input$atom_text)
        atomstoquery$data <- update
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))), options = list(autoWidth = TRUE, columnDefs = list(list(width='50px', targets=c(1,2)))))})
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

    # Actual Stuff
    # Flexible Sidebar Stuff
    output$flex1 <- renderUI({
        switch(input$tab,
            review = tagList(
                selectInput('protein', 'Select Protein', selected = '', choices=c('', sort(unique(as.character(review_data$target_name))))),
                div(
                    id = 'form',
                    # Ligand/Xtal Select????
                    textInput('name', 'Name/FedID', ''),
                    selectInput('ligand', 'Ligand', selected='', choices = rownames(isolate(r1())), multiple=FALSE),
                    selectInput("decision", "Decision", choices = possDec),
                    selectizeInput("reason", "Reason(s)", list(), multiple=TRUE),
                    textInput('comments', 'Additional  Comments', value = "", width = NULL,placeholder = NULL),
                    fluidRow(
                        column(6, actionButton('submit', 'Submit', class = 'btn-primary')),
                        column(6, actionButton('clear', 'Clear', class = 'btn-primary'))
                    )
                ),
                fluidRow(
                    column(4, checkboxInput('out4', 'Comp Chem Ready', value = TRUE)),
                    column(4, checkboxInput('out5', 'Deposition Ready', value = FALSE))
                ),
                fluidRow(
                    column(4, checkboxInput('out6', 'Deposited', value = FALSE))
                )
            ),
            fragview = tagList(
                    selectInput('fragSelect', 'Project Select', selected = '', choices=fragfolders),
                    checkboxInput('desync', 'Turn off automatic Updates', value = FALSE),
                    actionButton('goback', 'Prev Ligand'),
                    actionButton('gonext', 'Next Ligand'),
                    selectInput('goto', 'Go to Ligand', choices=list()),
                    textInput('crysname', 'Ligand Name', '' ),
                    textInput('smiles', 'Smiles String', ''),
                    textInput('new_smiles', 'New Smiles String', ''),
                    textInput('alternate_name', 'Alternate Fragment Name', ''),
                    selectizeInput('site_name', 'Site Label', list(), multiple=FALSE, options=list(create=TRUE)),
                    textInput('pdb_entry', 'PDB Entry', ''),
                    textOutput('metastatus'),
                    uiOutput('writeButton'),
                        #hr(),
                        #textInput('newCrystalName', 'New Fragment Name', ''),
                        #actionButton('changeName', 'Change Name of Fragment (will not assign site label)'),
                        #hr(),
                        #textOutput('massChange'),
                        #selectizeInput('site_name2', 'Old label', list(), multiple=FALSE),
                        #textInput('new_label', 'New label', ''),
                        #actionButton('mcl', 'Mass Convert Label'),
                        #hr(),
                    actionButton('updateTable', 'Refresh Metadata Table')
                #selectInput('b1', 'Selection', c('setosa', 'versicolor', 'virginica'))
            ),
            help = tagList(
                selectInput('c1', 'Selection', c('setosa', 'versicolor', 'virginica'))
            ),
            launchpad = tagList(
                selectInput('d1', 'Selection', c('setosa', 'versicolor', 'virginica'))
            ),
            summary = tagList(
                selectInput('protein_to_summarize', 'Selection', selected = '', choices=sort(unique(as.character(review_data$target_name))))
            )
        )
    })

    observeEvent(input$decision,{
        possAns <- possRes[[input$decision]]
        updateSelectizeInput(session,'reason', choices=possAns)
    })

    # Default Values for Control panel trick
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

    # Reads whatever the current input values are...
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

    # On session init, set control panel values to defaults.
    ngl_control_values <- reactiveValues()
    ngl_control_values$defaults <- loadDefaultParams()

    # Control Panel Listeners
    observeEvent(input$controls, ignoreNULL = TRUE, {
        showModal(
            controlPanelModal(
                values = isolate(ngl_control_values$defaults),
                title = 'NGL Viewer Controls'
            )
        )
    })

    observeEvent(input$tab, ignoreNULL=TRUE, {
        if(input$tab == 'review'){
            showModal(
                controlPanelModal(
                    values = isolate(ngl_control_values$defaults),
                    title = 'As part of setup please confirm NGL Viewer Controls'
                )
            )
        }
    })


    observeEvent(input$buster, ignoreNULL=TRUE, {
    	pdf_files = list.files(sessionlist$xtalroot, rec=T, pattern='report.pdf', full=T)
    	pdf_file = tail(pdf_files,1)
        message(pdf_file)
        copied <- file.copy(from=pdf_file, to=sprintf('/dls/science/groups/i04-1/software/xchemreview/www/report_%s.pdf', sID), overwrite=TRUE)
        message(copied)
        addResourcePath("www", '/dls/science/groups/i04-1/software/xchemreview/www')
    	output$pdfview <- renderUI({
		    #includeHTML(file.path('pdf_folder','index.html'))
      		tags$iframe(style="height:800px; width:100%", src=sprintf('www/report_%s.pdf', sID))
    	})
    	showModal(
    		customDraggableModalDialog(title=pdf_file, 
                if(length(pdf_file)>0) uiOutput("pdfview")
                else 'Unable to find Buster Report', size='l', easyClose=FALSE)
   	    )
    })

    observeEvent(input$updateParams, {
        removeModal()
        for(i in names(ngl_control_values$defaults)){
            ngl_control_values$defaults[[i]] <- input[[i]]
        }
    })

    updateParam <- function(which, what) session$sendCustomMessage('updateaparam', list(which, what))

    observeEvent(input$backgroundColor, { updateParam('backgroundColor', as.character(input$backgroundColor)) })
    observeEvent(input$cameraType, { updateParam('cameraType', as.character(input$cameraType)) })
    observeEvent(input$mousePreset, { updateParam('mousePreset', as.character(input$mousePreset)) })
    observeEvent(input$clipDist, { updateParam('clipDist', as.character(input$clipDist)) })
    observeEvent(input$fogging, {
        updateParam('fogNear', as.character(input$fogging[1]) )
        updateParam('fogFar' , as.character(input$fogging[2]) )
    })
    observeEvent(input$clipping, {
        updateParam('clipNear', as.character(input$clipping[1]) )
        updateParam('clipFar' , as.character(input$clipping[2]) )
    })


    # NGL Shiny Stage...
    output$nglShiny <- renderNglShiny(
        nglShiny(name = 'nglShiny', list(), width = NULL, height = NULL)
    )

    output$FragViewnglShiny <- renderNglShiny(
        nglShiny(name = 'nglShiny', list(), width=NULL, height=100)
    )

    # PBD uploader
    removeNamedComponent <- function(objectname) session$sendCustomMessage(type='removeNamedComponent', list(objectname))

    uploadPDB <- function(filepath, input){
        syscall <- sprintf('cat %s', filepath)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse = '\n')
        session$sendCustomMessage(
            type = 'setPDB2', # See TJGorrie/NGLShiny for details on setPDB2
            message = list(
                choice,
                isolate(input$clipDist),
                isolate(input$clipping)[1],
                isolate(input$clipping)[2],
                isolate(input$fogging)[1],
                isolate(input$fogging)[2]
            )
        )
    }


    tcl <- function(x) tolower(as.character(as.logical(x)))

    uploadApoPDB <- function(filepath, repr, focus){
        syscall <- sprintf('cat %s', filepath)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse = '\n')
        session$sendCustomMessage(
            type = 'setapoPDB',
            message = list(
                choice,
                repr,
                tcl(focus)
            )
        )
    }

    uploadMolAndFocus <- function(filepath, ext, focus){
        syscall <- sprintf('cat %s', filepath)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse = '\n')
        session$sendCustomMessage(
            type = 'addMolandfocus',
            list(choice,ext, tcl(focus))
        )
    }

    uploadUnfocussedMol <- function(filepath){
        syscall <- sprintf('cat %s', filepath)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse = '\n')
        session$sendCustomMessage(
            type='addMol',
            list(choice)
        )
    }


    getExt <- function(x) sapply(strsplit(x, '[.]'), tail, 1)

    # Map Uploader
    uploadVolumeDensity <- function(filepath, color, negateiso = FALSE, boxsize, isolevel, visable, windowname, isotype='value'){
        volume_bin <- readBin(filepath, what='raw', file.info(filepath)$size)
        volume_b64 <- base64encode(volume_bin, size=NA, endian=.Platform$endian)
        session$sendCustomMessage(
            type = 'addVolumeDensity',
            message = list(
                as.character(volume_b64),
                as.character(isolevel),
                as.character(color),
                tcl(negateiso),
                as.character(getExt(filepath)),
                as.character(boxsize),
                tcl(visable),
                as.character(windowname),
                as.character(isotype)
            )
        )
    }

    # Control how volume densities get toggled.
    updateVisability <- function(name, bool){
        session$sendCustomMessage(
            type = 'updateVolumeDensityVisability',
            list(
                as.character(name),
                tcl(bool)
            )
        )
    }

    # Map Listeners
    observeEvent(input$eventMap,   { updateVisability('eventmap', input$eventMap  ) })
    observeEvent(input$twofofcMap, { updateVisability('twofofc' , input$twofofcMap) })
    observeEvent(input$fofcMap,    {
        updateVisability('fofcpos', input$fofcMap)
        updateVisability('fofcneg', input$fofcMap)
    })

    updateDensityISO <- function(name, isolevel) session$sendCustomMessage('updateVolumeDensityISO', list(name, isolevel))

    observeEvent(input$isoEvent, {updateDensityISO('eventmap', input$isoEvent)})
    observeEvent(input$iso2fofc, {updateDensityISO('twofofc', input$iso2fofc)})
    observeEvent(input$isofofc , {
        updateDensityISO('fofcpos', input$isofofc)
        updateDensityISO('fofcneg', input$isofofc)
    })

    updateDensityBoxSize <- function(name, boxsize) session$sendCustomMessage('updateVolumeDensityBoxSize', list(name, boxsize))
    observeEvent(input$boxsize , {
        for(windowname in c('eventmap', 'twofofc', 'fofcpos', 'fofcneg')) updateDensityBoxSize(windowname, input$boxsize)
    })

    observeEvent(input$fitButton, {
        #try(uploadMolAndFocus(session = session, filepath = isolate(session_data$selected)$mol_file, ext = 'mol', focus=TRUE), silent = TRUE)
        try(session$sendCustomMessage(type='focus_on_mol', list()), silent=TRUE)
    })

    observeEvent(input$asuSwitch, {
        try(session$sendCustomMessage('updateAssembly', list(isolate(input$asuSwitch))))
    })

    output$isoEventSlider <- renderUI({
            sliderInput("isoEvent", "",
                    min = 0, max = 3,
                    value = 1, step = 0.1)
    })

    output$iso2fofcSlider <- renderUI({
            sliderInput("iso2fofc", "",
                    min = 0, max = 3,
                    value = 1.5, step = 0.1)
    })

    output$isofofcSlider <- renderUI({
            sliderInput("isofofc", "",
                min = 0, max = 3,
                value = 3, step = 0.1)
    })

    observeEvent(input$reviewtable_rows_selected, {
        rdat <- r1()[input$reviewtable_rows_selected,,drop=TRUE]
        print(rdat)
        sessionlist$rowname <- rownames(r1())[input$reviewtable_rows_selected]
        sessionlist$lig_name <- rdat$ligand_name
        sessionlist$lig_id <- rdat[[1]][1]
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
        if(is.null(input$views)) updateRadioButtons(session, 'views', selected = 'aligned')

        session$sendCustomMessage(type = 'setup', message = list())
        updateParam('mousePreset', as.character(input$mousePreset))
        updateParam('clipDist', as.character(input$clipDist))
        updateSelectInput(session, 'emap', choices = c('NotAMap.ccp4'), selected = c('NotAMap.ccp4'))
        updateSelectInput(session, 'asuSwitch', selected='AU', choices=c('AU', 'UNITCELL', 'SUPERCELL'))

        the_pdb_file <- isolate(sessionlist$apo_file)
        the_mol_file <- isolate(sessionlist$mol_file)
        the_emaps <- dir(dirname(isolate(sessionlist$apo_file)), pattern='event', full=TRUE)
        the_2fofc_map <- isolate(sessionlist$twofofc_file)
        the_fofc_map <- isolate(sessionlist$fofc_file)

        if(input$renderMisc & !isolate(sessionlist$apo_file) == ""){
            #spfile <- tail(dir(isolate(sessionlist$xtalroot), pattern='A-1101.png', full.names=T, rec=T),1)
            #output$spiderPlot <- renderImage({
            #    if(length(spfile) == 1){
            #        list(src = spfile, contentType = 'image/png', width=200, height=200)
            #    } else {
            #        list(src = '', contentType = 'image/png', width=200, height=200)
            #    }
            #}, deleteFile=FALSE)
            ligfile <- tail(dir(sprintf('%s/compound', isolate(sessionlist$xtalroot)), pattern = '.png', full.names=T),1)
            rendered_ligfile <- tail(dir(dirname(isolate(sessionlist$apo_file)), pattern='.png', full.names = TRUE), 1)
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
            output$rendered_ligimage <- renderImage({
                if(length(rendered_ligfile) == 1){
                    list(src=rendered_ligfile, contentType='image/png', width=200, height=200)
                } else {
                    list(src = '', contentType = 'image/png', width=200, height=200)
                }
            }, deleteFile=False)
        }

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
                        try(uploadApoPDB(the_pdb_file, 'line', focus=input$autocenter), silent=T)
                        try(uploadMolAndFocus(the_mol_file, 'mol', focus=input$autocenter), silent=T)
                        session$sendCustomMessage(type = 'restore_camera_pos', message = list())
                    },
                    'unaligned' = {
                        session$sendCustomMessage(type = 'save_camera_pos', message = list())
                        the_pdb_file <- gsub('staging_test', 'unaligned_test', the_pdb_file)
                        the_mol_file <- gsub('staging_test', 'unaligned_test', the_mol_file)
                        the_emaps <- dir(dirname(the_pdb_file), pattern='event', full=TRUE)
                        the_2fofc_map <- gsub('staging_test', 'unaligned_test', the_2fofc_map)
                        the_fofc_map <- gsub('staging_test', 'unaligned_test', the_fofc_map)
                        try(uploadApoPDB(the_pdb_file, 'line', focus=TRUE), silent=T)
                        try(uploadMolAndFocus(the_mol_file, 'mol', focus=TRUE), silent=T)
                    },
                    'crystallographic' = {
                        session$sendCustomMessage(type = 'save_camera_pos', message = list())
                        splitted <-  rsplit(the_pdb_file, '/')
                        the_folder <- dirname(gsub('aligned', 'crystallographic', splitted[1]))
                        the_xtal_name <- gsub('_[0-9][A-Z]_apo.pdb', '', splitted[2])
                        the_pdb_file <- sprintf('%s/%s.pdb', the_folder, the_xtal_name)
                        try(uploadPDB(the_pdb_file, input=input), silent=T)
                        the_2fofc_map <- sprintf('%s/%s_2fofc.map', the_folder, the_xtal_name)
                        the_fofc_map <- sprintf('%s/%s_fofc.map', the_folder, the_xtal_name)
                        the_emaps <- dir(the_folder, pattern=sprintf('%s_event', the_xtal_name), full=TRUE)
                    }
                )

                names(the_emaps) <- basename(the_emaps)
                sessionlist$current_emaps <- the_emaps
                print(the_emaps)
                incProgress(.2, detail = 'Uploading Event map')
                updateSelectInput(session, 'emap', choices = names(isolate(sessionlist$current_emaps)), selected = names(isolate(sessionlist$current_emaps))[1])
                # Move this to a different part?
                message('Upload fofcs')
                incProgress(.2, detail = 'Uploading 2fofc map')
                try(uploadVolumeDensity(the_2fofc_map,
                    color = 'blue', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$iso2fofc, visable=input$twofofcMap, windowname='twofofc', isotype='sigma'), silent=T)
                incProgress(.1, detail = 'Uploading fofc map')
                try(uploadVolumeDensity(the_fofc_map,
                    color = 'lightgreen', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcpos', isotype='sigma'), silent=T)
                incProgress(.1, detail = 'Uploading fofc map')
                try(uploadVolumeDensity(the_fofc_map,
                    color = 'tomato', negateiso = TRUE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcneg', isotype='sigma'), silent=T)
            }
            setProgress(1)
        })
    })

    observeEvent(input$ligand, ignoreNULL = TRUE, {
        atomstoquery$data <- data.frame(name=character(),
                 index=character(),
                 comment=character(),
                 stringsAsFactors=FALSE)
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data)})
        previous = isolate(input$views)
        if(previous == 'aligned'){
            session$sendCustomMessage(type = 'setup', message = list())
            updateParam('mousePreset', as.character(input$mousePreset))
            the_pdb_file <- isolate(sessionlist$apo_file)
            the_mol_file <- isolate(sessionlist$mol_file)
            the_emaps <- dir(dirname(isolate(sessionlist$apo_file)), pattern='event', full=TRUE)
            the_2fofc_map <- isolate(sessionlist$twofofc_file)
            the_fofc_map <- isolate(sessionlist$fofc_file)

            if(input$renderMisc & !isolate(sessionlist$apo_file) == ""){
                spfile <- tail(dir(isolate(sessionlist$xtalroot), pattern='A-1101.png', full.names=T, rec=T),1)
                output$spiderPlot <- renderImage({
                    if(length(spfile) == 1){
                        list(src = spfile, contentType = 'image/png', width=200, height=200)
                    } else {
                        list(src = '', contentType = 'image/png', width=200, height=200)
                    }
                }, deleteFile=FALSE)
                ligfile <- tail(dir(sprintf('%s/compound', isolate(sessionlist$xtalroot)), pattern = '.png', full.names=T),1)
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
            withProgress(message = sprintf('Loading %s Ligand', input$views), value = 0,{
                if(! isolate(sessionlist$apo_file) == ""){
                    incProgress(.2, detail = 'Uploading Crystal + Ligand')
                    try(uploadApoPDB(the_pdb_file, 'line', focus=input$autocenter), silent=T)
                    try(uploadMolAndFocus(the_mol_file, 'mol', focus=input$autocenter), silent=T)
                    names(the_emaps) <- basename(the_emaps)
                    sessionlist$current_emaps <- the_emaps
                    incProgress(.2, detail = 'Uploading Event map')
                    updateSelectInput(session, 'emap', choices = names(isolate(sessionlist$current_emaps)), selected = names(isolate(sessionlist$current_emaps))[1])
                    # Move this to a different part?
                    message('Upload fofcs')
                    incProgress(.2, detail = 'Uploading 2fofc map')
                    try(uploadVolumeDensity(the_2fofc_map,
                        color = 'blue', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$iso2fofc, visable=input$twofofcMap, windowname='twofofc'), silent=T)
                    incProgress(.1, detail = 'Uploading fofc map')
                    try(uploadVolumeDensity(the_fofc_map,
                        color = 'lightgreen', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcpos'), silent=T)
                    incProgress(.1, detail = 'Uploading fofc map')
                    try(uploadVolumeDensity(the_fofc_map,
                        color = 'tomato', negateiso = TRUE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcneg'), silent=T)
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

    })

    observeEvent(input$gotores, ignoreNULL=TRUE, {
        if(!input$gotores == ''){
            pos <- strsplit(input$gotores, '_')[[1]][2]
            print(pos)
            try(session$sendCustomMessage(type='go_to_residue', list(pos)))
        }
    })

    observeEvent(input$highlight_res, ignoreNULL=TRUE,{
        if(!input$highlight_res == ''){
            pos <- paste(sapply(strsplit(input$highlight_res, '_'), '[', 2), collapse=', ')
            print(pos)
            try(session$sendCustomMessage(type='highlight_residues', list(pos)))
        }
    })
    observeEvent(input$emap, ignoreNULL = TRUE, {
        message('Upload EMAPS')
        sel <- isolate(sessionlist$current_emaps)[input$emap]
        message(sel)
        try(uploadVolumeDensity(sel,
            color = 'orange', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$isoEvent, visable=input$eventMap, windowname='eventmap'), silent=T)
        message('Completed!')
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
            sessionlist$lig_id <- rdat[[1]][1]
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

    perc.rank <- function(x) trunc(rank(x))/length(x)
    perc.rank2 <- function(x, xo=NULL)  length(x[x <= xo])/length(x)

    plotData <- reactive({
        lapply(c('res', 'r_free', 'rcryst', 'ramachandran_outliers', 'rmsd_angles', 'rmsd_bonds'),
            function(x, alldata, extradata){
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
            }, alldata=review_data, extradata=isolate(sessionlist)
        )
    })

    output$plottoRender <- renderPlot({
        hmapbar(data=plotData(), title = isolate(sessionlist$lig_name), target_name=(isolate(sessionlist$target_name)))
    })


    ## Launch Pad Stuff:
    output$launchpad_stuff <- renderUI({
        fluidPage(
            selectInput('lp_selection','Select Target', selected = '', choices=fragfolders),
            actionButton('lp_launcher', "Launch!!!!"),
	    checkboxInput('lp_copymaps', 'Copy MapFiles?', value=TRUE),
            downloadButton("downloadFragData", "Download")
        )
    })

    observeEvent(input$lp_selection, {
        if(!isolate(input$lp_selection) == ''){
        sessionlist$fumeta <- createUniqueMetaData(db = db, host_db = host_db, db_port = db_port,
            db_user = db_user, db_password = db_password,
            target = isolate(input$lp_selection))
        }
        print(head(sessionlist$fumeta))
        message('Meta Compiled')
    })

    observeEvent(input$lp_launcher, {
        message('LAUNCH!!!')
        sessionlist$fullpath_frag <- createFragUploadFolder(meta=sessionlist$fumeta, target=isolate(input$lp_selection), copymaps=input$lp_copymaps, mtz=mtzzz)
    })

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
    autoInvalidate <- reactiveTimer(10000)
    observe({
        autoInvalidate()
        cat("")
    })
}

app <- shinyApp(ui = ui, server = server)
ip <- '0.0.0.0'
port <- '3838'
# Run App
runApp(app, host=ip, port = as.numeric(port), launch.browser = FALSE)
