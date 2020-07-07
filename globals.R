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

epochTime <- function() as.integer(Sys.time())
humanTime <- function() format(Sys.time(), "%Y%m%d-%H%M%OS")

getRootFP <- function(pdbpath){
    splits <- strsplit(pdbpath, split ='/')[[1]]
    n <- length(splits)
    paste0(c(splits[1:(n-2)],''), collapse='/')
}

#################################################################################
# Global Variables used in server/UI
#################################################################################
nglRepresentations = c('angle', 'axes', 'ball+stick', 'backbone', 'base', 'cartoon', 
    'contact','dihedral', 'distance', 'helixorient', 'licorice', 'hyperball', 'label',
    'line', 'surface', 'point', 'ribbon', 'rocket', 'rope', 'spacefill', 'trace', 'unitcell',
    'validation')

nglColorSchemes <-  c('atomindex','bfactor','chainid','chainindex','chainname','densityfit','electrostatic',
    'element','entityindex','entitytype','geoquality','hydrophobicity','modelindex','moleculetype','occupancy',
    'random','residueindex','resname','sstruc','uniform','value','volume')


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
    'Space_Group',
    'XCEoutcome',
    'CIF',
    'Latest.PDB',
    'Latest.MTZ',
    'Protein'
)

if(!local) source('/dls/science/users/mly94721/xchemreview/db_config.R') # Config file...
con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
refinement_data <- dbGetQuery(con, "SELECT id, crystal_name_id, r_free, rcryst, ramachandran_outliers, res, rmsd_angles, rmsd_bonds, lig_confidence_string, cif, pdb_latest, mtz_latest FROM refinement WHERE outcome=4 OR outcome=5 OR outcome=6")
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
if(!local) source('/dls/science/users/mly94721/xchemreview/mailing_list.R')

xtalList <-  sort(unique(jd[,'Xtal']))

rm(refinement_data, crystal_data, target_data, compound_data, jd, targs, comps)
gc()


defaultRepresentation <- "ball+stick"
defaultColorScheme <- "chainIndex"

possDec <- c("", "Release", "More Refinement", "More Experiments", "Reject")
possAns <- possAns2 <- c('Select Decision')
possRes <- list()
possRes[['Release']] <- c('High Confidence', 'Clear Density, Unexpected Ligand', 'Correct Ligand, Weak Density', 'Low Confidence', 'No Ligand Present')
possRes[["More Refinement"]] <- c('Check ligand conformation','Check sidechain rotamers','Check Rfactors','Check that refinement converged','Improve water model','Build alternate conformations','Fix geometry','Trim down ligand','Density did not load','Other')
possRes[["More Experiments"]] <- c('Get better density','Get better resolution','Confirm ligand identity','Check if ligand changed','Other')
possRes[["Reject"]] <- c('Density not convincing','Too few interactions','Binding site too noisy','Other')
possDec_int <- 1:4
names(possDec_int) <- c("Release", "More Refinement", "More Experiments", "Reject")

fragfolders <- dir('/dls/science/groups/i04-1/fragprep/staging/')
fragfolders <- fragfolders[!fragfolders %in% c('~', 'tmp')]

