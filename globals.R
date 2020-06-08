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
if(!local) source('/dls/science/users/mly94721/xchemreview/mailing_list.R')

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

