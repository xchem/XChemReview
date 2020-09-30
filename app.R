# App Refactor
# TODO: Scrub out filepaths add as command line options...

debug = TRUE
local = FALSE # Binding for testing on a local machine bodes no use other than me...

# Set Path: May need to add something later for files on /dls

# Server Bindings # Consider changing...
responsesDir <- '/dls/science/users/mly94721/xchemreview/Responses/' 
dataDir <- '/dls/science/users/mly94721/xchemreview/Data/'

# Load R packages...
library(devtools)
library(shiny)
library(DT)
library(htmlwidgets)
library(caTools)
library(DBI)

if(local){
    gpath <- '.'
    responsesDir <-file.path(sprintf('%s/%s', gpath, "Responses"))
    #source('./db_config.R')
    source('/dls/science/groups/i04-1/software/xchemreview/db_config.R') 
    install.packages('~/nglshiny', repos=NULL, type='source')
    library(nglShiny)
    library(httr)
} else {
    gpath <- '/srv/shiny-server/'
    # Need to move the lib location...
    #install.packages("/dls/science/groups/i04-1/software/nglshiny", repos=NULL, type='source', lib='/dls/science/groups/i04-1/software/xchemreview/xcrlib')
    library(nglShiny, lib.loc = '/dls/science/groups/i04-1/software/xchemreview/xcrlib')
    # Move this to docker...
    library(sendmailR) #, lib.loc = "/dls/science/users/mly94721/R/")
    library(httr)#, lib.loc = "/dls/science/users/mly94721/R/")
}

# Source the rest of the code
if(local){
    source('/dls/science/users/mly94721/xchemreview/slackconfig.R')
    source('./globals.R')
    source('./modules.R')
    source('./ui.R')
    source('./server.R')

} else {
    # This is not good!...
    source('/dls/science/groups/i04-1/software/xchemreview/slackconfig.R')
    source('/dls/science/groups/i04-1/software/xchemreview/globals.R')
    source('/dls/science/groups/i04-1/software/xchemreview/modules.R')
    source('/dls/science/groups/i04-1/software/xchemreview/ui.R')
    source('/dls/science/groups/i04-1/software/xchemreview/server.R')

}
# Create App

app <- shinyApp(ui = ui, server = server)
ip <- '0.0.0.0'
port <- '3838'
#cmd <- commandArgs(T)
#ip <- cmd[1]
#port <- cmd[2]

# Run App
runApp(app, host=ip, port = as.numeric(port), launch.browser = FALSE)
