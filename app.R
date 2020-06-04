# App Refactor
debug = TRUE
local = FALSE # Binding for testing on a local machine bodes no use other than me...

# Set Path: May need to add something later for files on /dls

# Server Bindings # Consider changing...
gpath <- '/srv/shiny-server/'
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

# Source the rest of the code

source('./globals.R')
source('./ui.R')
source('./server.R')

# Create App

app <- shinyApp(ui = ui, server = server)
if(local){
    ip <- '0.0.0.0'
    port <- '3838'
} else {
    cmd <- commandArgs(T)
    ip <- cmd[1]
    port <- cmd[2]
}

# Run App
runApp(app, host=ip, port = as.numeric(port), launch.browser = FALSE)
