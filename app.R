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
# Lib to read Binaries to encode to b64
library(caTools)
# Lib to read Pdb files
library(bio3d)
# Plotting Libs
library(ggplot2)
library(plotly)
# Weird forceful hack for the loading bar for XCR
library(waiter)
library(xfun)

read_toml <- function(file, x = xfun::read_utf8(file)){
  # extract the top-level key name, e.g., foo.bar.baz -> foo
  keys = function(x) {
    unlist(lapply(strsplit(x, '[.]'), `[[`, 1))
  }
  # generate list(name = x)
  named_list = function(x, name) setNames(list(x), name)
  # remove comments
  x = gsub('\\s+#.+', '', x)
  # arbitrary values of the form 'foo = bar' or '[foo]' or '[[foo.bar]]'
  r = '^(([[:alnum:]_]+?)\\s*=\\s*(.+)\\s*|\\[{1,2}([^]]+)\\]{1,2}(\\s*))$'
  m = regexec(r, x)
  z = lapply(regmatches(x, m), function(v) {
    if (length(v) < 6) return()
    # when data is '[foo]' instead of 'foo = bar', just return NULL
    if (v[3] == '') return(named_list(NULL, keys(v[5])))
    y = v[4]
    # strings
    if (grepl(r <- '^"([^"]*?)"$', y)) y = gsub(r, '\\1', y) else {
      # boolean
      if (y %in% c('true', 'false')) y = as.logical(y) else {
        # numbers
        if (grepl('^[0-9.]+$', y)) {
          y2 = as.numeric(y)
          if (!is.na(y2)) {
            y = y2
            y2 = as.integer(y)
            if (y2 == y) y = y2
          }
        }
      }
    }
    named_list(y, v[3])
  })
  do.call(c, z)
}

# Import a config file as specified in the shinyproxy application. 
configuration <- read_toml(tail(commandArgs(),1))

# Reinstall the latest version of nglShiny to propagate new functionality
#install.packages(configuration$nglshiny_lib, repos=NULL, type='source', lib=configuration$xcr_r_lib)
library('nglShiny', lib.loc=configuration$xcr_r_lib)

# Source components of XCR application:
# ./config.R contains key decryption function in addition to email lists...
source(file.path(configuration$xcr_path, 'config.R'))
# Stuff for XCR to run.
source(file.path(configuration$xcr_path, 'functions.R'))
source(file.path(configuration$xcr_path, 'ui.R'))
source(file.path(configuration$xcr_path, 'server.R'))

# Initialise XCR and fetch targets based on authenticated fedid.
fedid <- Sys.getenv(x = 'SHINYPROXY_USERNAME', unset = "", names = NA)

# It is likely that fetching visits is going to be the better granular unit instead of targets
# To stop people from looking at data that happens to have the misfortune of people poorly named. 
target_list <- fetchTargets(fedid=fedid, configuration = configuration)
# visit_list <- fetchVisits(fedid=fedid, configuration=configuration)
target_list2 <- target_list
fragfolders <- c('', target_list)
# Query fedid to get proposals...

waiting_screen <- tagList(
  spin_solar(),
  h4("Loading your request...")
) 

# Print Session info for debugging?
sessionInfo()

app <- shinyApp(ui = ui, server = server)
ip <- configuration$ip
xcrport <- configuration$xcrport 
# Run App
runApp(app, host=configuration$ip, port = as.numeric(xcrport), launch.browser = FALSE)
