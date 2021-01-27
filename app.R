local = FALSE

# Generic Shiny Libraries
library(httr)
library(shiny)
library(grDevices)
library(lubridate)
library(shinydashboard)
library(shinyjqui)
library(shinyWidgets)
# DB Lib
library(DBI)
library(sendmailR)
# Read Binaries to encode to b64
library(caTools)
# Render structures in ngl window
if(local) {
    library(nglShiny)
    source('./my_config.R')
} else {
    source('/dls/science/groups/i04-1/software/xchemreview/config.R')
    install.packages('/dls/science/groups/i04-1/software/nglshiny', repos=NULL, type='source', lib='/dls/science/groups/i04-1/software/xchemreview/xcrlib')
    library(nglShiny, lib.loc='/dls/science/groups/i04-1/software/xchemreview/xcrlib')
}

# Plotting Libs
library(ggplot2)
library(plotly)

sessionInfo()

# I can't believe this doesn't exist in R!
# Functional equivalent to string.rsplit('_', 1)
rsplit <- function(string, split_by, n=1){spl <- strsplit(string, split_by)[[1]]; c(paste(unlist(head(spl, length(spl)-n)), collapse=split_by), unlist(tail(spl, n)))}

getReviewData <- function(db, host_db, db_port, db_user, db_password){
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
        ligand_refinement_data[as.character(ligand_data$crystal_id),3:14],
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

    # Fix to handle duplicate row names... Use the latest modification date...
    #rids <- 1:nrow(output)
    #multi <- names(which(table(as.character(output$ligand_name)) > 1))
    #rids2 <- which(output$ligand_name %in% multi)
    #dupes <- cbind('id' =rids2, output[rids2, c('ligand_name','modification_date')])
    #todel <- unlist(sapply(unique(as.character(dupes$ligand_name)), function(x) {
    #    y <- dupes[dupes[,'ligand_name'] %in% x,]
    #    y[-which.max(as.numeric(y[,'modification_date'])),1]
    #}))
    #output <- output[-todel,]
    #rownames(output) <- output$ligand_name
    rownames(output) <- make.names(as.character(output$ligand_name), unique=TRUE)
    dbDisconnect(con)

    return(output)
}

getFragalysisViewData <- function(db, host_db, db_port, db_user, db_password){
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
    # Now we work on this exclusively...
    # This means we only annotate data which has been reviewed
    # AND has been marked as RELEASE...
    # This will also hide data that then becomes rejected down the line.
    annotatable_fv_dat <- fvdat[!fvdat$id %in% liganded_ligands[!ind, 1], ]
    targets <- dbGetQuery(con, sprintf("SELECT * from \"FragalysisTarget\" WHERE id IN (%s)", paste(unique(annotatable_fv_dat$fragalysis_target_id), collapse=',')))
    rownames(targets) <- as.character(targets$id)
    output <- cbind(annotatable_fv_dat, targetname=targets[as.character(annotatable_fv_dat$fragalysis_target_id), 'target'], md[as.character(annotatable_fv_dat$id),])

    #rids <- 1:nrow(output)
    #multi <- names(which(table(as.character(output$ligand_name))>1))
    #rids2 <- which(output$ligand_name %in% multi)
   #dupes <- cbind('id' =rids2, output[rids2, c('ligand_name','modification_date')])
    #todel <- unlist(sapply(unique(as.character(dupes$ligand_name)), function(x) {
    #    y <- dupes[dupes[,'ligand_name'] %in% x,]
    ##    y[-which.max(as.numeric(y[,'modification_date'])),1]
    #}))
    #output <- output[-todel,]
    dbDisconnect(con)
    #rownames(output) <- as.character(output$ligand_name)
    rownames(output) <- make.names(as.character(output$ligand_name), unique=TRUE)
    return(output)
}

updateOrCreateRow <- function(ligand_name_id, fragalysis_name, original_name, site_label='', new_smiles='', alternate_name='', pdb_id='',
    dbname, host, port, user, password){
    df = data.frame(Site_Label=as.character(site_label),
                    new_smiles=as.character(new_smiles),
                    alternate_name=as.character(alternate_name),
                    pdb_id=as.character(pdb_id),
                    fragalysis_name=as.character(fragalysis_name),
                    original_name=as.character(original_name),
                    Ligand_name_id=as.character(ligand_name_id))
    print(df)
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)
    id = dbGetQuery(con, sprintf("SELECT id from \"MetaData\" WHERE \"Ligand_name_id\"=%s", ligand_name_id))[1,1]
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
    draggableModalDialog(
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
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), data[[2]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, data[[2]][2], globalcolour)
    rasterImage(legend_image,.15,.6,.9,.7)

    text(x= 0.05, y = .5, labels = 'R Cryst')
    text(x= 0.95, y = .5, labels = data[[3]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), data[[3]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, data[[3]][2], globalcolour)
    rasterImage(legend_image,.15,.45,.9,.55)

    text(x= 0.05, y = .35, labels = 'Ram Outliers')
    text(x= 0.95, y = .35, labels = data[[4]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), data[[4]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, data[[4]][2], globalcolour)
    rasterImage(legend_image,.15,.3,.9,.4)

    text(x= 0.06, y = .2, labels = 'RMSD Angles')
    text(x= 0.95, y = .2, labels = data[[5]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), data[[5]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, data[[5]][2], globalcolour)
    rasterImage(legend_image,.15,0.15,.9,.25)

    text(x= 0.06, y = .05, labels = 'RMSD Bonds')
    text(x= 0.95, y = .05, labels = data[[6]][3])
    legend_image <- changeRasterRank(as.raster(rbind(colfunc(100))), data[[6]][1], experimentcolour )
    legend_image <- changeRasterRank(legend_image, data[[6]][2], globalcolour)
    rasterImage(legend_image,.15,0,.9,.1)

    text(x=0.15, y=-0.05, 'Worst')
    text(x=0.9, y=-0.05, 'Best')

    legend(x=0.10, y=-0.05, legend=c(sprintf('Percentile Relative to other %s crystals marked as comp. chem ready (or above)', target_name), 'Percentile Relative to crystals marked as comp. chem ready (or above) from all XChem Experiments '), fill=c(experimentcolour , globalcolour), bty='n')
}


header <- dashboardHeader(
    # What we see in the top bar
    title = 'XChemReview',
    # Commented out as idk what to do with these yet...
    #dropdownMenuOutput('messageMenu'),
    #dropdownMenuOutput('notifications'),
    #dropdownMenuOutput('tasks'),
    # The little cog in the top right for the controls.
    tags$li(class='dropdown', actionButton('controls', '', class = 'btn-primary', icon = icon('cog', lib = 'glyphicon')))
)

sidebar <- dashboardSidebar(
    sidebarMenu(
        id = 'tab',
        menuItem('Summary', tabName = 'summary', icon=icon('th'), badgeLabel = 'new', badgeColor = 'green'),
        menuItem('Review', tabName = 'review', icon = icon('dashboard')),
        menuItem('FragView', tabName = 'fragview', icon = icon('dashboard')),
        menuItem('LaunchPad', tabName = 'launchpad', icon = icon('th'), badgeLabel = 'new', badgeColor = 'green'),
        menuItem('Help', tabName = 'help', icon = icon('th')),
        hr(),
        # Flexible Sidebar options depending on which menuitem is selected.
        uiOutput('flex1')
    )
)

body <- dashboardBody(
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
                            actionButton("fitButton", "Center on Ligand"),
                            fluidRow(
                                chooseSliderSkin("Flat", color='#112446'),
                                column(6,
                                    fluidRow(
                                        column(2, checkboxInput('eventMap', 'Show Event Map', value = TRUE)),
                                        column(10, sliderInput("isoEvent", "", min = 0, max = 3, value = 1, step = 0.1))
                                        #column(10, uiOutput('isoEventSlider'))
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('twofofcMap', 'Show 2fofc Map', value = TRUE)),
                                        column(10, sliderInput("iso2fofc", "", min = 0, max = 3, value = 1.5, step = 0.1))
                                        #column(10, uiOutput('iso2fofcSlider'))
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('fofcMap', 'Show fofc Map', value = TRUE)),
                                        column(10, sliderInput("isofofc", "", min = 0, max = 3, value = 3, step = 0.1))
                                        #column(10, uiOutput('isofofcSlider'))
                                    )
                                ),
                                column(6,
                                    imageOutput('ligimage2')
                                )
                            )
                        ),
                        tabPanel(
                            title = 'Ligand Information',
                            div(style='overflow-y:scroll;height:600px;',
                            fluidRow(
                                column(8,
                                    column(6, imageOutput('ligimage')),
                                    column(6,imageOutput('spiderPlot'))
                                ),
                                column(4,
                                    div(style = "margin-top:-1em", checkboxInput('renderMisc', 'Render Image/Spider Plot', value = TRUE, width = NULL)),
                                    div(style = "margin-top:-1em", selectInput('emap', 'Select Eventmap', choices='', multiple=FALSE)),
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
                            title = 'Atom Selection',
                            DT::dataTableOutput('atoms')
                        ),
                        tabPanel(title = 'Slack',
                            fluidPage(
                                tags$head(
                                    tags$style("#chatpanel {overflow: auto;}")
                                ),
                                sidebarLayout(
                                    sidebarPanel(
                                        actionButton('updateSlackChannels', label = 'Update All Slack Channels'),
                                        selectizeInput("channelSelect", "", select='', choices = '', multiple=FALSE, width=-100),
                                        textAreaInput('TextInput', 'Message Body', value = "", width = NULL, height = NULL,
                                        cols = NULL, rows = NULL, placeholder = NULL, resize = 'both'),
                                        textInput('slackUser', label = 'Name', value =''),
                                        actionButton('slackSubmit', label = 'Submit')
                                    ), # sidebarpanel
                                    mainPanel(
                                            textOutput('chatURL'),
                                            textOutput('scrollDialog'),
                                            textOutput('chat')
                                    )
                                )
                            )
                        )
                    )
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
                    div(style='overflow-y:scroll;height:600px;',DT::dataTableOutput('therow'))))
            )
        ),
        tabItem(
            tabName = 'launchpad',
            h2('LaunchPad - Coming Soon...')
        )
    )
)

ui <- dashboardPage(header, sidebar, body)

server <- function(input, output, session){
    sID <- sample(1:100000, 1)
    debug = TRUE
    if(debug) debugMessage(sID=sID, sprintf('Session init'))
    session$allowReconnect(FALSE)
    sessionDisconnect <- function() debugMessage(sID=sID, 'Disconnected')
    session$onSessionEnded(sessionDisconnect)
    epochTime <- function() as.integer(Sys.time())
    humanTime <- function() format(Sys.time(), "%Y%m%d%H%M%OS")
    sessionTime <- reactive({epochTime()})


    slackFilters <- c('has joined the channel')
    getChannelList <- function(){
        #return(c('No Channels'))
        channellist <- list()
        channellist$response_metadata$next_cursor <- 'First'
        channels <- data.frame(name = 'welcome', id='abc', stringsAsFactors=F)
        while(!channellist$response_metadata$next_cursor == ''){
            if(channellist$response_metadata$next_cursor == 'First'){
                channellist <- httr::content(httr::POST(url='https://slack.com/api/conversations.list', body=list(token=api, limit=1000)))
            } else {
                channellist <- httr::content(httr::POST(url='https://slack.com/api/conversations.list',
                                                        body=list(token=api,
                                                                    limit=1000,
                                                                    cursor=channellist$response_metadata$next_cursor)))
            }
            data_to_join <- sapply(c('name', 'id'), function(x) sapply(channellist[[2]], '[[', x))
            channels <- rbind(channels, data_to_join)
        }
        channels <- channels[!channels[,1] %in% c('welcome', 'team', 'project', 'i04-1'),]
        rownames(channels) <- channels[,1]
        return(channels)
    }

    getChannel <- function(structure, channels){
        channelname <- tolower(gsub('[^[:alnum:]]', '', structure))
        ifelse(channelname %in% rownames(channels), channels[channelname, 2], NA)
    }

    parseConversation <- function(channel){
        history <- httr::content(httr::POST(url='https://slack.com/api/conversations.history',
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
        thread <- httr::content(httr::POST(url='https://slack.com/api/conversations.replies',
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
        users <- httr::content(httr::POST(url='https://slack.com/api/users.list',
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
        resp <- httr::POST(url='https://slack.com/api/conversations.create', body=list(token=api, name =channelname))
        return(httr::content(resp)$channel$id)
    }

    sendMessageToSlack <- function(channel, message, name){
        httr::content(httr::POST(url='https://slack.com/api/chat.postMessage',
            body=list(token=apiuser,
                channel=channel,
                text= message,
                as_user='true',
                username=name)))
    }

    sendEmail <- function(structure, user, decision, reason, comments){
        if(debug) debugMessage(sID=sID, sprintf('Communicating to Slack...'))
        #channels <- getChannelList()
        #channelID <- getChannel(structure, channels)
        #if(is.na(channelID)){
            # Create Channel, then get ID immediately
        #    channelID <- createChannel(structure=structure)
        #}
        # Post comment...
        #sendMessageToSlack(channel=channelID,
        #                    message= sprintf('%s has been labelled as %s by %s for the following reason(s): %s.
#With these additional comments:
#%s', structure, decision, user, reason, comments),
#                            name = 'xchemreview-bot')
	channelID <- 'NA'

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
If you disagree with this decision please discuss at the in the slack channel: (https://xchemreview.slack.com/archives/%s) or change the outcome by submitting a new response.
This email was automatically sent by The XChem Review app
If you have trouble joining the slack channel please use this invitation link: https://join.slack.com/t/xchemreview/shared_invite/zt-fpocaf6e-JQp~U6rcbGrre33E~7~faw
If you believe you have been sent this message in error, please email tyler.gorrie-stone@diamond.ac.uk',
            structure, decision, user, reason, comments, structure, protein, channelID), # replace NA with channelID
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

    # Loading Data Gubbins:
    restartSessionKeepOptions <- function(){
        message('Updating Data')
        dbdat <- getReviewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password)
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

    updateMainTable <- function(r1, pl=10){
        DT::renderDataTable({
            DT::datatable(
                r1(),
                selection = 'single',
                options = list(
                    pageLength = pl
                )
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

    updateMainTable2 <- function(r1, pl=10){
        DT::renderDataTable({
            DT::datatable(
                r1(),
                selection = 'single',
                options = list(
                    pageLength = pl
                )
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
    review_data <- getReviewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password)

    updateSelectInput(session, 'protein', selected = '', choices=c('', sort(unique(as.character(review_data$target_name)))))
    updateSelectInput(session, 'fpe_target', selected = '', choices=c('', sort(unique(as.character(review_data$target_name)))))

    inputData <- restartSessionKeepOptions()
    r1 <- reactiviseData(inputData=inputData, input=input)
    output$reviewtable <- updateMainTable(r1=r1)
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
                    data()[toshow, ]
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
                    data()[toshow, c('original_name', 'fragalysis_name', 'alternate_name','Site_Label', 'new_smiles', 'pdb_id')]
                } else {
                    data()[NULL,]
                }
            }
        })
    }

    reactivegetFragalysisViewData <- function(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password){
        reactive({getFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password)})
    }

    fvd <- getFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password)
    fragview_data <- reactivegetFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password)
    fragfolders <- c('', sort(unique(fvd$targetname)))
    print('Print FragFolders?')
    print(fragfolders)
    updateSelectInput(session, 'fragSelect', selected='', choices=fragfolders)

    fragview_input <- react_fv_data(fragview_data, input)
    fragview_table_data <- react_fv_data2(fragview_data, input)

    output$therow <- updateMainTable2(fragview_input, pl=100)


    fv_values <- reactiveValues()
    fv_values$apofiles <- c()
    fv_values$molfiles <- c()
    fv_values$molfil <- c()

    observeEvent(input$fragSelect,{
        if(debug) debugMessage(sID=sID, sprintf('Selecting: %s', input$fragSelect))
        #folderPath <- getAlignedStagingFolder()
        fv_values$apofiles <- as.character(isolate(fragview_table_data()$apo_pdb))
        fv_values$molfiles <- as.character(isolate(fragview_table_data()$lig_mol_file))
        fv_values$molfil <- gsub('.mol', '', basename(fv_values$molfiles))
        updateSelectInput(session, 'goto', choices = fv_values$molfil)
        fragview_input <- react_fv_data(fragview_data, input)
        output$therow <- updateMainTable2(fragview_input, pl=100)
        tryAddPDB <- try(uploadApoPDB(filepath=fv_values$apofiles[1], repr='cartoon'), silent=T)
        molout <- try(sapply(fv_values$molfiles, uploadUnfocussedMol), silent=T)

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
        folder <- dirname(molfiles[input$goto])
        mol_file <- molfiles[input$goto]
        smi_file <- gsub('.mol', '_smiles.txt', mol_file)
        smilestr <- system(sprintf('cat %s', smi_file), intern=T)
        # Fill Form as seen
        updateTextInput(session, 'crysname', value = input$goto)
        updateTextInput(session, 'smiles', value = smilestr)
        updateTextInput(session, 'new_smiles', value = as.character(isolate(fragview_table_data()[input$goto, 'new_smiles'])))
        updateTextInput(session, 'alternate_name', value = as.character(isolate(fragview_table_data()[input$goto, 'alternate_name'])))
        updateSelectizeInput(session, 'site_name', selected = as.character(isolate(fragview_table_data()[input$goto, 'Site_Label'])))
        updateTextInput(session, 'pdb_entry', value = as.character(isolate(fragview_table_data()[input$goto, 'pdb_id'])))
        # Go to specific ligand do not edit go next loop
        if(debug) debugMessage(sID=sID, sprintf('Selected: %s', input$goto))
        if(debug) debugMessage(sID=sID, sprintf('trying to view: %s', molfiles[input$goto]))
        gogogo <- try(uploadMolAndFocus2(mol_file), silent=T)
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
        fragview_data <- reactivegetFragalysisViewData(db=db, host_db=host_db, db_port=db_port, db_user=db_user, db_password=db_password)
        fragview_input <- react_fv_data(fragview_data, input)
        fragview_table_data <- react_fv_data2(fragview_data, input)
        output$therow <- updateMainTable2(fragview_input, pl=100)
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
            output$therow <- updateMainTable2(fragview_input, pl=100)
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
        output$reviewtable <- updateMainTable(r1=r1)
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
                    output$reviewtable <- updateMainTable(r1=r1)
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

    output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data)})

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
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))))})
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
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data, editable = list(target = 'cell', disable = list(columns = c(1,2))))})
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
                    textInput('comments', 'Additonal  Comments', value = "", width = NULL,placeholder = NULL),
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
            fogging = c(45,58),
            clipping = c(47,100),
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
    observeEvent(input$controls, ignoreNULL = FALSE, {
        if(is.null(input$controls)){
            title = 'As part of setup please confirm NGL Viewer Controls'
        } else {

        }
        showModal(
            controlPanelModal(
                values = isolate(ngl_control_values$defaults),
                title = 'NGL Viewer Controls'
            )
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

    uploadPDB <- function(filepath){
        syscall <- sprintf('cat %s', filepath)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse = '\n')
        session$sendCustomMessage(
            type = 'setPDB2', # See TJGorrie/NGLShiny for details on setPDB2
            message = list(
                choice,
                input$clipDist,
                input$clipping[1],
                input$clipping[2],
                input$fogging[1],
                input$fogging[2]
            )
        )
    }

    uploadApoPDB <- function(filepath, repr){
        syscall <- sprintf('cat %s', filepath)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse = '\n')
        session$sendCustomMessage(
            type = 'setapoPDB',
            message = list(
                choice,
                repr
            )
        )
    }

    uploadMolAndFocus <- function(filepath, ext){
        syscall <- sprintf('cat %s', filepath)
        pdbstrings <- system(syscall, intern = TRUE)
        choice <- paste0(pdbstrings, collapse = '\n')
        session$sendCustomMessage(
            type = 'addMolandfocus',
            list(choice,ext)
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

    tcl <- function(x) tolower(as.character(as.logical(x)))

    getExt <- function(x) sapply(strsplit(x, '[.]'), tail, 1)

    # Map Uploader
    uploadVolumeDensity <- function(filepath, color, negateiso = FALSE, boxsize, isolevel, visable, windowname){
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
                as.character(windowname)
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
        try(uploadMolAndFocus(isolate(sessionlist$mol_file), 'mol'), silent=T)
    })

    output$isoEventSlider <- renderUI({
            #chooseSliderSkin("Modern")
            sliderInput("isoEvent", "",
                    min = 0, max = 3,
                    value = 1, step = 0.1)
    })

    output$iso2fofcSlider <- renderUI({
            #chooseSliderSkin("Modern")
            sliderInput("iso2fofc", "",
                    min = 0, max = 3,
                    value = 1.5, step = 0.1)
    })

    output$isofofcSlider <- renderUI({
            #chooseSliderSkin("Modern")
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

    observeEvent(input$ligand, ignoreNULL = TRUE, {
        atomstoquery$data <- data.frame(name=character(),
                 index=character(),
                 comment=character(),
                 stringsAsFactors=FALSE)
        output$atoms <- DT::renderDataTable({DT::datatable(atomstoquery$data)})

        session$sendCustomMessage(type = 'setup', message = list())
        updateParam('mousePreset', as.character(input$mousePreset))

        updateSelectizeInput(session, 'channelSelect', select=tolower(gsub('[^[:alnum:]]', '', input$ligand)))
        withProgress(message = 'Loading Crystal', value = 0,{
            message('Crystal Updated')
            print(isolate(sessionlist$apo_file))
            if(! isolate(sessionlist$apo_file) == ""){
            try(uploadApoPDB(isolate(sessionlist$apo_file), 'ball+stick'), silent=T)
            incProgress(.1, detail = 'Uploading PDB file')
            message('Uploaded APO')
            try(uploadMolAndFocus(isolate(sessionlist$mol_file), 'mol'), silent=T)
            incProgress(.1, detail = 'Uploading mol file')
            message('Upload MOL')
            emaps <- dir(dirname(isolate(sessionlist$apo_file)), pattern='event', full=TRUE)
            names(emaps) <- basename(emaps)
            sessionlist$current_emaps <- emaps
            print(emaps)
            incProgress(.2, detail = 'Uploading Event map')
            updateSelectInput(session, 'emap', choices = names(isolate(sessionlist$current_emaps)), selected = names(isolate(sessionlist$current_emaps))[1])
            # Move this to a different part?
            message('Upload fofcs')
            incProgress(.2, detail = 'Uploading 2fofc map')
            try(uploadVolumeDensity(isolate(sessionlist$twofofc_file),
                color = 'blue', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$iso2fofc, visable=input$twofofcMap, windowname='twofofc'), silent=T)
            incProgress(.1, detail = 'Uploading fofc map')
            try(uploadVolumeDensity(isolate(sessionlist$fofc_file),
                color = 'lightgreen', negateiso = FALSE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcpos'), silent=T)
            incProgress(.1, detail = 'Uploading fofc map')
            try(uploadVolumeDensity(isolate(sessionlist$fofc_file),
                color = 'tomato', negateiso = TRUE, boxsize = input$boxsize, isolevel = input$isofofc, visable=input$fofcMap, windowname='fofcneg'), silent=T)
            }

            if(input$renderMisc){
                incProgress(.1, detail = 'Finding Misc Files...')
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
            setProgress(1)
        })
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

    observeEvent(input$updateSlackChannels,{
        channels <- getChannelList()
        message('ChannelList')
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
        message(input$channelSelect)
        channels <- getChannelList()
        channelSelect <- as.character(channels[input$channelSelect,2])
        sendMessageToSlack(channel = channelSelect,
                            message = sprintf('on behalf of: %s. \n %s', input$slackUser, input$TextInput),
                            name=input$slackUser)
        }
        refreshChat(channel = channelSelect[input$channelSelect])
    })

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
