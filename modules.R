# Not so much modules as they are modals that I frequently want to use since I cannot get the modules working this will work for now...
pictureModal <- function(){
  modalDialog(
    fluidRow(
      column(6,imageOutput('ligimage')),
      column(6,imageOutput('spiderPlot'))
    ), easyClose=TRUE
  )
}


contolPanelModal <- function(values){
  modalDialog(title='Sorry to bug you - Please update the controls',
    numericInput("boxsize", 'Box Size', value = values$boxsize, min = 0, max = 100, width='100px'),
    numericInput("clipDist", "Clip Dist", value = values$clipDist, min = 0, max = 100, width='100px'),
    sliderInput("fogging", "Fogging:",min = 0, max = 100, value = values$fogging),
    sliderInput("clipping", "Clipping:", min = 0, max = 100, value = values$clipping),
    selectInput('backgroundColor', 'Background Colour', selected=values$backgroundColor, choices=c('black', 'white')),
    selectInput('cameraType', 'Camera Type', selected = values$cameraType, choices=c('orthographic', 'perspective')),
    selectInput('mousePreset', 'Mouse Preset', selected = values$mousePreset, choices=c('default', 'pymol', 'coot')),
    easyClose=FALSE, footer = tagList(actionButton("updateParams", "Update Controls"))
  )
}

slackPanel <- function(input){
  modalDialog(title='SlackSlackSlack',
            fluidPage(
                tags$head(
                    tags$style("#chatpanel {overflow: auto;}")
                ),
                sidebarLayout(
                    sidebarPanel(
                        actionButton('updateSlackChannels', label = 'Update All Slack Channels'),
                        selectizeInput("channelSelect", "", select='', choices = names(channelSelect), multiple=FALSE, width=-100),
                        textAreaInput('TextInput', 'Message Body', value = "", width = NULL, height = NULL,
                        cols = NULL, rows = NULL, placeholder = NULL, resize = 'both'),
                        textInput('slackUser', label = 'Name', value =''),
                        actionButton('slackSubmit', label = 'Submit')
                    ), # sidebarpanel
                    mainPanel(
                        absolutePanel(id = 'chatpanel', fixed=T,
                            #tableOutput('chatTable')]
                            textOutput('chatURL'),
                            textOutput('scrollDialog'),
                            textOutput('chat'),
                            tags$style(type="text/css", "#chat {white-space: pre-wrap; max-height: 500px}")
                        )
                    )
                )
            ),
  easyClose=TRUE, size='l')
}

