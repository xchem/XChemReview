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
    easyClose=FALSE, footer = tagList(actionButton("updateParams", "Update Controls"))
  )
}