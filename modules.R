# Modal module UI
modalModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    textOutput(ns("myText")),
    actionButton(ns("openModalBtn"), "Open Modal")
  )
}

myModal <- function(session) {
  ns <- session$ns
  modalDialog(
    textInput(ns("modalTextInput"), "Show Me What You Got!", "Get Schwifty!"),
    actionButton(ns("closeModalBtn"), "Close Modal")
  )
}

# Modal module server
modalModule <- function(input, output, session) {
  
  # open modal on button click
  observeEvent(input$openModalBtn,
               ignoreNULL = TRUE,   # Show modal on start up
               showModal(myModal(session))
  )
  
  # show user text input in main window (not working)
  output$myText <- renderText({ req(input$modalTextInput) })
  
  # close modal on button click (not working)
  observeEvent(input$closeModalBtn, { 
    removeModal() 
  })
}

# Main app UI
ui <- fluidPage(modalModuleUI("controlPanel"))

# Main app server
server <- function(input, output, session) {
  callModule(modalModule, "foo")
}
