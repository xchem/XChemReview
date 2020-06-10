# UI Code
ui <- navbarPage("XChem Review", id='beep',
	# First Page
	tabPanel('Main',
		fluidPage(
            tags$head(
                tags$style("#panel1 {position: fixed;}"),
                tags$style("#panel2 {overflow: auto;}")
            ),
            sidebarLayout(
                sidebarPanel(width=2,
                    uiOutput('proteinselect'),
                    div(
                        id = "form",
                        textInput("name", "FedID", ""),
                        uiOutput('xtalselect'),
                        selectInput("decision", "Decision", choices = possDec),
                        selectizeInput("reason", "Reason(s)", list(), multiple=TRUE, options= list(create=TRUE)),
                        textInput('comments', 'Additonal  Comments', value = "", width = NULL,placeholder = NULL),
                        textOutput('msg'),
                        actionButton("submit", "Submit", class = "btn-primary"),
                        actionButton('clear', 'Clear', class = 'btn-primary'),
                        selectInput('columns', 'Select Columns to View? (delete/add more values as needed)', choices=colss, selected= defOrder, multiple = TRUE)
                    ),
                    imageOutput('spiderPlot'),
                    fluidRow(
                        column(4, checkboxInput('out4', 'Structures in Refinement', value = TRUE)),
                        column(4, checkboxInput('out5', 'CompChem Ready Structures', value = TRUE)),
                        column(4, checkboxInput('out6', 'Deposited Structures', value = FALSE))
                    ),
                    textOutput('missingFiles')
                    #uiOutput('spiderPlot')
                ), #sidebarpanel
                mainPanel(
                    absolutePanel(id = 'panel1', top='6.5%', bottom='0%', width='90vw', height='50vw',fixed=T,
                        fluidRow(
                            column(8,nglShinyOutput('nglShiny', height='600px')),
                            column(2,
                                textOutput('msg3'),
                                actionButton("fitButton", "Fit to Ligand"),
                                actionButton("defaultViewButton", "Restart Viewer"),
                                hr(),          
                                checkboxInput('eventMap', 'Event map', value = TRUE),
                                uiOutput('isoEventSlider'), 
                                checkboxInput('twofofcMap', '2fofc map', value = FALSE),
                                uiOutput('iso2fofcSlider'),
                                checkboxInput('fofcMap', 'fofc Map', value = FALSE),
                                uiOutput('isofofcSlider'),
                                hr(),
                                fluidRow(column(6, numericInput("boxsize", 'Box Size', value = 10, min = 0, max = 100, width='100px')),
                                        column(6, numericInput("clipDist", "Clip Dist", value=5, min = 0, max = 100, width='100px'))
                                ),
                                sliderInput("fogging", "Fogging:",
                                    min = 0, max = 100,
                                    value = c(45,58)),
                                sliderInput("clipping", "Clipping:",
                                    min = 0, max = 100,
                                    value = c(47,100)),  
                                hr()      
                            ) # column
                        ) # Fluid row
                    ), 
                    absolutePanel(id = 'panel2', top='70%', bottom='0%', fixed=T,
                        fluidRow(column(8, DT::dataTableOutput("table")), column(2))
                    )
                ), # main panel
            ) # sidebarlayout
        ) # Fluid Page
    ), # Tab Panel
	tabPanel('Help',
		includeMarkdown(sprintf('%s/%s', gpath, "Pages/include.md"))
	) # Tab Panel
) # Nav Bar Page
# End of UI