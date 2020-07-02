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
                        textInput("name", "Name/FedID", ""),
                        uiOutput('xtalselect'),
                        selectInput("decision", "Decision", choices = possDec),
                        selectizeInput("reason", "Reason(s)", list(), multiple=TRUE, options= list(create=TRUE)),
                        textInput('comments', 'Additonal  Comments', value = "", width = NULL,placeholder = NULL),
                        textOutput('msg'),
                        actionButton("submit", "Submit", class = "btn-primary"),
                        actionButton('clear', 'Clear', class = 'btn-primary'),
                        selectInput('columns', 'Select Columns to View? (delete/add more values as needed)', choices=colss, selected= defOrder, multiple = TRUE)
                    ),
                    fluidRow(
                        column(4, checkboxInput('out4', 'Comp Chem Ready', value = TRUE)),
                        column(4, checkboxInput('out5', 'Deposition Ready', value = TRUE)),
                        column(4, checkboxInput('out6', 'Deposited', value = FALSE))
                    ),
                    imageOutput('spiderPlot'),
                    imageOutput('ligimage')
                ), #sidebarpanel
                mainPanel(
                    absolutePanel(id = 'panel1', top='6.5%', bottom='0%', width='90vw', height='50vw',fixed=T,
                        fluidRow(
                            column(8,
                                textOutput('progtext'),
                                textOutput('missingFiles'),
                                nglShinyOutput('nglShiny', height='600px')
                                ),
                            column(2,
                                textOutput('msg3'),
                                actionButton("fitButton", "Fit to Ligand"),
                                actionButton("defaultViewButton", "Restart Viewer"), 
                                selectInput('assembly2', 'Assembly', choices=c('AU', 'UNITCELL', 'SUPERCELL'), selected='AU', multiple=FALSE),       
                                checkboxInput('eventMap', 'Event map', value = TRUE),
                                uiOutput('isoEventSlider'), 
                                checkboxInput('twofofcMap', '2fofc map', value = FALSE),
                                uiOutput('iso2fofcSlider'),
                                checkboxInput('fofcMap', 'fofc Map', value = FALSE),
                                uiOutput('isofofcSlider'), 
                                fluidRow(column(6, numericInput("boxsize", 'Box Size', value = 10, min = 0, max = 100, width='100px')),
                                         column(6, numericInput("clipDist", "Clip Dist", value=5, min = 0, max = 100, width='100px'))
                                ),
                                sliderInput("fogging", "Fogging:",
                                    min = 0, max = 100,
                                    value = c(45,58)),
                                sliderInput("clipping", "Clipping:",
                                    min = 0, max = 100,
                                    value = c(47,100)
                                    )       
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
    tabPanel('FragView',
        tabPanel('Main',
            sidebarLayout(
                sidebarPanel(
                    actionButton("restartViewer", "Restart Viewer"),
                    selectInput('fragSelect', 'Project Select', choices=c('Select', fragfolders)),
                    actionButton('gonext', 'Next Ligand'),
                    actionButton('goback', 'Prev Ligand'),
                    selectizeInput('goto', 'Go to Ligand', list(), multiple=FALSE, options= list(create=FALSE)),
                    selectizeInput('sitelabel', 'Site Label (no commas)', list(), multiple=FALSE, options=list(create=TRUE))
                ), # sidebarpanel
                mainPanel(
                    nglShinyOutput('FragViewnglShiny', height='600px')
                ) # Main Panel
            ) # sidebar layout
        ) # mainPanel
    ), # tabPanel
    tabPanel('FragChat',
        tabPanel('Main',
            sidebarLayout(
                sidebarPanel(

                ), # sidebarpanel
                mainPanel('Coming Soon...'

                ) # Main Panel
            ) # sidebar layout
        ) # mainPanel
    ), # tabPanel
	tabPanel('Help',
		includeMarkdown(sprintf('%s/%s', gpath, "Pages/include.md"))
	) # Tab Panel
) # Nav Bar Page
# End of UI