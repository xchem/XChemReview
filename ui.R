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
                        selectizeInput("reason", "Reason(s)", list(), multiple=TRUE),
                        textInput('comments', 'Additonal  Comments', value = "", width = NULL,placeholder = NULL),
                        textOutput('msg'),
                        actionButton("submit", "Submit", class = "btn-primary"),
                        actionButton('clear', 'Clear', class = 'btn-primary'),
                        selectInput('columns', 'Select Columns to View? (delete/add more values as needed)', choices=colss, selected= defOrder, multiple = TRUE)
                    ),
                    fluidRow(
                        column(4, checkboxInput('out4', 'Comp Chem Ready', value = TRUE)),
                        column(4, checkboxInput('out5', 'Deposition Ready', value = FALSE)),
                        column(4, checkboxInput('out6', 'Deposited', value = FALSE))
                    ),
                    actionButton('pictureModal', 'Show Images')
                    #imageOutput('ligimage'),
                    #imageOutput('spiderPlot')
                    
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
                                actionButton("fitButton", "Center on Ligand"),
                                actionButton("defaultViewButton", "Restart Viewer"),
                                textOutput('selAtoms'), 
                                selectInput('assembly2', 'Assembly', choices=c('AU', 'UNITCELL', 'SUPERCELL'), selected='AU', multiple=FALSE),       
                                checkboxInput('eventMap', 'Event map', value = TRUE),
                                uiOutput('isoEventSlider'), 
                                checkboxInput('twofofcMap', '2fofc map', value = FALSE),
                                uiOutput('iso2fofcSlider'),
                                checkboxInput('fofcMap', 'fofc Map', value = FALSE),
                                uiOutput('isofofcSlider'), 
                                actionButton('controlPanel', 'Show Control Panel')#,
                                #fluidRow(
                                #    column(6, numericInput("boxsize", 'Box Size', value = 10, min = 0, max = 100, width='100px')),
                                #    column(6, numericInput("clipDist", "Clip Dist", value=5, min = 0, max = 100, width='100px'))
                                #),
                                #sliderInput("fogging", "Fogging:",
                                #    min = 0, max = 100,
                                #    value = c(45,58)
                                #),
                                #sliderInput("clipping", "Clipping:",
                                #    min = 0, max = 100,
                                #    value = c(47,100)
                                #)       
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
            fluidPage(
                tags$head(
                    tags$style("#panel3 {position: fixed;}"),
                    tags$style("#panel4 {overflow: auto;}")
                ),
                sidebarLayout(
                    sidebarPanel( width = 2,
                        selectInput('fragSelect', 'Project Select', choices=c('Select', fragfolders)),
                        checkboxInput('desync', 'Turn off automatic Updates', value = FALSE),
                        actionButton('goback', 'Prev Ligand'),
                        actionButton('gonext', 'Next Ligand'),
                        selectInput('goto', 'Go to Ligand', choices=list()),
                        textInput('crysname', 'Crystal Name', '' ),
                        textInput('smiles', 'Smiles String', ''),
                        textInput('new_smiles', 'New Smiles String', ''),
                        textInput('alternate_name', 'Alternate Fragment Name', ''),
                        selectizeInput('site_name', 'Site Label', list(), multiple=FALSE, options=list(create=TRUE)),
                        textInput('pdb_entry', 'PDB Entry', ''),
                        textOutput('metastatus'),
                        uiOutput('writeButton'),
                        hr(),
                        textInput('newCrystalName', 'New Fragment Name', ''),
                        actionButton('changeName', 'Change Name of Fragment (will not assign site label)'),             
                        hr(),
                        textOutput('massChange'),
                        selectizeInput('site_name2', 'Old label', list(), multiple=FALSE),
                        textInput('new_label', 'New label', ''),
                        actionButton('mcl', 'Mass Convert Label'),
                        hr(),
                        actionButton('updateTable', 'Refresh Metadata Table')
                    ), # sidebarpanel
                    mainPanel(
                        absolutePanel(id = 'panel3', top='6.5%', bottom='0%', width='90vw', height='50vw',fixed=T,
                            nglShinyOutput('FragViewnglShiny', height='600px')
                        ),
                        absolutePanel(id = 'panel4', top='70%', bottom='0%', fixed=T,
                            DT::dataTableOutput('therow')
                        )
                    ) # Main Panel
                ) # sidebar layout
            ) # Fluid Page
        ) # mainPanel
    ), # tabPanel
    tabPanel('FragChat',
        fluidPage(
                tags$head(
                    tags$style("#chatpanel {overflow: auto;}")
                ),
                sidebarLayout(
                    sidebarPanel(
                        actionButton('updateSlackChannels', label = 'Update All Slack Channels'),
                        selectizeInput("channelSelect", "Channel/Crystal", select='', choices = names(channelSelect), multiple=FALSE),
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
                    ) # Main Panel
                ) # sidebar layout
            ) # Fluid Page
    ), # tabPanel
	tabPanel('Help',
		includeMarkdown(sprintf('%s/%s', gpath, "Pages/include.md"))
	) # Tab Panel
) # Nav Bar Page
# End of UI
