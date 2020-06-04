# UI Code
ui <- navbarPage("XChem Review", id='beep',
	# First Page
	tabPanel('Main',
		fluidRow(
			column(2,
                uiOutput('proteinselect'),
                div(
                    id = "form",
                    textInput("name", "FedID", ""),
                    uiOutput('xtalselect'),
                    #selectizeInput('Xtal', 'Which Structure?', choices = xtalList, multiple = FALSE),
                    selectInput("decision", "Decision", choices = possDec),
                    selectizeInput("reason", "Reason(s)", list(), multiple=TRUE, options(create=TRUE)),
                    textOutput('msg'),
                    actionButton("submit", "Submit", class = "btn-primary"),
                    actionButton('clear', 'Clear', class = 'btn-primary'),
                    #selectInput('protein', 'Select Specific Protein', choices = proteinList, selected= uiOutput("inVar"), multiple=TRUE),
                    selectInput('columns', 'Select Columns to View? (delete/add more values as needed)', choices=colss, selected= defOrder, multiple = TRUE)
                )
			),
			column(8,
				nglShinyOutput('nglShiny', height = '400px'),
                div(style="height: 50px;", br()),
            	DT::dataTableOutput("table")
			),
			column(2,
                textOutput('msg3'),
                #actionButton("fitButton", "Fit"),
                actionButton("defaultViewButton", "View/Update Parameters"),
                actionButton("clearRepresentationsButton", "Clear Representations"),
                #actionButton("updateView", "Update Parameters"),
                hr(),          
                checkboxInput('eventMap', 'Event map', value = FALSE),
                uiOutput('isoEventSlider'), 
                #sliderInput("isoEvent", "Event ISO",
                #    min = 0, max = 10,
                #    value = 1, step = 0.1),
                checkboxInput('twofofcMap', '2fofc map', value = FALSE),
                uiOutput('iso2fofcSlider'),
                #sliderInput("iso2fofc", "2fofc ISO",
                #    min = 0, max = 10,
                #    value = 1.5, step = 0.1),
                checkboxInput('fofcMap', 'fofc Map', value = FALSE),
                uiOutput('isofofcSlider'),
                #sliderInput("isofofc", "fofc ISO",
                #    min = 0, max = 10,
                #    value = 3, step = 0.1),
                hr(),
                numericInput("boxsize", 'Box Size', value = 10, min = 0, max = 100),
                numericInput("clipDist", "Clipping Distance", value=5, min = 0, max = 100),
                sliderInput("fogging", "Fogging:",
                  min = 0, max = 100,
                  value = c(45,58)),
                sliderInput("clipping", "Clipping:",
                  min = 0, max = 100,
                  value = c(47,100)),  
                hr()#,      
                #selectInput("representationSelector", "", nglRepresentations, selected=defaultRepresentation, width=0),
                #selectInput("colorSchemeSelector", "", nglColorSchemes, selected=defaultColorScheme,width=0)
			)
		) # Fluid row
	), # Tab Panel
	tabPanel('Help',
		includeMarkdown(sprintf('%s/%s', gpath, "Pages/include.md"))
	) # Tab Panel
) # Nav Bar Page
# End of UI