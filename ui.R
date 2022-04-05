# Main UI contains 3 elements, a header, a siderbar and a body.

# The header will remain static, there is a button called 'Additional NGL Controls' 
# that provide global options for the NGL viewer.

# The Sidebar contains navigation for the tabs within the application. As well as contain unique form inputs for 
# each tab depending what each tab requires in terms of functionality. This is specified in the flex1 UI element that 
# gets constructed in the server.R functionality. So if you see an UI element that is beneath all the tab options, you 
# need to look at a different piece of code.

# The Body elements content will change depending on which tab is selected. For reference
# = Summary is the splash page and will give a count of all the ligands they have access to review/annotate
# = FragView provides an NGL window to annotate the biological meaning of all ligands - classifying the site labels for
# fragalysis
# - Atom Quality Zone provides an NGL window to ascribe per atom quality to both ligand and protein atoms. It used to be
#   part of the review tab but Frank felt it was important enough to require it's own tab.
# - Review provides an NGL window to finally review and give the final checkover for what will appear in Fragalysis.
# - LaunchPad provides a simple UI to upload a set of data to Fragalysis
# - Help, provides documentation within the app to assist users.

# In theory, and with better organisation, this the FragView, AQZ and Review could all be handled in the same screen with
# some exception clever handling of data. For the time being it is organised this way as Frank expects the review to do 
# Stage 1 and 2 themselves, leaving the reviewer to do stage 3.

# Documentation for all shiny functionality: https://rdrr.io/cran/shiny/man/
# Documentation for all shinydashboard functionality: https://rdrr.io/cran/shinydashboard/man/

header <- dashboardHeader(
    title = 'XChemReview',
    tags$li(class='dropdown', actionButton('controls', 'Additional NGL Controls', class = 'btn-primary', icon = icon('cog', lib = 'glyphicon')))
)

sidebar <- dashboardSidebar(
    sidebarMenu(
        id = 'tab',
        # This creates an object input$tab in the server.R code, which will have the values indicated be `tabName`
        # If you need behaviour to occur when a tab is switched-to for switched-from use the existing observeEvent(input$tab, {...}) in server.R
        menuItem('Summary', tabName = 'summary', icon=icon('th')),
        menuItem('FragView', tabName = 'fragview', icon = icon('dashboard'), badgeLabel = 'Stage 1'),
        menuItem('Atom Quality Zone', tabName = 'aqz', icon = icon('dashboard'), badgeLabel = 'Stage 2'),
        menuItem('Review', tabName = 'review', icon = icon('dashboard'), badgeLabel = 'Stage 3'),
        menuItem('LaunchPad', tabName = 'launchpad', icon = icon('th'), badgeLabel = 'Stage 4'),
        menuItem('Help', tabName = 'help', icon = icon('th')),
        menuItem('Pipeline Config', tabName='config', icon=icon('th'))
        hr(),
        # This is the flexible UI element. Referenced in server.R using the input$tab variable incase elements
        # need to be changed later.
        uiOutput('flex1') # output$flex1
    )
)


body <- dashboardBody(
    # Initialize any addition js we may need.
    useWaiter(),
	tags$head(tags$script("$(function() {$.fn.dataTableExt.errMode = 'throw';});")),
    tags$head(shiny::tags$style(HTML("
    .modal-backdrop{
      display: none;
    }
    .modal {
      pointer-events: none;
    }
    .modal-content {
      pointer-events: all;
    }"))),
    # dashboard body can create tabs based on the tabNames described in the sidebar
    # All body elements are contained within the tabItems(...) call, 
    # with each new page body contained within a tabItem(...). One sidebar tab = one tabItem()
    tabItems(
        tabItem(
            tabName = 'review',
            # Just to describe some level of semantics for positioning things in shiny.
            # If you want a set of features to appear together with the same vertical positioning you wrap them in a fluidRow.
            # Calling two fluidRow(), in sequence will create two non-overlapping rows.
            # Each row in a fluidRow as a unitless with of 12, which can be paritioned into whole numbers using column(...),
            fluidRow( 
                nglShinyOutput('nglShiny', height = '500px'), # output$nglShiny
                jqui_draggable(
                    tabBox(
                        tabPanel(
                            title = 'NGL Controls',
                            fluidRow(
                                column(6, actionButton(
                                    "fitButton", # input$fitButton in server.R
                                    "Center on Ligand" # Send custom message to center on ligand in NGLviewer if moved away
                                )),
                                column(6, checkboxInput(
                                    'autocenter', # input$autocenter in server.R
                                    'Automatically Center on load', # Toggle center on ligand for crystals.
                                value=TRUE))
                            ),
                            fluidRow(
                                shinyWidgets::chooseSliderSkin("Flat", color='#112446'),
                                column(6,
                                    # 3 elements that contain a toggle and slider for electron density maps.
                                    fluidRow(
                                        column(2, checkboxInput('eventMap', 'Show Event Map', value = TRUE)), # input$eventMap
                                        column(10, sliderInput("isoEvent", "", min = 0, max = 3, value = 1, step = 0.1)) # input$isoEvent
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('twofofcMap', 'Show 2fofc Map', value = TRUE)), # input$twofofcMap
                                        column(10, sliderInput("iso2fofc", "", min = 0, max = 3, value = 1.5, step = 0.1)) # input$iso2fofc
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('fofcMap', 'Show fofc Map', value = TRUE)), # input$fofcMap
                                        column(10, sliderInput("isofofc", "", min = 0, max = 3, value = 3, step = 0.1)) # input$isofofc
                                    ),
                            	    selectInput('gotores', 'Go to Residue:', choices = '', multiple=FALSE), # input$gotores
                            	    selectizeInput('highlight_res', 'Highlight Residues:', choices = '', multiple=TRUE) # input$highlight_res
                                ),
                                column(6,
                                    helpText('Modelled Ligand'),
                                    # Render image as calculated by fragalysis-api */*_*.png
                                    imageOutput('ligimage2', height='300px'), # input$ligimage2
                                    # Toggles to flip between aligned, unaligned and crystallographic views.
                                    radioButtons('views', 'View Type', selected = 'aligned', inline = FALSE, width = NULL, # input$views
                                        choiceNames = c('Aligned (what will be in Fragalysis)', 'Unaligned (to check if the api alignment introduces problems)', 'Raw Input Files (What you should see in coot, maps may take long time to load)'),
                                        choiceValues = c('aligned', 'unaligned', 'crystallographic')
                                    ),
                                    # Toggle the Assembly (the symmetry in the files is mostly busted - but can be useful.)
                                    selectInput('asuSwitch', 'Assembly Type (Only in Raw and Unalign)', selected='AU', choices=c('AU', 'UNITCELL', 'SUPERCELL')) # input$asuSwitch
                                )
                            )
                        ),
                        tabPanel(
                            title = 'Ligand Information',
                            div(style='overflow-y:scroll;height:600px;',
                            fluidRow(
                                column(8,
                                    column(6, 
                                        helpText('Input Ligand'),
                                        imageOutput('ligimage') # input$$ligimage
                                    ),
                                    column(6, 
                                        helpText('Modelled Ligand'),
                                        imageOutput('rlimage') # input$rlimage
                                    )
                                ),
                                column(4,
                                    div(style = "margin-top:-1em", checkboxInput('renderMisc', 'Render Ligand Images', value = TRUE, width = NULL)), # input$renderMisc
                                    div(style = "margin-top:-1em", selectInput('emap', 'Select Eventmap', choices='', multiple=FALSE)), # input$emap
                                    fluidRow(                                                                             
                                        actionButton('buster', 'Buster Report'), # input$buster
                                        materialSwitch(inputId = "bfactor", label = "Render B Factors", status='success', value = FALSE) # input$bfactor
                                    ),
                                    fluidRow(actionButton('interactions', 'Visualise Interactions')) # input$interactions
                                )
                            ),
                            column(12,div(style = "margin-top:-15em",
                                fluidRow(
                                    uiOutput('plotElement') # output$plotElement
                                )
                            ))
                            ),
                        ),
                        tabPanel(
                            title = 'Ligand Relationships + Sites',
                            fluidRow(
                                tabBox(
                                    tabPanel(
                                        title = 'Relationships',
                                        fluidRow(column(12, DT::dataTableOutput('relationship_table'))) # output$relationship_table
                                    ),
                                    tabPanel(
                                        title = 'Sites (WIP)',
                                        fluidRow(column(12, DT::dataTableOutput('site_table'))) # output$site_table
                                    )
                                )
                            )
                        )    
                    ), options = list(delay = '1000', cancel = '.selectize-control')
                ),
                jqui_draggable(
                    tabBox(
                        tabPanel(
                            title='Review Table',
                            div(
                                style='overflow-y:scroll;height:600px;',
                                DT::DTOutput('reviewtable') # output$reviewtable
                            )
                        ),
                        tabPanel( # This section is largely redundant since no one uses it, but some people may use it. It controls the panelling that controls a plotly plot.
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
            tabName = 'config',
            h2('N.B. Fiddling with this will trigger changes to how data will look in XCR - please give 24 hours after changing parameters')
            selectInput('config_target', 'Target', choices=c('')),
            checkboxInput('monomeric', 'Do not run PISA/Gemmi convert on target?' value = FALSE),
            checkboxInput('reduce', 'Reduce Reference Frame to single Chain (usually Chain A)' value = TRUE),
            checkboxInput('covalent', 'Convert covalently attached mol files', value = TRUE),
            checkboxInput('active', 'Actively keep data up to date', value = TRUE),
            actionButton('config_change', 'Submit')
        ),
        tabItem(
            tabName = 'summary',
            # Summary Page needs more work.
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
                nglShinyOutput('FragViewnglShiny', height = '500px'), # output$FragViewnglShiny
                jqui_draggable(tabBox(
                    div(style='overflow-y:scroll;height:600px;',
                        DT::dataTableOutput('therow') # output$therow
                    ), width=10)), textOutput('fv_warn'), textOutput('fv_hover')
            )
        ),
        tabItem(
            tabName = 'aqz',
            fluidRow(
                nglShinyOutput('AVnglShiny', height = '500px'), # output$AVnglShiny
            ),
            fluidRow(
                tabBox(width=4,
                    tabPanel(
                        title = 'Atom Selection',
                        htmlOutput('as_message'),
                        fluidRow(
                            column(3, fluidRow(actionButton('as_clear', label = 'Clear Atoms'), # input$as_clear
                            actionButton('submit_atoms', label='Submit Atom Qualities'))), # input$submit_atoms
                            column(3, fluidRow(
                            actionButton('write_all', 'Write to All Atoms?', value=FALSE), # input$write_all
                            actionButton('write_selected', label = 'Write to selected rows') # input$write_selected
                            )),
                            column(3, 
                                selectizeInput('atom_text', 'Comment', choices=c('', 'Weak Density', 'No Density Evidence', 'Unexpected Atom', 'Multiple Conformations', 'High b-factor'), options=list(create=TRUE)), # input$atom_text
                                materialSwitch(inputId = "aq_bfactor", label = "Render B Factors", status='success', value = FALSE) # input$aq_bfactor
                            )
                        ),
                        DT::dataTableOutput('atoms') # output$atoms
                    )
                ),
                tabBox(width=3, id='ctab',
                    tabPanel(
                        title = 'Controls',
                        fluidRow(
                                shinyWidgets::chooseSliderSkin("Flat", color='#112446'),
                                    actionButton(
                                        "aq_fitButton", # input$aq_fitButton
                                        "Center on Ligand"
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('aq_eventMap', 'Show Event Map', value = TRUE)), # input$aq_eventMap
                                        column(10, sliderInput("aq_isoEvent", "", min = 0, max = 3, value = 1, step = 0.1)) # input$aq_isoEvent
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('aq_twofofcMap', 'Show 2fofc Map', value = TRUE)), # input$aq_twofofcMap
                                        column(10, sliderInput("aq_iso2fofc", "", min = 0, max = 3, value = 1.5, step = 0.1)) # input$aq_iso2fofc
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('aq_fofcMap', 'Show fofc Map', value = TRUE)), # input$aq_fofcMap
                                        column(10, sliderInput("aq_isofofc", "", min = 0, max = 3, value = 3, step = 0.1)) # input$isofofc
                                    )
                        )   
                    ), tabPanel(title = actionButton('aq_buster', 'Buster Report')) # input$aq_buster
                ),
                tabBox(width = 5,
                    tabPanel(title='Ligands',
                    div(style='overflow-y:scroll;height:600px;', DT::dataTableOutput('aqp')) # output$aqp
                    )
                )
            )
        ),
        tabItem(
            tabName = 'launchpad',
            uiOutput('launchpad_stuff') # output$launchpad_stuff - defined in server.R.
        )
    )
)

# Create the UI, used in app.R to launch the application.
ui <- dashboardPage(header, sidebar, body)