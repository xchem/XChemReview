header <- dashboardHeader(
    title = 'XChemReview',
    tags$li(class='dropdown', actionButton('controls', 'Additional NGL Controls', class = 'btn-primary', icon = icon('cog', lib = 'glyphicon')))
)

sidebar <- dashboardSidebar(
    sidebarMenu(
        id = 'tab',
        menuItem('Summary', tabName = 'summary', icon=icon('th')),
        menuItem('FragView', tabName = 'fragview', icon = icon('dashboard'), badgeLabel = 'Stage 1'),
        menuItem('Atom Quality Zone', tabName = 'aqz', icon = icon('dashboard'), badgeLabel = 'Stage 2'),
        menuItem('Review', tabName = 'review', icon = icon('dashboard'), badgeLabel = 'Stage 3'),
        menuItem('LaunchPad', tabName = 'launchpad', icon = icon('th'), badgeLabel = 'Stage 4'),
        menuItem('Help', tabName = 'help', icon = icon('th')),
        hr(),
        # Flexible Sidebar options depending on which menuitem is selected.
        uiOutput('flex1')
    )
)

# Copied from SO
customDraggableModalDialog <- function(..., title = NULL,
                                 footer = shiny::modalButton("Dismiss"),
                                 size = c("m", "s", "l"),
                                 easyClose = FALSE, fade = FALSE) {
  size <- match.arg(size)
  cls <- if (fade) { "modal fade" } else { "modal" }
  shiny::div(
    id = "shiny-modal",
    class = cls,
    # tabindex = "-1", This line should be commented out or removed
    #`data-backdrop` = if (!easyClose) { "static" } ,
    #`data-keyboard` = if (!easyClose) { "false" } ,
    shiny::div(
      class = "modal-dialog",
      class = switch(size, s = "modal-sm", m = NULL, l = "modal-lg"),
      jqui_draggable(shiny::div(
        class = "modal-content",
        if (!is.null(title)) {
          shiny::div(
            class = "modal-header",
            shiny::tags$h4(class = "modal-title",  title)
          )
        },
        shiny::div(class = "modal-body", ...),
        if (!is.null(footer)) {
          shiny::div(class = "modal-footer", footer)
        }
      ))
    ),
    shiny::tags$script("$('#shiny-modal').modal().focus();"),
    shiny::tags$style(HTML("
    .modal-backdrop{
      display: none;
    }
    .modal {
      pointer-events: none;
    }
    .modal-content {
      pointer-events: all;
    }"))
  )
}

body <- dashboardBody(
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
                            fluidRow(
                                column(6, actionButton(
                                    "fitButton",
                                    "Center on Ligand"
                                )),
                                column(6,checkboxInput('autocenter', 'Automatically Center on load', value=TRUE))
                            ),
                            fluidRow(
                                shinyWidgets::chooseSliderSkin("Flat", color='#112446'),
                                column(6,
                                    fluidRow(
                                        column(2, checkboxInput('eventMap', 'Show Event Map', value = TRUE)),
                                        column(10, sliderInput("isoEvent", "", min = 0, max = 3, value = 1, step = 0.1))
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('twofofcMap', 'Show 2fofc Map', value = TRUE)),
                                        column(10, sliderInput("iso2fofc", "", min = 0, max = 3, value = 1.5, step = 0.1))
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('fofcMap', 'Show fofc Map', value = TRUE)),
                                        column(10, sliderInput("isofofc", "", min = 0, max = 3, value = 3, step = 0.1))
                                    ),
                            	    selectInput('gotores', 'Go to Residue:', choices = '', multiple=FALSE),
                            	    selectizeInput('highlight_res', 'Highlight Residues:', choices = '', multiple=TRUE)
                                ),
                                column(6,
                                    helpText('Modelled Ligand'),
                                    imageOutput('ligimage2', height='300px'),
                                    radioButtons('views', 'View Type', selected = 'aligned', inline = FALSE, width = NULL,
                                        choiceNames = c('Aligned (what will be in Fragalysis)', 'Unaligned (to check if the api alignment introduces problems)', 'Raw Input Files (What you should see in coot, maps may take long time to load)'),
                                        choiceValues = c('aligned', 'unaligned', 'crystallographic')
                                    ),
                                    selectInput('asuSwitch', 'Assembly Type (Only in Raw and Unalign)', selected='AU', choices=c('AU', 'UNITCELL', 'SUPERCELL'))
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
                                        imageOutput('ligimage')
                                    ),
                                    column(6, 
                                        helpText('Modelled Ligand'),
                                        imageOutput('rlimage')
                                    )
                                ),
                                column(4,
                                    div(style = "margin-top:-1em", checkboxInput('renderMisc', 'Render Ligand Images', value = TRUE, width = NULL)),
                                    div(style = "margin-top:-1em", selectInput('emap', 'Select Eventmap', choices='', multiple=FALSE)),
                                    fluidRow(actionButton('buster', 'Buster Report'), materialSwitch(inputId = "bfactor", label = "Render B Factors", status='success', value = FALSE)),
                                    fluidRow(actionButton('interactions', 'Visualise Interactions'))
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
                            title = 'Ligand Relationships + Sites',
                            fluidRow(
                                tabBox(
                                    tabPanel(
                                        title = 'Relationships',
                                        fluidRow(column(12, DT::dataTableOutput('relationship_table')))
                                    ),
                                    tabPanel(
                                        title = 'Sites (WIP)',
                                        fluidRow(column(12, DT::dataTableOutput('site_table')))
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
                    div(style='overflow-y:scroll;height:600px;',DT::dataTableOutput('therow')), width=10)), textOutput('fv_warn'), textOutput('fv_hover')
            )
        ),
        tabItem(
            tabName = 'aqz',
            fluidRow(
                nglShinyOutput('AVnglShiny', height = '500px'),
            ),
            fluidRow(
                tabBox(width=4,
                    tabPanel(
                        title = 'Atom Selection',
                        htmlOutput('as_message'),
                        fluidRow(
                            column(3, fluidRow(actionButton('as_clear', label = 'Clear Atoms'), 
                            actionButton('submit_atoms', label='Submit Atom Qualities'))),
                            column(3, fluidRow(
                            actionButton('write_all', 'Write to All Atoms?', value=FALSE),
                            actionButton('write_selected', label = 'Write to selected rows')
                            )),
                            column(3, 
                                selectizeInput('atom_text', 'Comment', choices=c('', 'Weak Density', 'No Density Evidence', 'Unexpected Atom', 'Multiple Conformations', 'High b-factor'), options=list(create=TRUE)),
                                materialSwitch(inputId = "aq_bfactor", label = "Render B Factors", status='success', value = FALSE)
                            )
                        ),
                        DT::dataTableOutput('atoms')
                    )
                ),
                tabBox(width=3, id='ctab',
                    tabPanel(
                        title = 'Controls',
                        fluidRow(
                                shinyWidgets::chooseSliderSkin("Flat", color='#112446'),
                                    actionButton(
                                        "aq_fitButton",
                                        "Center on Ligand"
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('aq_eventMap', 'Show Event Map', value = TRUE)),
                                        column(10, sliderInput("aq_isoEvent", "", min = 0, max = 3, value = 1, step = 0.1))
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('aq_twofofcMap', 'Show 2fofc Map', value = TRUE)),
                                        column(10, sliderInput("aq_iso2fofc", "", min = 0, max = 3, value = 1.5, step = 0.1))
                                    ),
                                    fluidRow(
                                        column(2, checkboxInput('aq_fofcMap', 'Show fofc Map', value = TRUE)),
                                        column(10, sliderInput("aq_isofofc", "", min = 0, max = 3, value = 3, step = 0.1))
                                    )
                        )   
                    ), tabPanel(title = actionButton('aq_buster', 'Buster Report'))
                ),
                tabBox(width = 5,
                    tabPanel(title='Ligands',
                    div(style='overflow-y:scroll;height:600px;', DT::dataTableOutput('aqp'))
                    )
                )
            )
        ),
        tabItem(
            tabName = 'launchpad',
            uiOutput('launchpad_stuff')
        )
    )
)

ui <- dashboardPage(header, sidebar, body)