shinyUI(navbarPage("CellCycler",
tabPanel("Cells",
         #img(src = "SystemsForecastingBanner.png", height = '10%', width = '100%'),
         a(img(src = "SystemsForecastingBanner.png", height = '10%', width = '100%'),
           target="_blank",href="http://www.systemsforecasting.com"),
         sidebarLayout(
           sidebarPanel(
             helpText("Press button below to begin simulation. Plot shows volume of cells
                in each phase, with total normalised to 1. Use PK 
                tabs to set drug action, Tumor tab to plot tumor radius."),
             fluidRow(
               column(6,
                      h5('Simulation time')
               ),
               column(6,
                      numericInput("tmax", label = NULL, value = 100)
               )
             ),
             fluidRow(
               column(6,
                      h5('Doubling time')
               ),
               column(6,
                      numericInput("tdoub", label = NULL, value = 24)
               )
             ),
            fluidRow(
               column(6,
                      h5('G1 phase')
               ),
               column(6,
                      numericInput('G1phase', label = NULL, value = 0.2, step=0.1)
               )
             ),
             fluidRow(
               column(6,
                      h5('S phase')
              ),
              column(6,
                    numericInput('Sphase', label = NULL, value = 0.3, step=0.1)
              )
            ),
            fluidRow(
              column(6,
                     h5('M phase')
              ),
              column(6,
                    numericInput('Mphase', label = NULL, value = 0.1, step=0.1)
             )
            ),
            fluidRow(
             column(12,
                    h5(textOutput('G2phase'))
             )
           ),
           p(' '),
           actionButton("runButton", "Run simulation", width='250px'),
           hr(),
           downloadButton('savefile', label='Save settings to file'),
           hr(),
           fileInput("readfile", label = "Read settings from file")
           ),
           mainPanel(
             fluidRow(
               column(4,
                      h5("Plot options"),
                      checkboxGroupInput("plotVariables", 
                                         label = NULL, 
                                         choices = list("growing volumes" = 1, 
                                                        "damaged volumes" = 2),
                                         selected = 1)
                      # checkboxInput("phaseSum", label = "sum in each phase", value = TRUE)

               ),
               column(4,
                      h5("Axis options"),
                      checkboxInput("logAxis", label = "semilog", value = TRUE)
               ),
               column(4,
                      sliderInput("rangeAxis", "Time range",
                                  min = 0, max = 100, value = c(0,100), width = '200px')
                      )
             ),
             plotOutput("volPlot"),
             br()
           )
           )
         ),

tabPanel("PK1",
         a(img(src = "SystemsForecastingBanner.png", height = '10%', width = '100%'),
           target="_blank",href="http://www.systemsforecasting.com"),
         sidebarLayout(
           sidebarPanel(
             helpText('Run PK model to create time series. Drug effect parameters can be 
                      adjusted without recalculation.'),
             #              selectInput('pknum', label=NULL, c(1,2), selectize=FALSE),
             #              verbatimTextOutput('out1'),
             fluidRow(
               column(6,
                      h5('Type')
               ),
               column(6,
                      selectInput("pktype1", label = NULL, 
                                  choices=c('Step','K-PD','1-cmpt','2-cmpt','R file'), 
                                  selected='K-PD')
               )
             ),
             conditionalPanel(
               condition = "input.pktype1 != 'File'",
               fluidRow(
                 column(6,
                        h5('Dose/V')
                 ),
                 column(6,
                        numericInput('dose1', label = NULL, value = 0, step=0.1)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.pktype1 != 'File' & input.pktype1 != 'Step'",
               fluidRow(
                 column(6,
                        h5('Elimination rate')
                 ),
                 column(6,
                        numericInput('pkke1', label = NULL, value = 0.1, step=0.01)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.pktype1 == '1-cmpt' | input.pktype1 == '2-cmpt'",
               fluidRow(
                 column(6,
                        h5('Absorption rate')
                 ),
                 column(6,
                        numericInput('pkka1', label = NULL, value = 0.4, step=0.01)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.pktype1 == '2-cmpt'",
               fluidRow(
                 column(6,
                        h5('Peripheral rate')
                 ),
                 column(6,
                        numericInput('pkkp1', label = NULL, value = 0.4, step=0.01)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.pktype1 == 'File'",
               fileInput("sourcepkfile", label = NULL)
             ),
             p(' '),
             actionButton("runpkButton1", "Run PK 1 model", width='250px'),
             hr(),
             fluidRow(
               column(6,
                      h5('Phase')
               ),
               column(6,
                      selectInput("drugphase1", label = NULL, choices=c('All','G1','S','G2','M'), 
                                  selected='All')
               )
             ),
             fluidRow(
               column(6,
                      h5('Apoptosis rate')
               ),
               column(6,
                      numericInput("kap", label = NULL, value = 1, step=0.1)
               )
             ),
             fluidRow(
               column(6,
                      h5('Damage rate')
               ),
               column(6,
                      numericInput("kdmg", label = NULL, value = 1, step=0.1)
               )
             ),
             fluidRow(
               column(6,
                      h5('Repair rate')
               ),
               column(6,
                      numericInput("krep", label = NULL, value = 0.1, step=0.1)
               )
             )
           ),
           mainPanel(
             fluidRow(
               column(4,
                      p('')
                      #                       h5("Plot options"),
                      #                       checkboxGroupInput("pkplotVariables", 
                      #                                          label = NULL, 
                      #                                          choices = list("PK model 1" = 1, 
                      #                                                         "PK model 2" = 2),
                      #                                          selected = 1)
               ),
               column(4,
                      h5("Axis options"),
                      checkboxInput("pklogAxis", label = "semilog", value = FALSE)
               ),
               column(4,
                      sliderInput("pkrangeAxis", "Time range",
                                  min = 0, max = 100, value = c(0,100), width = '200px')
               )
             ),
             plotOutput("pkPlot"),
             br(),
             textInput('sch1', 'Schedule 1', '0,48', width='100%'),
             numericInput('wks1', label='Weeks repeated', value = 1, step=1, width='20%')
           )
         )
         
         
         
),

tabPanel("PK2",
         a(img(src = "SystemsForecastingBanner.png", height = '10%', width = '100%'),
           target="_blank",href="http://www.systemsforecasting.com"),
         sidebarLayout(
           sidebarPanel(
             helpText('Run PK model to create time series. Drug effect parameters can be 
                      adjusted without recalculation.'),
             fluidRow(
               column(6,
                      h5('Type')
               ),
               column(6,
                      selectInput("pktype2", label = NULL, 
                                  choices=c('Step','K-PD','1-cmpt','2-cmpt','R file'), 
                                  selected='K-PD')
               )
             ),
             conditionalPanel(
               condition = "input.pktype2 != 'File'",
               fluidRow(
                 column(6,
                        h5('Dose/V')
                 ),
                 column(6,
                        numericInput('dose2', label = NULL, value = 0, step=0.1)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.pktype2 != 'File' & input.pktype2 != 'Step'",
               fluidRow(
                 column(6,
                        h5('Elimination rate')
                 ),
                 column(6,
                        numericInput('pkke2', label = NULL, value = 0.1, step=0.01)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.pktype2 == '1-cmpt' | input.pktype2 == '2-cmpt'",
               fluidRow(
                 column(6,
                        h5('Absorption rate')
                 ),
                 column(6,
                        numericInput('pkka2', label = NULL, value = 0.4, step=0.01)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.pktype2 == '2-cmpt'",
               fluidRow(
                 column(6,
                        h5('Peripheral rate')
                 ),
                 column(6,
                        numericInput('pkkp2', label = NULL, value = 0.4, step=0.01)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.pktype2 == 'File'",
               fileInput("sourcepkfile2", label = NULL)
             ),
             p(' '),
             actionButton("runpkButton2", "Run PK 2 model", width='250px'),
             hr(),
             fluidRow(
               column(6,
                      h5('Phase')
               ),
               column(6,
                      selectInput("drugphase2", label = NULL, choices=c('All','G1','S','G2','M'), 
                                  selected='All')
               )
             ),
             fluidRow(
               column(6,
                      h5('Apoptosis rate')
               ),
               column(6,
                      numericInput("kap2", label = NULL, value = 1, step=0.1)
               )
             ),
             fluidRow(
               column(6,
                      h5('Damage rate')
               ),
               column(6,
                      numericInput("kdmg2", label = NULL, value = 1, step=0.1)
               )
             ),
             fluidRow(
#                column(6,
#                       h5('Repair rate')
#                ),
               column(12,
                      h5(textOutput('repairrate'))
                      #numericInput("krep2", label = NULL, value = 0.1, step=0.1)
               )
             )
           ),
           mainPanel(
             fluidRow(
               column(4,
                      p('')
               ),
               column(4,
                      h5("Axis options"),
                      checkboxInput("pklogAxis2", label = "semilog", value = FALSE)
               ),
               column(4,
                      sliderInput("pkrangeAxis2", "Time range",
                                  min = 0, max = 100, value = c(0,100), width = '200px')
               )
             ),
             plotOutput("pkPlot2"),
             br(),
             textInput('sch2', 'Schedule 2', '0', width='100%'),
             numericInput('wks2', label='Weeks repeated', value = 1, step=1, width='20%')
           )
           )
         

         ),

tabPanel("Tumor",
         a(img(src = "SystemsForecastingBanner.png", height = '10%', width = '100%'),
           target="_blank",href="http://www.systemsforecasting.com"),
         sidebarLayout(
           sidebarPanel(
             helpText("Plot shows growth in tumour radius for a given initial radius and
                      growing layer thickness. Select drug effect for radius lost due to cell death
                      from drugs 1 and 2 (red and blue) and damage (green). Can also overlay
                      experimental data filtered by treatment, default is data from ACC 
                      study, see Notes."),
             fluidRow(
               column(6,
                      h5('Initial radius')
               ),
               column(6,
                      numericInput("rini", NULL, value = 6.8, step=0.1)
               )
             ),
             fluidRow(
               column(6,
                      h5('Growing layer')
               ),
               column(6,
                      numericInput("dgr", label = NULL, value = 0.24, step=0.1)
               )
             ),
             downloadButton('downloadData', 'Save simulated radius'),
             hr(),
             fileInput("dataFile", label = h5("Read data for overlay from file")),
             selectInput("treatment",label='Treatment',choices='untreated',selected='untreated'),
             h5(tableOutput('resultsText')), 
             h5(tableOutput('resultsTextDrug1')), 
             h5(tableOutput('resultsTextDrug2')), 
             hr()
           ),
           mainPanel(
             fluidRow(
               column(4,
                      h5("Plot options"),
                      checkboxGroupInput("tumplotVariables", 
                                         label = NULL,
                                         choices = list('control' = 1,
                                                        'drug effect' = 2),
                                         selected = NULL)
               ),
               column(4,
                      h5("Overlay options"),
                      checkboxGroupInput("tumplotOverlay", 
                                         label = NULL, 
                                         choices = list("overlay", 
                                                        "linear fit"),
                                         selected = NULL)
               ),
               column(4,
                      sliderInput("radrangeAxis", "Time range",
                                  min = 0, max = 100, value = c(0,100), width = '500px')
               )
             ),
             plotOutput("radPlot"),
             fluidRow(
               column(6,
                      h5(tableOutput('fitcoefsText'))
               )
             )
           )
         )
),
tabPanel("Advanced",
         a(img(src = "SystemsForecastingBanner.png", height = '10%', width = '100%'),
           target="_blank",href="http://www.systemsforecasting.com"),
         shinyUI(fluidPage(
           h4('Number of compartments'),
           p('The CellCycler divides the cell cycle into N discrete compartments. As discussed
             in the Notes, the model simulates a population of cells whose doubling time
             is variable, with a normalised standard deviation equal to 1/sqrt(N). 
             So 25 compartments has a 
             normalised standard deviation of about 0.2 or 20%, while 100 compartments reduces 
             this to 0.1 or 10%. The number of compartments also controls the discretisation of 
             the cell cycle. For example with 50 compartments the proportional phase
             times, set in the Cells tab, will be rounded off internally to the nearest 1/50=0.02.
             Note that run time for the solver scales roughly with the number of compartments.
             Results for tumor growth are usually quite insensitive to the number of compartments
             used. The default setting is 50 compartments. The faster setting of 25 can be used
             when cell synchronisation effects are not important. The slower setting of 100 is
             primarily useful for exploring model sensitivity when synchronisation is significant.'),
           fluidRow(
             column(3,
                    h5('Compartments N')
             ),
             column(3,
                    selectInput('ncompstr', label = NULL, 
                                choices=c(25,50,100), 
                                selected=50)
             )
           ),
           fluidRow(
             column(6,
                    h5(textOutput('tdoubstdev'))
             )
           ),
           fluidRow(
             column(6,
                    h5(textOutput('discret'))
             )
           ),
           hr() 
           )
         )
),
tabPanel("Notes",
         a(img(src = "SystemsForecastingBanner.png", height = '10%', width = '100%'),
           target="_blank",href="http://www.systemsforecasting.com"),
         shinyUI(fluidPage(
           h4('Description'),
           p('The CellCycler is a parsimonius model of a growing tumour. The aim is to 
            capture the key dynamics of tumour behaviour so as to simulate and predict the 
            effect of treatment, either in vitro, in vivo, or in the clinic. The model uses 
            circa 100 ordinary differential equations to simulate cells as they pass through 
            the phases of the cell cycle. However the guiding philosophy of the model is that 
            it should only use parameters that can be observed or reasonably well approximated.
            The level of analysis is limited to cell state observables such as phase, apoptosis,
            and damage.'),
            h4('Cells'),
            p('The starting point in the Shiny program is the Cells page, which is used to 
            model the dynamics of a growing cell population. The key parameters are the average cell 
            doubling time in hours, and the fraction spent in each phase (G2 is set automatically since 
            the proportions must add to 1). The doubling time is assumed to be variable, with
            a range that depends on the number of model compartments. This can be adjusted in 
            the Advanced tab.'),
            p('In addition the user selects the simulation time in hours, 
            and plotting choices such as growing or damaged cells. The simulation is not 
            reactive, so the user must press the Run simulation button to run the model. 
            The plot will then show the volume of cells in each phase, as well as the total 
            volume, normalised to an initial volume of 1. Model settings can be saved to or 
            loaded from a csv file using the buttons at the bottom of the page.'),
            h4('PK'),
            p('The next pages, PK1 and PK2, are used to parameterise the PK models and drug effects. 
            The program has a choice of simple decay (K-PD), step, or one or two-compartment models, with 
            adjustable parameters such as Dose/Volume. The dosing schedule can be input as
            a comma separated sequence of times in hours. The step model has constant dose
            which turns on and off at the schedule times. In addition 
            the phase of action (choices are G1, S, G2, M, or all), and rates for death, 
            damage, and repair can be adjusted.'),
            p('The schematic below shows the different compartment models. The parameters are
            ke (elimination rate), ka (absorption by central compartment), and kp (transfer
            to peripheral compartment).'),
            img(src = "PKschematic.png", height = '70%', width = '70%'),
            h4('Tumor'),
            p('Finally, the Tumor page uses the model simulation to generate a plot of 
            tumor radius, given an initial radius and growing layer. A table is shown giving 
            total radius gain; the maximum gain that would be obtained in the absence of 
            drug; the radius loss due to drug; and the proportions of this loss that are 
            due to death or cell damage.'),
            p('The results can be compared with experimental results by using the Read Data 
            From File option. This can be in one of two forms. The first is a text file 
            containing a matrix, where the top row is times in hours, and the other rows are 
            individual experiments with the radius at each time. The second option is a text
            file with four columns: experiment ID, treatment time in days, treatment label 
            e.g. "untreated", and volume. Data from the adenoid cystic carcinoma
            study cited below is loaded by default.'),
            p('If the user selects show linear fit, then the plot will include a linear 
            interpolation, with estimates for initial volume and growing layer thickness given 
            in the table. If the data corresponds to control curves with no drug, then these 
            values can be used to parameterise the model.'),
            p('For cases with drug, a two-piece linear fit to the data should be done 
            separately. The first line segment corresponds to growth before drug, with a slope 
            corresponding to that of the control, and the second segment corresponds to when 
            growth has resumed after treatment, which again should have the same slope but a 
            lower y-intercept. The vertical distance between the two line segments can then 
            be used to parameterise the radius lost due to drug.'),
           hr(),
           a("Model equations and background (pdf)",target="_blank",href="CellCyclerDetails.pdf"),
           hr(),
           p('For information on the data file ACCX16_TRT.csv, see:'),
           a("Development and characterization of xenograft model systems for adenoid cystic carcinoma 
             (Moskaluk et al., 2011)",target="_blank",
             href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4151120/"),
           hr(),
           a("Contact Systems Forecasting",target="_blank",href="http://systemsforecasting.com/contact/"),
           hr() 
            )
         )
)


)
)