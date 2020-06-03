# Rshiny app for exploring metrics of tree shape across many different simulation models.

library(dplyr)
library(shiny)
library(shinyBS)
library(corrplot)
library(stringr)

source("code/treeMetrics.R")
source("code/pcaFunctions.r")

# Read in tree metrics output
treeOutput = read.table("treeOutput.txt", header = T, sep = '\t', stringsAsFactors = FALSE)
empiricalData = read.table("amphibianAndMammalCladeMetrics.txt", header = T, sep = '\t', stringsAsFactors = FALSE)

# Read in model parameter names and values for all models
paramfiles = list.files("parameters")
allModelParameterNames = list()
allModelParameterValues = list()
i = 0
for (pfile in paramfiles) {
  i = i + 1
  paramVals = read.csv(paste("parameters/", pfile, sep = ""), header = T, stringsAsFactors = FALSE)
  allModelParameterValues[[i]] = paramVals
  allModelParameterNames[[i]] = names(paramVals)[3:ncol(paramVals)]
  names(allModelParameterValues)[[i]] = str_extract(pfile, "^[A-Za-z]*")
  names(allModelParameterNames)[[i]] = str_extract(pfile, "^[A-Za-z]*")
}



# Function for creating descriptive text over radio buttons upon hovering
radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){
  
  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}




# Define UI for app  ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Exploring simulated trees"),
  
  # Create tabs
  # 1 Between simulation variable plots
  # 2 Between simulation PCA plots
  # 3 Within simulation variable plots
  
  tabsetPanel(
    
    # TAB 1: between simulation comparisons of individual tree metrics
    tabPanel("Between Sims: metric space",
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 helpText("Explore how various tree shape metrics vary by model and with model categories"),
                 
                 selectInput(inputId = "xvar",
                             label = "X-axis",
                             choices = c("Beta",
                                         "Gamma",
                                         "log10 Richness",
                                         "PD",
                                         "Colless",
                                         "Sackin",
                                         "Ratio of Yule / PDA likelihoods",
                                         "Mean Root Distance",
                                         "Variance in Root Distance",
                                         "Phylogenetic Species Variability",
                                         "mean I'",
                                         "Mean Pairwise Distance",
                                         "RPANDA spectral: principal eigenvalue",
                                         "RPANDA spectral: asymmetry",
                                         "RPANDA spectral: peakedness",
                                         "RPANDA spectral: eigengap",
                                         "nLTT_stat")),
                 
                 selectInput(inputId = "yvar",
                             label = "Y-axis",
                             choices = c("Gamma",
                                         "Beta", 
                                         "log10 Richness",
                                         "PD",
                                         "Colless",
                                         "Sackin",
                                         "Ratio of Yule / PDA likelihoods",
                                         "Mean Root Distance",
                                         "Variance in Root Distance",
                                         "Phylogenetic Species Variability",
                                         "mean I'",
                                         "Mean Pairwise Distance",
                                         "RPANDA spectral: principal eigenvalue",
                                         "RPANDA spectral: asymmetry",
                                         "RPANDA spectral: peakedness",
                                         "RPANDA spectral: eigengap",
                                         "nLTT_stat")),
                 
                 radioButtons(inputId = "colorBy",
                              label = "Trees color coded by:",
                              choices = c("Model", 
                                          "Model Family",
                                          "Entity Modeled",
                                          "Diversity Dependence",
                                          "Sympatric",
                                          "Allopatric",
                                          "Point Mutation Speciation")),
                 
                 radioTooltip(id = "colorBy", choice = "Model", title = "Model author or abbreviation", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Model Family", title = "Statistical, Conceptual, or Realistical", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Entity Modeled", title = "Individuals, Populations, or Species", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Diversity Dependence", title = "Is diversity dependence modeled (either implicitly or explicitly) or not?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Sympatric", title = "Can speciation occur sympatrically in this model?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Allopatric", title = "Can speciation occur allopatrically in this model?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Point Mutation Speciation", title = "Can speciation occur via point mutation in this model?", placement = "right", trigger = "hover"),
                 
                 radioButtons(inputId = "pchBy",
                              label = "Tree symbols coded by:",
                              choices = c("Model",
                                          "Model Family",
                                          "Entity Modeled",
                                          "Diversity Dependence",
                                          "Sympatric",
                                          "Allopatric",
                                          "Point Mutation Speciation")),
                 
                 radioTooltip(id = "pchBy", choice = "Model", title = "Model author or abbreviation", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Model Family", title = "Statistical, Conceptual, or Realistical", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Entity Modeled", title = "Individuals, Populations, or Species", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Diversity Dependence", title = "Is diversity dependence modeled (either implicitly or explicitly) or not?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Sympatric", title = "Can speciation occur sympatrically in this model?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Allopatric", title = "Can speciation occur allopatrically in this model?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Point Mutation Speciation", title = "Can speciation occur via point mutation in this model?", placement = "right", trigger = "hover"),
                 
                 
                 checkboxGroupInput(inputId = "modelsToInclude",
                                    label = "Models to include in the comparison",
                                    choices = list("Yule (yu)" = 'yu',
                                                   "PDA (pda)" = 'pda',
                                                   "TreeSim (oh)" = 'oh',
                                                   "DAISIE (ve)" = 've',
                                                   "Hurlbert-Stegen (hs)" = 'hs',
                                                   "Pontarp (pontarp)" = 'pontarp',
                                                   "Hartig (fh)" = 'fh',
                                                   "Coelho et al. (mt)" = 'mt',
                                                   "Leprieur et al. (split)" = 'split',
                                                   "Rangel (ra)" = 'ra',
                                                   "gen3sis (ga)" = 'ga',
                                                   "Cabral et al. (ca)" = 'ca',
                                                   "Xu-Etienne (XE)" = 'XE'),
                                    selected = list("Yule (yu)" = 'yu',
                                                    "PDA (pda)" = 'pda',
                                                    "TreeSim (oh)" = 'oh',
                                                    "DAISIE (ve)" = 've',
                                                    "Hurlbert-Stegen (hs)" = 'hs',
                                                    "Pontarp (pontarp)" = 'pontarp',
                                                    "Hartig (fh)" = 'fh',
                                                    "Coelho et al. (mt)" = 'mt',
                                                    "Leprieur et al. (split)" = 'split',
                                                    "Rangel (ra)" = 'ra',
                                                    "gen3sis (ga)" = 'ga',
                                                    "Cabral et al. (ca)" = 'ca',
                                                    "Xu-Etienne (XE)" = 'XE')),
                 
                 sliderInput("alphaSlider", "Transparency",
                             min = 0, max = 255, value = 200),
                 
                 radioButtons(inputId = "empiricalDataSwitch1",
                              label = "Overlay empirical data?",
                              choices = c("No", 
                                          "Yes")),
                 
                 sliderInput("minRichness1", "Minimum size empirical tree",
                             min = 20, max = max(500), value = 20)
                 
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 plotOutput("varPlot", height = 600)
                 
               )
             ) #end sidebarLayout
    ), #end tab 1
    
    # TAB 2
    tabPanel("Between Sims: PCA space",
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 helpText("First, choose which variables to include in the PCA"),
                 
                 checkboxGroupInput(inputId = "pcaVars",
                                    label = "Tree metrics",
                                    choices = list("Beta" = "Beta",
                                                   "Gamma" = "Gamma",
                                                   "log10 Richness" = "log10S",
                                                   "PD" = "PD",
                                                   "Colless" = "Colless",
                                                   "Sackin" = "Sackin",
                                                   "Ratio of Yule / PDA likelihoods" = "Yule.PDA.ratio",
                                                   "Mean Root Distance" = "MRD",
                                                   "Variance in Root Distance" = "VRD",
                                                   "Phylogenetic Species Variability" = "PSV",
                                                   "mean I'" = "mean.Iprime",
                                                   "Mean Pairwise Distance" = "MPD",
                                                   "RPANDA spectral: principal eigenvalue" = "MGL_principal_eigenvalue",
                                                   "RPANDA spectral: asymmetry" = "MGL_asymmetry",
                                                   "RPANDA spectral: peakedness" = "MGL_peakedness",
                                                   "RPANDA spectral: eigengap" = "MGL_eigengap",
                                                   "nLTT_stat" = "nLTT_stat"),
                                    selected = list("Beta" = "Beta",
                                                    "Gamma" = "Gamma",
                                                    "Mean Root Distance" = "MRD",
                                                    "Variance in Root Distance" = "VRD",
                                                    "Mean Pairwise Distance" = "MPD",
                                                    "RPANDA spectral: peakedness" = "MGL_peakedness")),
                 
                 selectInput(inputId = "xvar2",
                             label = "X-axis",
                             choices = c("PC1",
                                         "PC2",
                                         "PC3",
                                         "PC4")),
                 
                 selectInput(inputId = "yvar2",
                             label = "Y-axis",
                             choices = c("PC2",
                                         "PC1",
                                         "PC3",
                                         "PC4")),
                 
                 radioButtons(inputId = "PCAcolorBy",
                              label = "Trees color coded by:",
                              choices = c("Model", 
                                          "Model Family",
                                          "Entity Modeled",
                                          "Diversity Dependence",
                                          "Sympatric",
                                          "Allopatric",
                                          "Point Mutation Speciation")),
                 
                 radioTooltip(id = "colorBy", choice = "Model", title = "Model author or abbreviation", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Model Family", title = "Statistical, Conceptual, or Realistical", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Entity Modeled", title = "Individuals, Populations, or Species", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Diversity Dependence", title = "Is diversity dependence modeled (either implicitly or explicitly) or not?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Sympatric", title = "Can speciation occur sympatrically in this model?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Allopatric", title = "Can speciation occur allopatrically in this model?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "colorBy", choice = "Point Mutation Speciation", title = "Can speciation occur via point mutation in this model?", placement = "right", trigger = "hover"),
                 
                 radioButtons(inputId = "PCApchBy",
                              label = "Tree symbols coded by:",
                              choices = c("Model",
                                          "Model Family",
                                          "Entity Modeled",
                                          "Diversity Dependence",
                                          "Sympatric",
                                          "Allopatric",
                                          "Point Mutation Speciation")),
                 
                 radioTooltip(id = "pchBy", choice = "Model", title = "Model author or abbreviation", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Model Family", title = "Statistical, Conceptual, or Realistical", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Entity Modeled", title = "Individuals, Populations, or Species", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Diversity Dependence", title = "Is diversity dependence modeled (either implicitly or explicitly) or not?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Sympatric", title = "Can speciation occur sympatrically in this model?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Allopatric", title = "Can speciation occur allopatrically in this model?", placement = "right", trigger = "hover"),
                 radioTooltip(id = "pchBy", choice = "Point Mutation Speciation", title = "Can speciation occur via point mutation in this model?", placement = "right", trigger = "hover"),
                 
                 
                 checkboxGroupInput(inputId = "modelsToIncludePCA",
                                    label = "Models to include in the comparison",
                                    choices = list("Yule (yu)" = 'yu',
                                                   "PDA (pda)" = 'pda',
                                                   "TreeSim (oh)" = 'oh',
                                                   "DAISIE (ve)" = 've',
                                                   "Hurlbert-Stegen (hs)" = 'hs',
                                                   "Pontarp (pontarp)" = 'pontarp',
                                                   "Hartig (fh)" = 'fh',
                                                   "Coelho et al. (mt)" = 'mt',
                                                   "Leprieur et al. (split)" = 'split',
                                                   "Rangel (ra)" = 'ra',
                                                   "gen3sis (ga)" = 'ga',
                                                   "Cabral et al. (ca)" = 'ca',
                                                   "Xu-Etienne (XE)" = 'XE'),
                                    selected = c('yu', 'pda', 'oh', 've', 'hs', 'pontarp', 'fh', 'mt', 'split', 'ra', 'ga', 'ca', 'XE')),
                 
                 sliderInput("alphaSlider2", "Transparency",
                             min = 0, max = 255, value = 200),
                 
                 radioButtons(inputId = "empiricalDataSwitch2",
                              label = "Overlay empirical data?",
                              choices = c("No", 
                                          "Yes")),
                 
                 sliderInput("minRichness2", "Minimum size empirical tree",
                             min = 20, max = 500, value = 20)
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 fluidRow( 
                   verticalLayout( 
                     
                     plotOutput("pcaPlot", height = 600),
                     
                     helpText("Refer to the correlation matrix below for how strongly variables covary."),
                     
                     plotOutput("corPlot", height = 600)
                     )
                 ) 
               )
               
             ) #end sidebarLayout
    ), #end tab 2
    
    # TAB 3
    tabPanel("Individual Model Exploration",
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 helpText("Choose an individual model to explore"),
                 
                 selectInput(inputId = "model",
                             label = "Select a model",
                             choices = c("DAISIE (ve)",
                                         "TreeSim (oh)",
                                         "Hurlbert-Stegen (hs)",
                                         "Pontarp (pontarp)",
                                         "Hartig (fh)",
                                         "Coelho et al. (mt)",
                                         "Leprieur et al. (split)",
                                         "Rangel (ra)",
                                         "gen3sis (ga)",
                                         "Cabral et al. (ca)",
                                         "Xu-Etienne (XE)")),
                 
                 helpText("Explore how various tree shape metrics vary with model parameters"),
                 
                 selectInput(inputId = "xvar3",
                             label = "X-axis",
                             choices = c("Beta",
                                         "Gamma",
                                         "log10 Richness",
                                         "PD",
                                         "Colless",
                                         "Sackin",
                                         "Ratio of Yule / PDA likelihoods",
                                         "Mean Root Distance",
                                         "Variance in Root Distance",
                                         "Phylogenetic Species Variability",
                                         "mean I'",
                                         "Mean Pairwise Distance",
                                         "RPANDA spectral: principal eigenvalue",
                                         "RPANDA spectral: asymmetry",
                                         "RPANDA spectral: peakedness",
                                         "RPANDA spectral: eigengap",
                                         "nLTT_stat")),
                 
                 selectInput(inputId = "yvar3",
                             label = "Y-axis",
                             choices = c("Gamma",
                                         "Beta", 
                                         "log10 Richness",
                                         "PD",
                                         "Colless",
                                         "Sackin",
                                         "Ratio of Yule / PDA likelihoods",
                                         "Mean Root Distance",
                                         "Variance in Root Distance",
                                         "Phylogenetic Species Variability",
                                         "mean I'",
                                         "Mean Pairwise Distance",
                                         "RPANDA spectral: principal eigenvalue",
                                         "RPANDA spectral: asymmetry",
                                         "RPANDA spectral: peakedness",
                                         "RPANDA spectral: eigengap",
                                         "nLTT_stat")),
                 
                 uiOutput("parColorBy"),
                 
                 uiOutput("parPchBy"),
                 
                 sliderInput("alphaSlider3", "Transparency",
                             min = 0, max = 255, value = 200),
                 
                 radioButtons(inputId = "empiricalDataSwitch3",
                              label = "Overlay empirical data?",
                              choices = c("No", 
                                          "Yes")),
                 
                 sliderInput("minRichness3", "Minimum size empirical tree",
                             min = 20, max = 500, value = 20)
                 
                 
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 fluidRow( 
                   verticalLayout( 
                     
                     plotOutput("varParPlot", height = 600)
                     
                   )
                 ) 
               )
               
             ) #end sidebarLayout
    ), #end tab 3
    
    # TAB 4
    tabPanel("Metric boxplots",
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 helpText("Choose the metric to compare. Remember that some of these metrics are richness-dependent."),
                 
                 radioButtons(inputId = "metric",
                             label = "Metric",
                             choices = c("Beta",
                                         "Gamma",
                                         "log10 Richness",
                                         "PD",
                                         "Colless",
                                         "Sackin",
                                         "Ratio of Yule / PDA likelihoods",
                                         "Mean Root Distance",
                                         "Variance in Root Distance",
                                         "Phylogenetic Species Variability",
                                         "mean I'",
                                         "Mean Pairwise Distance",
                                         "RPANDA spectral: principal eigenvalue",
                                         "RPANDA spectral: asymmetry",
                                         "RPANDA spectral: peakedness",
                                         "RPANDA spectral: eigengap",
                                         "nLTT_stat")),
                 
                 radioButtons(inputId = "log",
                              label = "Log10 transform?",
                              choices = c("No",
                                          "Yes"))
                                          
                 
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 fluidRow( 
                   verticalLayout( 
                     
                     plotOutput("boxPlot", height = 500)
                     
                   )
                 ) 
               )
               
             ) #end sidebarLayout
    ) #end tab 4
    
    
  ) #end tabSetPanel
) #end UI





# Define server logic ----
server <- function(input, output, session) {
  
  # Plot trees based on two selected shape metrics
  output$varPlot <- renderPlot({
    
    xvar <- switch(input$xvar, 
                   "log10 Richness" = 'log10S',
                   "PD" = 'PD',
                   "Gamma" = 'Gamma',
                   "Beta" = 'Beta', 
                   "Colless" = 'Colless',
                   "Sackin" = 'Sackin',
                   "Ratio of Yule / PDA likelihoods" = 'Yule.PDA.ratio',
                   "Mean Root Distance" = 'MRD',
                   "Variance in Root Distance" = 'VRD',
                   "Phylogenetic Species Variability" = 'PSV',
                   "mean I'" = 'mean.Iprime',
                   "Mean Pairwise Distance" = 'MPD',
                   "RPANDA spectral: principal eigenvalue" = 'MGL_principal_eigenvalue',
                   "RPANDA spectral: asymmetry" = 'MGL_asymmetry',
                   "RPANDA spectral: peakedness" = 'MGL_peakedness',
                   "RPANDA spectral: eigengap" = 'MGL_eigengap',
                   "nLTT_stat" = 'nLTT_stat')
    
    yvar <- switch(input$yvar, 
                   "log10 Richness" = 'log10S',
                   "PD" = 'PD',
                   "Gamma" = 'Gamma',
                   "Beta" = 'Beta', 
                   "Colless" = 'Colless',
                   "Sackin" = 'Sackin',
                   "Ratio of Yule / PDA likelihoods" = 'Yule.PDA.ratio',
                   "Mean Root Distance" = 'MRD',
                   "Variance in Root Distance" = 'VRD',
                   "Phylogenetic Species Variability" = 'PSV',
                   "mean I'" = 'mean.Iprime',
                   "Mean Pairwise Distance" = 'MPD',
                   "RPANDA spectral: principal eigenvalue" = 'MGL_principal_eigenvalue',
                   "RPANDA spectral: asymmetry" = 'MGL_asymmetry',
                   "RPANDA spectral: peakedness" = 'MGL_peakedness',
                   "RPANDA spectral: eigengap" = 'MGL_eigengap',
                   "nLTT_stat" = 'nLTT_stat')
    
    colorBy <- switch(input$colorBy, 
                      "Model" = 'model', 
                      "Model Family" = 'ModelFamily',
                      "Entity Modeled" = 'EntityModeled',
                      "Diversity Dependence" = 'DiversityDependence',
                      "Sympatric" = 'Sympatric',
                      "Allopatric" = 'Allopatric',
                      "Point Mutation Speciation" = 'PointMutation')
    
    pchBy <- switch(input$pchBy, 
                    "Model" = 'model', 
                    "Model Family" = 'ModelFamily',
                    "Entity Modeled" = 'EntityModeled',
                    "Diversity Dependence" = 'DiversityDependence',
                    "Sympatric" = 'Sympatric',
                    "Allopatric" = 'Allopatric',
                    "Point Mutation Speciation" = 'PointMutation')
    
    
    plottingOutput = filter(treeOutput, model %in% input$modelsToInclude) %>%
      dplyr::select(model, simID, !!xvar, !!yvar) %>%
      na.omit()
    
    par(mar = c(5, 5, 1, 1), cex.lab = 1.75)
    betweenModelVarPlot(treeOutput = plottingOutput, 
                        xvar = xvar, 
                        yvar = yvar, 
                        colorBy = colorBy, 
                        pchBy = pchBy, 
                        alpha = input$alphaSlider,
                        cex = 2)  
    
    if (input$empiricalDataSwitch1 == "Yes") {
      empData = filter(empiricalData, S >= input$minRichness1)
      points(empData[, xvar], empData[, yvar], pch = 15, col = rgb(.3, .3, .3, .3, maxColorValue = 1))
    }
  })
  

  # Plot trees in PCA space
  output$pcaPlot <- renderPlot({
    
    xvar2 <- switch(input$xvar2, 
                   "PC1" = 1,
                   "PC2" = 2,
                   "PC3" = 3,
                   "PC4" = 4)
    
    yvar2 <- switch(input$yvar2, 
                    "PC1" = 1,
                    "PC2" = 2,
                    "PC3" = 3,
                    "PC4" = 4)
    
    PCAcolorBy <- switch(input$PCAcolorBy, 
                  "Model" = 'model', 
                  "Model Family" = 'ModelFamily',
                  "Entity Modeled" = 'EntityModeled',
                  "Diversity Dependence" = 'DiversityDependence',
                  "Sympatric" = 'Sympatric',
                  "Allopatric" = 'Allopatric',
                  "Point Mutation Speciation" = 'PointMutation')

    PCApchBy <- switch(input$PCApchBy, 
                    "Model" = 'model', 
                    "Model Family" = 'ModelFamily',
                    "Entity Modeled" = 'EntityModeled',
                    "Diversity Dependence" = 'DiversityDependence',
                    "Sympatric" = 'Sympatric',
                    "Allopatric" = 'Allopatric',
                    "Point Mutation Speciation" = 'PointMutation')
    
    
    plottingOutput = filter(treeOutput, model %in% input$modelsToIncludePCA)
    
    pcaOutput = treeMetricsPCA(plottingOutput, models = input$modelsToInclude, vars = input$pcaVars)

    par(mar = c(5, 5, 1, 1), cex.lab = 1.75)
    betweenModelPCAPlot(pcaOutput, 
                        xscore = xvar2, 
                        yscore = yvar2, 
                        colorBy = PCAcolorBy, 
                        pchBy = PCApchBy, 
                        alpha = input$alphaSlider2,
                        cex = 2)  
    
    if (input$empiricalDataSwitch2 == "Yes") {
      loadings = pcaOutput$pcaLoadings
      empDataMatrix = empiricalData[, input$pcaVars] %>%
        data.matrix() %>%
        na.omit()
      
      
      empScores = empDataMatrix %*% loadings
      
      points(empScores[, xvar2], empScores[, yvar2], pch = 15, col = rgb(.3, .3, .3, .3, maxColorValue = 1))  
      
    }
    
  })
  
  
  # Plot a correlation matrix of selected tree metric variables
  output$corPlot <- renderPlot({
    
    corData = treeOutput[, input$pcaVars]
    
    varCor = cor(corData, use = "pairwise.complete.obs")
    corrplot(varCor)
    
  })
  
  
  # Parameter values for coding depend on the model selected
  output$parColorBy <- renderUI({
    model <- switch(input$model, 
                    "TreeSim (oh)" = 'oh',
                    "DAISIE (ve)" = 've',
                    "Hurlbert-Stegen (hs)" = 'hs',
                    "Pontarp (pontarp)" = 'pontarp',
                    "Hartig (fh)" = 'fh',
                    "Coelho et al. (mt)" = 'mt',
                    "Leprieur et al. (split)" = 'split',
                    "Rangel (ra)" = 'ra',
                    "gen3sis (ga)" = 'ga',
                    "Cabral et al. (ca)" = 'ca',
                    "Xu-Etienne (XE)" = 'XE')
    
    selectInput(inputId = "parColor",
                label = "Model parameter to code by color",
                choices = allModelParameterNames[[model]],
                selected = allModelParameterNames[[model]][1])
  })

  
  output$parPchBy <- renderUI({
    model <- switch(input$model, 
                    "TreeSim (oh)" = 'oh',
                    "DAISIE (ve)" = 've',
                    "Hurlbert-Stegen (hs)" = 'hs',
                    "Pontarp (pontarp)" = 'pontarp',
                    "Hartig (fh)" = 'fh',
                    "Coelho et al. (mt)" = 'mt',
                    "Leprieur et al. (split)" = 'split',
                    "Rangel (ra)" = 'ra',
                    "gen3sis (ga)" = 'ga',
                    "Cabral et al. (ca)" = 'ca',
                    "Xu-Etienne (XE)" = 'XE')
    
    selectInput(inputId = "parPch",
                label = "Model parameter to code by symbol",
                choices = allModelParameterNames[[model]],
                selected = allModelParameterNames[[model]][2])
  })
  
  
  # Plot trees based on two selected shape metrics WITHIN AN INDIVIDUAL MODEL
  output$varParPlot <- renderPlot({
    
    model <- switch(input$model, 
                    "TreeSim (oh)" = 'oh',
                    "DAISIE (ve)" = 've',
                    "Hurlbert-Stegen (hs)" = 'hs',
                    "Pontarp (pontarp)" = 'pontarp',
                    "Hartig (fh)" = 'fh',
                    "Coelho et al. (mt)" = 'mt',
                    "Leprieur et al. (split)" = 'split',
                    "Rangel (ra)" = 'ra',
                    "gen3sis (ga)" = 'ga',
                    "Cabral et al. (ca)" = 'ca',
                    "Xu-Etienne (XE)" = 'XE')
    
    xvar3 <- switch(input$xvar3, 
                   "log10 Richness" = 'log10S',
                   "PD" = 'PD',
                   "Gamma" = 'Gamma',
                   "Beta" = 'Beta', 
                   "Colless" = 'Colless',
                   "Sackin" = 'Sackin',
                   "Ratio of Yule / PDA likelihoods" = 'Yule.PDA.ratio',
                   "Mean Root Distance" = 'MRD',
                   "Variance in Root Distance" = 'VRD',
                   "Phylogenetic Species Variability" = 'PSV',
                   "mean I'" = 'mean.Iprime',
                   "Mean Pairwise Distance" = 'MPD',
                   "RPANDA spectral: principal eigenvalue" = 'MGL_principal_eigenvalue',
                   "RPANDA spectral: asymmetry" = 'MGL_asymmetry',
                   "RPANDA spectral: peakedness" = 'MGL_peakedness',
                   "RPANDA spectral: eigengap" = 'MGL_eigengap',
                   "nLTT_stat" = 'nLTT_stat')
    
    yvar3 <- switch(input$yvar3, 
                   "log10 Richness" = 'log10S',
                   "PD" = 'PD',
                   "Gamma" = 'Gamma',
                   "Beta" = 'Beta', 
                   "Colless" = 'Colless',
                   "Sackin" = 'Sackin',
                   "Ratio of Yule / PDA likelihoods" = 'Yule.PDA.ratio',
                   "Mean Root Distance" = 'MRD',
                   "Variance in Root Distance" = 'VRD',
                   "Phylogenetic Species Variability" = 'PSV',
                   "mean I'" = 'mean.Iprime',
                   "Mean Pairwise Distance" = 'MPD',
                   "RPANDA spectral: principal eigenvalue" = 'MGL_principal_eigenvalue',
                   "RPANDA spectral: asymmetry" = 'MGL_asymmetry',
                   "RPANDA spectral: peakedness" = 'MGL_peakedness',
                   "RPANDA spectral: eigengap" = 'MGL_eigengap',
                   "nLTT_stat" = 'nLTT_stat')
    
    
    par(mar = c(5, 5, 1, 1), cex.lab = 1.75)
    withinModelVarPlot(treeOutput = treeOutput,
                       modelAbbrev = model,
                       modelParams = allModelParameterValues[[model]],
                        xvar = xvar3, 
                        yvar = yvar3, 
                        colorBy = input$parColor, 
                        pchBy = input$parPch, 
                        alpha = input$alphaSlider3,
                        cex = 2) 
    
    if (input$empiricalDataSwitch3 == "Yes") {
      empData = filter(empiricalData, S >= input$minRichness1)
      points(empData[, xvar3], empData[, yvar3], pch = 15, col = rgb(.3, .3, .3, .3, maxColorValue = 1))
    }
  })
  
  output$boxPlot <- renderPlot({
    metric <- switch(input$metric, 
                   "log10 Richness" = 'log10S',
                   "PD" = 'PD',
                   "Gamma" = 'Gamma',
                   "Beta" = 'Beta', 
                   "Colless" = 'Colless',
                   "Sackin" = 'Sackin',
                   "Ratio of Yule / PDA likelihoods" = 'Yule.PDA.ratio',
                   "Mean Root Distance" = 'MRD',
                   "Variance in Root Distance" = 'VRD',
                   "Phylogenetic Species Variability" = 'PSV',
                   "mean I'" = 'mean.Iprime',
                   "Mean Pairwise Distance" = 'MPD',
                   "RPANDA spectral: principal eigenvalue" = 'MGL_principal_eigenvalue',
                   "RPANDA spectral: asymmetry" = 'MGL_asymmetry',
                   "RPANDA spectral: peakedness" = 'MGL_peakedness',
                   "RPANDA spectral: eigengap" = 'MGL_eigengap',
                   "nLTT_stat" = 'nLTT_stat')
    
    
    empirData = empiricalData
    empirData$model = 'empirical'
    allData = rbind(treeOutput, empirData)
    colors = c('turquoise', 'red', 'yellow2', 'darkblue', 'limegreen', 'hotpink', 'blue', 
               'purple', 'brown', 'seagreen', 'darkorange', 'pink', 'firebrick4', 'olivedrab1')
    nModels = length(unique(treeOutput$model))
    
    if (input$log == "Yes") {
      boxplot(log10(allData[, metric]) ~ allData$model, ylab = paste("log10", metric), 
              col = c('gray70', colors[1:nModels]), las = 2)
    } else {
      boxplot(allData[, metric] ~ allData$model, ylab = metric, 
              col = c('gray70', colors[1:nModels]), las = 2)
    }
  })
  
}


shinyApp(ui = ui, server = server)