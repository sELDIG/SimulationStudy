# Rshiny app for exploring metrics of tree shape across many different simulation models.

library(dplyr)
library(shiny)
library(shinyBS)

source("code/treeMetrics.R")
source("code/pcaFunctions.r")

treeOutput = read.table("treeOutput.txt", header = T, sep = '\t', stringsAsFactors = FALSE)

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
                                         "Richness",
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
                                         "Richness",
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
                                                   "DAISIE (etienne)" = 'etienne',
                                                   "Hurlbert-Stegen (hs)" = 'hs',
                                                   "Pontarp (pontarp)" = 'pontarp',
                                                   "Hartig (fh)" = 'fh',
                                                   "Coelho et al. (mt)" = 'mt',
                                                   "Leprieur et al. (split)" = 'split',
                                                   "Rangel (ra)" = 'ra',
                                                   "GaSM (ga)" = 'ga'),
                                    selected = list("Yule (yu)" = 'yu',
                                                    "PDA (pda)" = 'pda',
                                                    "TreeSim (oh)" = 'oh',
                                                    "DAISIE (etienne)" = 'etienne',
                                                    "Hurlbert-Stegen (hs)" = 'hs',
                                                    "Pontarp (pontarp)" = 'pontarp',
                                                    "Hartig (fh)" = 'fh',
                                                    "Coelho et al. (mt)" = 'mt',
                                                    "Leprieur et al. (split)" = 'split',
                                                    "Rangel (ra)" = 'ra',
                                                    "GaSM (ga)" = 'ga')),
                 
                 sliderInput("alphaSlider", h3("Transparency"),
                             min = 0, max = 255, value = 200)
                 
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 plotOutput("varPlot", height = 600)
                 
               )
             ) #end sidebarLayout
    ), #end tab
    
    # TAB 2
    tabPanel("Between Sims: PCA space",
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 helpText("Explore how principal components of tree shape vary with model categories"),
                 
                 selectInput(inputId = "xvar",
                             label = "X-axis",
                             choices = c("Beta",
                                         "Gamma",
                                         "Richness",
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
                                         "Richness",
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
                 
                 
                 checkboxGroupInput(inputId = "modelsToIncludePCA",
                                    label = "Models to include in the comparison",
                                    choices = list("Yule (yu)" = 'yu',
                                                   "PDA (pda)" = 'pda',
                                                   "TreeSim (oh)" = 'oh',
                                                   "DAISIE (etienne)" = 'etienne',
                                                   "Hurlbert-Stegen (hs)" = 'hs',
                                                   "Pontarp (pontarp)" = 'pontarp',
                                                   "Hartig (fh)" = 'fh',
                                                   "Coelho et al. (mt)" = 'mt',
                                                   "Leprieur et al. (split)" = 'split',
                                                   "Rangel (ra)" = 'ra',
                                                   "GaSM (ga)" = 'ga'),
                                    selected = c('yu', 'pda', 'oh', 'etienne', 'hs', 'pontarp', 'fh', 'mt', 'split', 'ra', 'ga')),
                 
                 sliderInput("alphaSlider", h3("Transparency"),
                             min = 0, max = 255, value = 200)
                 
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 fluidRow( 
                   verticalLayout( 
                     
                     plotOutput("pcaPlot"))
                 ) 
               )
               
             ) #end sidebarLayout
    ) #end tab
  ) #end tabSetPanel
) #end fluidPage



# Define server logic ----
server <- function(input, output, session) {
  
  #observeEvent(input$modelsToInclude, {
  #  modelsReactive = input$modelsToInclude
  #})
  
  output$varPlot <- renderPlot({
    
    xvar <- switch(input$xvar, 
                   "Richness" = 'S',
                   "PD" = 'PD',
                   "Gamma" = 'gamma',
                   "Beta" = 'beta', 
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
                   "Richness" = 'S',
                   "PD" = 'PD',
                   "Gamma" = 'gamma',
                   "Beta" = 'beta', 
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
    
    
    plottingOutput = filter(treeOutput, model %in% input$modelsToInclude)

    par(mar = c(4, 4, 1, 1), cex.lab = 1.75)
    betweenModelVarPlot(treeOutput = plottingOutput, 
                        xvar = xvar, 
                        yvar = yvar, 
                        colorBy = colorBy, 
                        pchBy = pchBy, 
                        alpha = input$alphaSlider,
                        cex = 2)  
  })
  

  #output$model_text = renderText({
  # paste("The set of models to plot are", paste(input$modelsToInclude, collapse = ", ")) 
  #})
  
}


shinyApp(ui = ui, server = server)