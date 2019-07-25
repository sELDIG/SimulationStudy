# Script for calculating metrics on all simulated phylogenies and conducting PCA
# Function that conducts PCA on treeOutput set of specified tree metrics
# -- vars is a vector of column names to include, default is all vars
# -- models is a vector of model abbreviations
treeMetricsPCA = function(treeOutput, models = 'all', vars = 'all') {

  if (models != 'all') {
    treeOutput = filter(treeOutput, model %in% models)
  }
  
  if (vars == 'all') {
    vars = names(treeOutput[, 3:ncol(treeOutput)])
  }
  outputSubsetNoNAs = na.omit(treeOutput[!names(treeOutput) == "VPD"])  
  
  pc = princomp(outputSubsetNoNAs[, names(outputSubsetNoNAs) %in% vars], cor = TRUE)
  pcaScores = cbind(outputSubsetNoNAs[, c("model", "simID")], pc$scores) 
  
  return(list(pcaScores = pcaScores, pcaLoadings = pc$loadings))
}




betweenModelPCAPlot = function(pcaOutput,          # dataframe with model, simID, and PC scores
                               xscore = 1,        # PC score plotted on the x-axis
                               yscore = 2,        # PC score plotted on the y-axis
                               colorBy = 'model', # any variable/parameter name by which to color points
                               pchBy = 'model'         # categorical variable/parameter name by which to specify point shape 
                               ) {
  
  # Reads in model classification codes and assigns colors and pchs
  require(gsheet)
  url = "https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=2047946073"
  
  modelClassification = gsheet2tbl(url)
  
  pcaScores = pcaOutput$pcaScores
  pcaLoadings = pcaOutput$pcaLoadings
  
  colors = c('turquoise', 'orangered', 'yellow2', 'darkblue', 'limegreen', 'magenta', 'blue', 'purple', 'brown', 'seagreen')
  pch = c(15, 16, 17, 18, 1, 7, 8, 10)

  colorCode = data.frame(val = unique(modelClassification[, colorBy]), 
                         color = colors[1:nrow(unique(modelClassification[, colorBy]))])
  colorCode$color = as.character(colorCode$color)
  names(colorCode)[1] = colorBy
  
  pchCode = data.frame(val = unique(modelClassification[, pchBy]),
                       pch = pch[1:nrow(unique(modelClassification[, pchBy]))])
  names(pchCode)[1] = pchBy
  
  plotOutput = left_join(pcaScores, modelClassification, by = c("model", "simID")) %>%
    left_join(colorCode, by = unname(colorBy)) %>%
    left_join(pchCode, by = unname(pchBy))

  plot(plotOutput[,paste("Comp.", xscore, sep = "")], plotOutput[, paste("Comp.", yscore, sep = "")], 
       col = plotOutput$color, cex = 2, 
       pch = plotOutput$pch, xlab = paste("PC", xscore), ylab = paste("PC", yscore), 
       ylim = c(-max(abs(range(plotOutput[,paste("Comp.", yscore, sep = "")]))), 
                max(abs(range(plotOutput[,paste("Comp.", yscore, sep = "")])))),
       xlim = c(-max(abs(range(plotOutput[,paste("Comp.", xscore, sep = "")]))), 
                max(abs(range(plotOutput[,paste("Comp.", xscore, sep = "")])))))
  legend("topleft", 
         legend = c(toupper(colorBy), colorCode[,1]), bty = "n",
         col = c('white', colorCode[,2]), pch = 16)
  
  legend("topright", 
         legend = c(toupper(pchBy), pchCode[,1]), bty = "n",
         col = c('white', rep('black', nrow(pchCode))), pch = c(16, pchCode[,2]))
  
  par(new = TRUE)
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       xlim = c(-1, 1), ylim = c(-1, 1))
  
  # Only print the top variable loadings since it otherwise becomes messy
  mainLoadings = pcaLoadings[order(pcaLoadings[, xscore], decreasing = TRUE) <= 4 | 
                               order(pcaLoadings[,yscore], decreasing = TRUE) <= 4, 
                             c(xscore, yscore)]
  text(mainLoadings[, 1], mainLoadings[, 2], row.names(mainLoadings))  
}



withinModelPCAPlot = function(pcaOutput,          # dataframe with model, simID, and PC scores
                              modelAbbrev,             # specify the abbreviation of the model to explore 
                              xscore = 1,        # PC score plotted on the x-axis
                               yscore = 2,        # PC score plotted on the y-axis
                               colorBy,        # any variable/parameter name by which to color points
                               pchBy            # categorical variable/parameter name by which to specify point shape 
) {

  pcaScores = pcaOutput$pcaScores %>%
    filter(model == modelAbbrev)
  
  modelParams = read.csv(paste('parameters/', modelAbbrev, '_parameters.csv', sep = ''), header = TRUE)
  
  modelScores = left_join(pcaScores, modelParams, by = "simID")
  
  if (is.numeric(modelScores[, colorBy])) {
    
    shades <- rainbow(130)[100:1]
    percents <- as.integer(cut(modelScores[, colorBy], 100, include.lowest = TRUE, ordered = TRUE))
    modelScores$color = shades[percents]
    modelOutput = modelScores
    
  } else {
    
    colors = c('turquoise', 'orangered', 'yellow2', 'darkblue', 'limegreen', 'magenta', 'blue', 'purple', 'brown', 'seagreen')
    colorCode = data.frame(val = unique(modelScores[, colorBy]), 
                           color = colors[1:length(unique(modelScores[, colorBy]))])
    colorCode$color = as.character(colorCode$color)
    names(colorCode)[1] = colorBy
    
    modelOutput = left_join(modelScores, colorCode, by = unname(colorBy))
  }
  
  pch = c(15, 16, 17, 18, 1, 7, 8, 10)
  if (class(pchBy) == 'numeric') {
    warning("Symbol size is best used for visualizing categorical variables")
  } else {
    pchCode = data.frame(val = unique(modelScores[, pchBy]),
                         pch = pch[1:length(unique(modelScores[, pchBy]))])
    names(pchCode)[1] = pchBy
    
  }
  
  plotOutput = left_join(modelOutput, pchCode, by = unname(pchBy))
  
  plot(plotOutput[,paste("Comp.", xscore, sep = "")], plotOutput[, paste("Comp.", yscore, sep = "")], 
       col = plotOutput$color, cex = 2, 
       pch = plotOutput$pch, xlab = paste("PC", xscore), ylab = paste("PC", yscore), 
       ylim = c(-max(abs(range(plotOutput[,paste("Comp.", yscore, sep = "")]))), 
                max(abs(range(plotOutput[,paste("Comp.", yscore, sep = "")])))),
       xlim = c(-max(abs(range(plotOutput[,paste("Comp.", xscore, sep = "")]))), 
                max(abs(range(plotOutput[,paste("Comp.", xscore, sep = "")])))))
  
  # color legend
  if (is.numeric(modelScores[, colorBy])) {
    maxvar = max(modelScores[, colorBy], na.rm = TRUE)
    minvar = min(modelScores[, colorBy], na.rm = TRUE)
    inc <- (maxvar - minvar) / 4
    legend.text <- round(c(minvar, minvar + inc, minvar + 2 * inc, minvar + 3 * inc, maxvar))
    
    legend("topleft", 
           legend = c(toupper(colorBy), legend.text), bty = "n",
           fill = c('white', shades[c(1, 25, 50, 75, 100)]))

  } else {

    legend("topleft", legend = c(toupper(colorBy), as.character(colorCode[,1])), bty = "n",
           col = c('white', as.character(colorCode[,2])), pch = 16)
  }

  legend("topright", 
         legend = c(toupper(pchBy), as.character(pchCode[,1])), bty = "n",
         col = c('white', rep('black', nrow(pchCode))), pch = c(16, pchCode[,2]))
  
  par(new = TRUE)
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       xlim = c(-1, 1), ylim = c(-1, 1))
  
  # Only print the top variable loadings since it otherwise becomes messy
  mainLoadings = pcaLoadings[order(pcaLoadings[, xscore], decreasing = TRUE) <= 4 | 
                               order(pcaLoadings[,yscore], decreasing = TRUE) <= 4, 
                             c(xscore, yscore)]
  text(mainLoadings[, 1], mainLoadings[, 2], row.names(mainLoadings))  
}




