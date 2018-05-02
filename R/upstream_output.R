#' @export
scorePlot = function(aResult, plotorder = "score"){
  aFull = ldply(aResult, .fun = function(.)return(data.frame(.$aResult, mxRank = .$mxRank) ))
  sp = makeScorePlot(aFull, plotorder)
}

#' @export
makeScorePlot = function (aFrame, plotorder = "score"){

  aFrame = aFrame %>% group_by(ClassName) %>% mutate(meanScore = mean(combinedScore, na.rm = TRUE), meanStat = mean(NormalizedSetStat, na.rm = TRUE))
  aFrame = aFrame %>% group_by(mxRank) %>% mutate(ClassRank = rank(-combinedScore))

  if (plotorder == "score"){
    orderx = aFrame$combinedScore
  }

  if (plotorder == "stat"){
    orderx = -aFrame$NormalizedSetStat
  }

  prt2 = ggplot(aFrame, aes(x= reorder(ClassName, orderx, median), y = NormalizedSetStat)) + geom_boxplot()
  prt2 = prt2 + geom_point(aes(colour = rescale(pFeatureScore, from = c(0,2), clip = TRUE),
                               size   = rescale(nFeatures, from = c(1,30), to = c(1,30), clip = TRUE)))
  prt2 = prt2 + geom_abline(slope = 0)
  prt2 = prt2 + coord_flip() + theme_bw()
  prt2 = prt2 + scale_colour_gradientn(name = "Specificity Score",space = "rgb",
                                       colours = c("black", "white", "red"),
                                       values = c(0,1.3, 2)/2,
                                       breaks = c(0,1.3, 2)/2,
                                       labels = c(0, 1.3, 2),
                                       limits = c(0,1))
  prt2 = prt2 + scale_size(breaks = c(10,20,31), labels = c("10", "20", ">30"))
  prt2 = prt2 + ylab("Normalized kinase statistic") + xlab("Kinase name")
  yRange = layer_scales(prt2)$y$range$range
  prt2 = prt2 + scale_y_continuous(position = "top", limits = c( min(0, yRange[1]), max(0, yRange[2])) )
  prt2 = prt2 + guides( size = guide_legend("Peptide set size"), color = guide_colorbar("Specificity score"))
  return(prt2)
}

#' @export
makeVolcanoPlot = function(aSummary){
  vp = ggplot(aSummary, aes(
                              x = meanStat,
                              y = medianScore,
                              colour = rescale(meanFeatScore, from = c(0,2), clip = TRUE),
                              label = ClassName,
                              size = meanSetSize,
                              alpha =  1-rescale(sdStat, from = c(0,1), clip = TRUE )
                            )
              )

  vp = vp + geom_text() + scale_colour_gradientn(name = "Specificity Score",space = "rgb",
                                   colours = c("black", "white", "red"),
                                   values = c(0,1.3, 2)/2,
                                   breaks = c(0,1.3, 2)/2,
                                   labels = c(0, 1.3, 2),
                                   limits = c(0,1))

  vp = vp + scale_alpha("Consistency", limits = c(0,1))
  vp = vp + scale_size("Number of peptides")
  vp = vp + theme_bw() + theme(panel.background = element_rect(fill = 'lightyellow', colour = 'black'))
  return(vp)
}

ubMelt = function(aList){
  aMelt = melt(aList$M)
  aMelt = subset(aMelt, value > 0)
  aMelt$value = aList$mxRank
  aMelt = aMelt[,1:3]
  colnames(aMelt) = c("Peptide.ID", "kinase.name", "mxrank")
  return(aMelt)
}

#' @export
makePerPeptidePlot = function(df, dbFrame, scanRank = NULL, minPScore = NULL){
  if(!is.null(minPScore)){
    dbFrame= subset(dbFrame, Database != "phosphoNET" |  Kinase_PKinase_PredictorVersion2Score > minPScore)
  }
  if(!is.null(scanRank)){
    dbFrame = subset(dbFrame, Kinase_Rank <= scanRank)
  }
  ixList = intersectById(dbFrame, df)
  dbFrame = ixList[[1]]
  df = ixList[[2]]
  dbFrame = dbFrame %>% group_by(ID) %>% dplyr::summarise(lowestRank = min(Kinase_Rank))
  if (!is.null(df[["grp"]])){
      perPepStats = df %>% group_by(ID) %>% do({
      thisPep = subset(dbFrame, ID == .$ID[1])
      aTest = t.test(value ~ grp, data = ., var.equal = TRUE)
      aResult = data.frame(pes = diff(aTest$estimate),
                           ciu = -aTest$conf.int[1],
                           cil = -aTest$conf.int[2],
                           lowRank = thisPep$lowestRank
      )
    })
    ppp = ggplot(perPepStats, aes(x = reorder(ID, -lowRank), colour = as.factor(lowRank), y = pes, ymin = cil, ymax = ciu)) + geom_point() + geom_errorbar()
    ppp = ppp + coord_flip()
    ppp = ppp + xlab("Peptide ID") + ylab("Group difference")
    ppp = ppp + geom_abline(slope = 0)
    ppp = ppp + guides( color = guide_legend("Lowest Kinase Rank"))

  } else {
    perPepStats = df %>% group_by(ID) %>% do({
      thisPep = subset(dbFrame, ID == .$ID[1])
      aResult = data.frame(value = .$value,
                           lowRank = thisPep$lowestRank
      )
    })
    ppp = ggplot(perPepStats, aes(x = reorder(ID, -lowRank), colour = as.factor(lowRank), y = value)) + geom_point()
    ppp = ppp + coord_flip()
    ppp = ppp + xlab("Peptide ID") + ylab("Value")
    ppp = ppp + geom_abline(slope = 0)
    ppp = ppp + guides( color = guide_legend("Lowest Kinase Rank"))
  }
  yRange = layer_scales(ppp)$y$range$range

  ppp = ppp + scale_y_continuous(position = "top", limits = c(-max(abs(yRange)), max(abs(yRange))))
  return(ppp)
}

#' @export
makeDetailsTable = function(df, dbFrame, scanRank = NULL, minPScore = NULL){
  if (!is.null(minPScore)){
    dbFrame= subset(dbFrame, Database != "phosphoNET" |  Kinase_PKinase_PredictorVersion2Score > minPScore)
  }
  if (!is.null(scanRank)){
    dbFrame = subset(dbFrame, Kinase_Rank <= scanRank)
  }
  ixList = intersectById(dbFrame, df)
  dbFrame = ixList[[1]]
  outFrame = dbFrame %>% group_by(ID, PepProtein_PhosLink, PepProtein_UniprotName) %>%
    dplyr::summarise(Database = paste(Database, collapse = " / "), Kinase_Rank = min(Kinase_Rank) )
  outFrame = outFrame %>% arrange(Kinase_Rank, -as.integer(ID))
  return(outFrame)
}

#' @export
makeSummary = function(df){
  aSum = df %>% group_by(ClassName) %>% dplyr::summarise(meanFeatScore = mean(pFeatureScore, na.rm = TRUE),
                                                           meanPhenoScore = mean(pPhenoScore, na.rm = TRUE),
                                                           meanScore = mean(combinedScore, na.rm = TRUE),
                                                           medianScore = median(combinedScore, na.rm = TRUE),
                                                           meanStat = mean(NormalizedSetStat, na.rm = TRUE),
                                                           medianStat = median(NormalizedSetStat, na.rm = TRUE),
                                                           sdStat   = sd(NormalizedSetStat, na.rm = TRUE),
                                                           meanSetSize =mean(nFeatures, na.rm = TRUE))
  aSum = aSum%>% arrange(-medianScore)
  return(aSum)
}

#' @export
getSettingsInfo = function(settings){
  return(settings)
}


