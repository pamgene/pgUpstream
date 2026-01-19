#' Function for Confidence score based UKA
#' Performs UKA for both comparison and fold-change analysis types
#'
#' @param df long format data frame with columns ID, colSeq, value, and grp (latter only for comparison mode)
#' @param dbFrame data frame with the UKA database
#' @param scores_df data frame with C-scores for normalization, containing columns Database, Kinase_Rank, cscore
#' @param nPermutations number of permutations to use, default 500
#' @param norm_func normalization function, default "log10"
#' @param ukaType analysis type, either "Comparison" for two-group comparison or "FC" for fold-change analysis
#' @return a data frame with UKA results ordered by combined score
#' @import pgFCS
#' @import reshape2
#' @export
pgCscoreAnalysis <- function(df, dbFrame, scores_df,
                             nPermutations = 500,
                             norm_func = "log10",
                             ukaType) {
  # unified function for both comparison and FC analysis
  dbFrame <- cscore00(db = dbFrame, scores_df = scores_df, norm_func = norm_func)
  # only keep those IDs that are found in both db and data
  ixList <- intersectById(dbFrame, df)
  dbFrame <- ixList[[1]]
  df <- ixList[[2]]
  # Always use colSeq for data matrix construction
  X <- acast(df, colSeq ~ ID, fun.aggregate = mean, value.var = "value")
  M <- acast(dbFrame, ID ~ Kinase_Name, value.var = "wn", fun.aggregate = max)
  M[!is.finite(M)] <- 0
  inx2 <- intersect(colnames(X), rownames(M))
  Xi <- X[, colnames(X) %in% inx2]
  M <- M[rownames(M) %in% inx2, ]
  M <- M[order(rownames(M)), ]
  if (!all(rownames(M) == colnames(Xi))) stop("Mismatch between data matrix and upstream kinase matrix")
  # run functional scoring function with Comparison mode (with group) or FC mode
  if (ukaType == "Comparison") {
    Xi <- Xi[, order(colnames(Xi))] # Xi is a matrix
    grp <- acast(df, colSeq ~ ID, value.var = "grp")[, 1]
    grp <- as.factor(grp)
    result_df <- fcs(Xi, M, grp, nPerms = nPermutations)
  } else if (ukaType == "FC") {
    Xi <- Xi[order(names(Xi)), drop = FALSE] # Xi is a vector
    result_df <- fcs(Xi, M, statFun = stat.identity, phenoGrp = NULL, phenoPerms = FALSE, nPerms = nPermutations)
  } else {
    stop("ukaType must be either 'Comparison' or 'FC'")
  }

  result_df <- result_df[order(result_df$combinedScore, decreasing = TRUE), ]

  return(result_df)
}


#' perform a UKA analysis while scanning the PNET rank of upstream kinases
#' Intended for use with BN upstream apps
#' pgScanAnalysis2g(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis using binary grouping based on the grp column in df, implements parallel processing
#' pgScanAnalysis0(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis based on a single column (LFC) in df, no grouping, implements parallel processing
#' pgScanAnalysis2g_np(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis using binary grouping based on the grp column in df, no parallel processing
#' pgScanAnalysis0_np(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis based on a single column (LFC) in df, no grouping, no parallel processing
#'
#' @param df long format data frame with columns ID, colSeq, value, grp
#' @param dbFrame data frame with the UK database
#' @param dbWeights named vector of weights for each database, must match the names in dbFrame
#' @param scanRank vector of ranks to scan
#' @param nPermutations number of permutations to use
#' @param bootstrap logical, if TRUE, a bootstrap sample of the observations is used.
#' @return a list of data frames, one for each scanRank
#' @import pgFCS
#' @import dplyr
#' @import plyr
#' @import ggplot2
#' @import reshape2
#' @import foreach
#' @import doParallel
#' @import pgscales
#' @import data.table
#' @export
pgScanAnalysis2g <- function(df, dbFrame,
                             dbWeights,
                             scanRank,
                             nPermutations = 500,
                             bootstrap = FALSE) {
  # run two group.
  # add dbWeight
  dbFrame <- dbFrame %>%
    addDbWeights(dbWeights)

  ixList <- intersectById(dbFrame, df)
  dbFrame <- ixList[[1]]
  df <- ixList[[2]]
  X <- acast(df, colSeq ~ ID, fun.aggregate = mean, value.var = "value")
  grp <- acast(df, colSeq ~ ID, value.var = "grp")[, 1]
  if (bootstrap) {
    bIdx <- sample(1:length(grp), replace = TRUE)
    grp <- grp[bIdx]
    X <- X[bIdx, ]
  }
  grp <- as.factor(grp)
  nCores <- detectCores()
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  # print(paste("Please wait, processing on ", nCores, " CPU cores ..."))
  aPackagesList <- c("pgFCS", "reshape2")

  aScanResult <- foreach(i = scanRank, .packages = aPackagesList) %dopar% {
    aTop <- subset(dbFrame, Kinase_Rank <= i)

    M <- acast(aTop, ID ~ Kinase_Name, value.var = "dbWeight", fun.aggregate = max)
    M[!is.finite(M)] <- 0
    inx2 <- intersect(colnames(X), rownames(M))
    Xi <- X[, colnames(X) %in% inx2]
    M <- M[rownames(M) %in% inx2, ]
    # cannot depend on dcast for correct order
    M <- M[order(rownames(M)), ]
    Xi <- Xi[, order(colnames(Xi))]
    if (!all(rownames(M) == colnames(Xi))) stop("Mismatch between data matrix and upstream kinase matrix")
    aResult <- fcs(Xi, M, grp, nPerms = nPermutations)
    aResult <- aResult[order(aResult$combinedScore, decreasing = TRUE), ]
    return(data.table(mxRank = i, aResult = list(aResult), X = list(Xi), M = list(M)))
  }
  stopCluster(cl)
  return(aScanResult)
}
#' perform a UKA analysis while scanning the PNET rank of upstream kinases
#' Intended for use with BN upstream apps
#' pgScanAnalysis2g(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis using binary grouping based on the grp column in df, implements parallel processing
#' pgScanAnalysis0(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis based on a single column (LFC) in df, no grouping, implements parallel processing
#' pgScanAnalysis2g_np(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis using binary grouping based on the grp column in df, no parallel processing
#' pgScanAnalysis0_np(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis based on a single column (LFC) in df, no grouping, no parallel processing
#'
#' @param df long format data frame with columns ID, colSeq, value, grp
#' @param dbFrame data frame with the UK database
#' @param dbWeights named vector of weights for each database, must match the names in dbFrame
#' @param scanRank vector of ranks to scan
#' @param nPermutations number of permutations to use
#' @return a list of data frames, one for each scanRank
#' @export
pgScanAnalysis0 <- function(df, dbFrame,
                            dbWeights,
                            scanRank,
                            nPermutations = 500) {
  # run a sinlge column, without grouping
  # add dbWeight
  dbFrame <- dbFrame %>%
    addDbWeights(dbWeights)

  ixList <- intersectById(dbFrame, df)
  dbFrame <- ixList[[1]]
  df <- ixList[[2]]

  nCores <- detectCores()
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  # print(paste("Please wait, processing on ", nCores, " CPU cores ..."))
  aPackagesList <- c("pgFCS", "reshape2")

  aScanResult <- foreach(i = scanRank, .packages = aPackagesList) %dopar% {
    aTop <- subset(dbFrame, Kinase_Rank <= i)
    M <- acast(aTop, ID ~ Kinase_Name, value.var = "dbWeight", fun.aggregate = max)
    M[!is.finite(M)] <- 0
    inx2 <- intersect(df$ID, rownames(M))
    dfx <- df %>% filter(ID %in% inx2)
    X <- matrix(ncol = dim(dfx)[1], nrow = 1, data = dfx$value)
    colnames(X) <- dfx$ID
    M <- M[rownames(M) %in% inx2, ]
    M <- M[order(rownames(M)), ]
    X <- X[, order(colnames(X))]
    if (!all(rownames(M) == colnames(X))) stop("Mismatch between data matrix and upstream kinase matrix")
    aResult <- fcs(X, M, statFun = stat.identity, phenoGrp = NULL, phenoPerms = FALSE, nPerms = nPermutations)
    aResult <- aResult[order(aResult$combinedScore, decreasing = TRUE), ]
    return(data.table(mxRank = i, aResult = list(aResult), X = list(X), M = list(M)))
  }
  stopCluster(cl)
  return(aScanResult)
}

#' perform a UKA analysis while scanning the PNET rank of upstream kinases
#' Intended for use with BN upstream apps
#' pgScanAnalysis2g(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis using binary grouping based on the grp column in df, implements parallel processing
#' pgScanAnalysis0(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis based on a single column (LFC) in df, no grouping, implements parallel processing
#' pgScanAnalysis2g_np(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis using binary grouping based on the grp column in df, no parallel processing
#' pgScanAnalysis0_np(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis based on a single column (LFC) in df, no grouping, no parallel processing
#'
#' @param df long format data frame with columns ID, colSeq, value, grp
#' @param dbFrame data frame with the UK database
#' @param dbWeights named vector of weights for each database, must match the names in dbFrame
#' @param scanRank vector of ranks to scan
#' @param nPermutations number of permutations to use
#' @return a list of data frames, one for each scanRank
#' @export
pgScanAnalysis2g_np <- function(df, dbFrame,
                                dbWeights,
                                scanRank,
                                nPermutations = 500) {
  # run two group, non parallel
  # add dbWeight
  dbFrame <- dbFrame %>%
    addDbWeights(dbWeights)

  ixList <- intersectById(dbFrame, df)
  dbFrame <- ixList[[1]]
  df <- ixList[[2]]
  X <- acast(df, colSeq ~ ID, fun.aggregate = mean, value.var = "value")
  grp <- acast(df, colSeq ~ ID, value.var = "grp")[, 1]
  grp <- as.factor(grp)

  result <- lapply(scanRank, FUN = function(i) {
    aTop <- dbFrame %>%
      filter(Kinase_Rank <= i)

    M <- acast(aTop, ID ~ Kinase_Name, value.var = "dbWeight", fun.aggregate = max)
    M[!is.finite(M)] <- 0
    inx2 <- intersect(colnames(X), rownames(M))
    Xi <- X[, colnames(X) %in% inx2]
    M <- M[rownames(M) %in% inx2, ]
    # cannot depend on dcast for correct order
    M <- M[order(rownames(M)), ]
    Xi <- Xi[, order(colnames(Xi))]
    if (!all(rownames(M) == colnames(Xi))) stop("Mismatch between data matrix and upstream kinase matrix")
    aResult <- fcs(Xi, M, grp, nPerms = nPermutations)
    aResult <- aResult[order(aResult$combinedScore, decreasing = TRUE), ]
    return(data.table(mxRank = i, aResult = list(aResult), X = list(Xi), M = list(M)))
  })
  result
}

#' perform a UKA analysis while scanning the PNET rank of upstream kinases
#' Intended for use with BN upstream apps
#' pgScanAnalysis2g(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis using binary grouping based on the grp column in df, implements parallel processing
#' pgScanAnalysis0(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis based on a single column (LFC) in df, no grouping, implements parallel processing
#' pgScanAnalysis2g_np(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis using binary grouping based on the grp column in df, no parallel processing
#' pgScanAnalysis0_np(df, dbFrame, dbWeights, scanRank, nPermutations = 500) UKA analysis based on a single column (LFC) in df, no grouping, no parallel processing
#'
#' @param df long format data frame with columns ID, colSeq, value, grp
#' @param dbFrame data frame with the UK database
#' @param dbWeights named vector of weights for each database, must match the names in dbFrame
#' @param scanRank vector of ranks to scan
#' @param nPermutations number of permutations to use
#' @return a list of data frames, one for each scanRank
#' @export
pgScanAnalysis0_np <- function(df, dbFrame,
                               dbWeights,
                               scanRank,
                               nPermutations = 500) {
  # run a sinlge column, without grouping, non parallel
  # add dbWeight
  dbFrame <- dbFrame %>%
    addDbWeights(dbWeights)

  ixList <- intersectById(dbFrame, df)
  dbFrame <- ixList[[1]]
  df <- ixList[[2]]

  result <- lapply(scanRank, FUN = function(i)
  # aScanResult = for(i in scanRank)#debug
  {
    aTop <- subset(dbFrame, Kinase_Rank <= i)
    M <- acast(aTop, ID ~ Kinase_Name, value.var = "dbWeight", fun.aggregate = max)
    M[!is.finite(M)] <- 0
    inx2 <- intersect(df$ID, rownames(M))
    dfx <- df %>% filter(ID %in% inx2)
    X <- matrix(ncol = dim(dfx)[1], nrow = 1, data = dfx$value)
    colnames(X) <- dfx$ID
    M <- M[rownames(M) %in% inx2, ]
    M <- M[order(rownames(M)), ]
    X <- X[, order(colnames(X))]
    if (!all(rownames(M) == colnames(X))) stop("Mismatch between data matrix and upstream kinase matrix")
    aResult <- fcs(X, M, statFun = stat.identity, phenoGrp = NULL, phenoPerms = FALSE, nPerms = nPermutations)
    aResult <- aResult[order(aResult$combinedScore, decreasing = TRUE), ]
    return(data.table(mxRank = i, aResult = list(aResult), X = list(X), M = list(M)))
  })
  result
}


#' Find intersection of peptide IDs between two data frames
#' Helper function to intersect data frames by ID column and ensure matching peptide IDs
#'
#' @param df1 first data frame with ID column
#' @param df2 second data frame with ID column
#' @return a list containing two data frames filtered to common IDs
intersectById <- function(df1, df2) {
  df1$ID <- droplevels(df1$ID)
  df2$ID <- droplevels(df2$ID)
  isct <- intersect(df1$ID, df2$ID)
  if (length(isct) == 0) {
    stop("Error: Kinase family (PTK / STK) is incorrectly selected!")
  }
  df1 <- subset(df1, ID %in% isct)
  df2 <- subset(df2, ID %in% isct)
  df1$ID <- droplevels(df1$ID)
  df2$ID <- droplevels(df2$ID)
  return(list(df1, df2))
}

#' Add database weights to upstream kinase database
#' Helper function to add database-specific weights to the kinase database
#'
#' @param db data frame with upstream kinase database
#' @param dbWeights named vector of weights for each database, must match the Database column in db
#' @return the input database with dbWeight column added
#' @import dplyr
addDbWeights <- function(db, dbWeights) {
  dfWeights <- dbWeights %>%
    as.data.frame() %>%
    setNames("dbWeight") %>%
    mutate(Database = names(dbWeights))

  db <- db %>%
    left_join(dfWeights, by = "Database")

  if (any(is.na(db$dbWeight))) {
    stop("Missing or incorrect database weights!")
  }
  return(db)
}
