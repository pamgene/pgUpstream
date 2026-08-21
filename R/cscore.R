# methods for c-scoring

#' Function for simple c-scoring on standard upstream DB
#' @param score_var which column to match scores_df against: "Kinase_Rank"
#'   (default) does an exact-match join on Database x Kinase_Rank via
#'   add_scores(); any other value is treated as a continuous score column
#'   (e.g. "Kinase_Score") matched against ascending per-Database thresholds
#'   in scores_df via add_scores_by_score()
#' @import dplyr
#' @export
cscore00 <- function(db, scores_df, norm_func, score_var = "Kinase_Rank") {
  db_dedup <- db %>%
    distinct(ID, Kinase_Name, PepProtein_PhosLink, Database, .keep_all = TRUE)

  db_scored <- if (score_var == "Kinase_Rank") {
    db_dedup %>% add_scores(scores_df = scores_df, add_by = c("Database", "Kinase_Rank"))
  } else {
    db_dedup %>% add_scores_by_score(scores_df = scores_df, score_col = score_var)
  }

  result <- db_scored %>%
    cscore2w(
      sub.var = "ID",
      up.var = "Kinase_Name",
      value.var = "cscore",
      na.wt = 0,
      norm_func = norm_func
    )

  return(result)
}

#' Function for calculating (normalized) weights based on c-scores
#' @import dplyr
#' @export
cscore2w <- function(dbframe, sub.var = "ID", up.var = "Kinase_Name", value.var = "cscore", na.wt = 0, norm_func) {
  dbframe$cvar <- dbframe %>%
    pull(value.var)
  dbframe$up.var <- dbframe %>%
    pull(up.var)
  dbframe$sub.var <- dbframe %>%
    pull(sub.var)
  dbw <- dbframe %>%
    group_by(sub.var, up.var) %>%
    dplyr::summarise(.w = 1 - prod(1 - cvar)) %>%
    ungroup() %>%
    pivot_wider(id_cols = sub.var, names_from = up.var, values_from = ".w") %>%
    pivot_longer(-sub.var, values_to = ".w", names_to = up.var) %>%
    mutate(w = ifelse(is.na(.w), na.wt, .w))

  subsum <- dbw %>%
    group_by(sub.var) %>%
    dplyr::summarise(
      sup_count = sum(!is.na(.w)),
      sup_mean_weight = mean(.w, na.rm = TRUE),
      sup_sum_weight = sum(w)
    ) %>%
    ungroup()

  if (norm_func == "sumw") {
    dbw <- dbw %>%
      left_join(subsum, by = "sub.var") %>%
      mutate(wn = w / sup_sum_weight)
  } else if (norm_func == "log10") {
    dbw <- dbw %>%
      left_join(subsum, by = "sub.var") %>%
      mutate(wn = w / log(sup_sum_weight + 1, 10)) %>% # log for dampening normalization effect
      mutate(wn = pmin(wn, 1)) # cap at 1
  }
  colnames(dbw)[1:2] <- c(sub.var, up.var)
  return(dbw)
}


#' Check that every Database in scores_df is present in db, used by both
#' add_scores() and add_scores_by_score() before they attempt any matching
#' @noRd
validate_scores_df_databases <- function(db, scores_df) {
  dbs <- scores_df %>%
    pull(Database) %>%
    unique()
  dbs_found <- dbs %in% db$Database
  if (!all(dbs_found)) {
    stop(paste("Some databases in scores not found in upstream db, check spelling: ", dbs[!dbs_found]))
  }
}


#' Function for adding cscores to an Upstream DB using a list with database and rank columns
#' @import dplyr
#' @export
add_scores <- function(db, scores_df, add_by = c("Database", "Kinase_Rank")) {
  validate_scores_df_databases(db, scores_df)

  ndistinct <- scores_df %>%
    distinct(Database, Kinase_Rank) %>%
    nrow()

  if (ndistinct != nrow(scores_df)) {
    stop("Not all entries in scores have unique mapping")
  }

  db %>%
    left_join(scores_df, by = add_by)
}


#' Function for converting a rank metric to a c-score based on rank_table
#' @import dplyr
#' @export
rank2cscore <- function(rank, rank_table) {
  score <- rep(rank_table$score[1])
  for (i in 2:nrow(rank_table)) {
    score <- ifelse(rank >= rank_table$rank[i], rank_table$score[i], score)
  }
  return(score)
}


#' Function for converting a per-Database continuous score metric to a
#' c-score using ascending threshold breakpoints within each Database. A
#' row's cscore is that of the highest threshold its score meets or exceeds,
#' for the Database it belongs to. A Database with a single threshold row
#' gets a flat cscore regardless of the actual score value (e.g. for
#' databases with no real score concept).
#' @param score numeric vector of raw per-row scores
#' @param database character vector of Database values, same length as score
#' @param score_table data.frame with columns Database, <score_col>, cscore
#' @param score_col name of the score column in score_table, default "Kinase_Score"
#' @import dplyr
#' @export
score2cscore <- function(score, database, score_table, score_col = "Kinase_Score") {
  result <- rep(NA_real_, length(score))
  for (db in unique(score_table$Database)) {
    tbl <- score_table[score_table$Database == db, , drop = FALSE]
    tbl <- tbl[order(tbl[[score_col]]), ]
    idx <- database == db
    if (!any(idx)) next
    if (nrow(tbl) == 1) {
      result[idx] <- tbl$cscore[1]
      next
    }
    s <- rep(tbl$cscore[1], sum(idx))
    for (i in 2:nrow(tbl)) {
      s <- ifelse(score[idx] >= tbl[[score_col]][i], tbl$cscore[i], s)
    }
    result[idx] <- s
  }
  return(result)
}


#' Function for adding cscores to an Upstream DB using a continuous
#' per-Database score metric (e.g. Kinase_Score), matched against ascending
#' thresholds instead of an exact Kinase_Rank match (see add_scores())
#' @param db data frame with the upstream database, must contain a Database
#'   column and the column named by score_col
#' @param scores_df data.frame with columns Database, <score_col>, cscore
#' @param score_col name of the score column to match on, default "Kinase_Score"
#' @import dplyr
#' @export
add_scores_by_score <- function(db, scores_df, score_col = "Kinase_Score") {
  validate_scores_df_databases(db, scores_df)

  db$cscore <- score2cscore(db[[score_col]], db$Database, scores_df, score_col)
  return(db)
}


#' Function for simple c-scoring on standard upstream DB
#' @import dplyr
#' @export
cscore0 <- function(db, rank_table = data.frame(rank = c(0, 1, 6, 12), score = c(0.9, 0.7, 0.5, 0))) {
  db <- db %>%
    distinct(ID, Kinase_Name, PepProtein_PhosLink, Database, .keep_all = TRUE) %>%
    mutate(cscore = rank_table$score[1])

  for (i in 2:nrow(rank_table)) {
    db <- db %>%
      mutate(cscore = ifelse(Kinase_Rank >= rank_table$rank[i], rank_table$score[i], cscore))
  }
  dbw <- db %>%
    group_by(ID, Kinase_Name) %>%
    dplyr::summarise(w = 1 - prod(1 - cscore)) %>%
    ungroup()

  pepcount <- dbw %>%
    filter(w > 0) %>%
    group_by(ID) %>%
    dplyr::summarise(
      ks_count = n(),
      ks_mean_weight = mean(w),
      ks_sum_weight = sum(w)
    ) %>%
    ungroup()

  dbw <- dbw %>%
    left_join(pepcount, by = "ID") %>%
    mutate(wn = w / ks_sum_weight)
}
