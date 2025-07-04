library(dplyr)
rank2cscore = function(rank, rank_table){
  score = rep(rank_table$score[1])
  for(i in 2:nrow(rank_table)){
    score = ifelse(rank >= rank_table$rank[i], rank_table$score[i], score)
  }
  return(score)
}

rtable = data.frame( rank = c(0,1,6,12), score = c(0.9, 0.7, 0.5, 0))

mr = data.frame(ranks = c(0,0,0,1,1,1,6,7,8,9,10,11,12,13))

mr = mr %>%
  mutate(c = rank2cscore(ranks, rtable))



