rm(list = ls())
library(pgUpstream)
library(tidyverse)


db = readRDS("250519_86312_86402_86412_87102_87202_UpstreamDb (1).RDS")

scores = list( data.frame(database = "iviv", rank = 0, cscore =0.9),
               expand.grid(database = "PhosphoNET", rank = 1:5) %>%  mutate(cscore = 0.7),
               expand.grid(database = "PhosphoNET", rank = 6:12) %>%  mutate(cscore = 0.5),
               expand.grid(database = "KL", rank = 1:5) %>%  mutate(cscore = 0.7),
               expand.grid(database = "KL", rank = 6:12) %>% mutate(cscore = 0) )

score_df =  scores %>%
  bind_rows()

db = db %>%
  filter(Kinase_group == "TYR")

nk = db %>%
  pull(Kinase_Name) %>%
  unique() %>%
  length()


dbw = db %>%
  add_scores(score_df) %>%
  cscore2w(na.wt = 1/nk)


dbw %>%
  filter(!is.na(.w)) %>%
  ggplot(aes(x = wn)) +
  geom_histogram(binwidth = 0.005)


dbw %>%
  #filter(!is.na(.w)) %>%
  ggplot(aes(x = Kinase_Name, y = reorder(ID, sup_sum_weight), fill = sqrt(wn))) +
  viridis::scale_fill_viridis() +
  geom_tile()





