rm(list = ls())
library(pgUpstream)
library(tidyverse)
library(reshape2)

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
  cscore2w(na.wt = 0)

dbw = dbw %>%
  arrange(Kinase_Name) %>%
  filter(!is.na(wn)) %>%
  filter(wn > 0) %>%
  select(Kinase_Name, ID, wn)

data_in  = read.csv("FGFR_example_train_2.csv" )
X = data_in %>%
  acast(Sample.name ~ ID, value.var = "value")

grp = data_in %>%
  acast(Sample.name ~ ID, value.var = "group")

grp = grp[,1] %>%  as.factor()

gresult = gtUKA(X, grp, dbw)

gresult = gresult %>%
  filter(n_substrates > 2) %>%
  arrange(p_value)






