### run sim plot ###

library(tidyverse)
library(ggthemr)
library(cowplot)

# read in file with helper functions ----
source("code/utils.R") 
today <- "8feb22"

# set up common parameters ----
set.seed(234)
n <- 8000 #sample size
J <- 10000 #number of sims
valsize = 2000

true <- 0.1 #true p(y) in study sample
pzvals <- c(.15,  .25,  .35, .45, .55, .65, .75, .85)

# run sims ----
sc0res <- runfigsc(0, n, valsize, J, pzvals)
sc1res <- runfigsc(1, n, valsize, J, pzvals)
sc2res <- runfigsc(2, n, valsize, J, pzvals)
sc3res <- runfigsc(3, n, valsize, J, pzvals)
sc4res <- runfigsc(4, n, valsize, J, pzvals)

## table ----


tab1 <- sc1res %>% filter(pz %in% c(0.35, 0.55, 0.75)) %>% mutate(sc = 1) 
tab2 <- sc2res %>% filter(pz %in% c(0.35, 0.55, 0.75)) %>% mutate(sc = 2)
tab3 <- sc3res %>% filter(pz %in% c(0.35, 0.55, 0.75)) %>% mutate(sc = 3)
tab4 <- sc4res %>% filter(pz %in% c(0.35, 0.55, 0.75)) %>% mutate(sc = 4)
tab <- bind_rows(tab1, tab2, tab3, tab4) %>% 
  select(pz, target, estimator, mu, bias, se, rmse, sc) %>% 
  mutate(bias = round(bias*100, 2), mu = round(mu, 2), se = round(se * 100, 2), 
         rmse = round(sqrt(bias^2 + se^2), 2)) %>% 
  pivot_wider(id_cols = c(sc, pz, estimator), names_from = target, values_from = c(mu, bias, se, rmse)) %>%
  filter(estimator %in% c(2,3) | pz == 0.35)   #limit to rows i want in table
  
  
tab <- tab[, c(1,2,3,4,6,8,10,5,7,9,11)]  
write.csv(tab, file = paste0("output/simres_", today, ".csv"))

## appendix table 1 ----
tab0 <- sc0res %>% filter(pz %in% c(0.35, 0.55, 0.75)) %>% mutate(sc = 0) 
tab0 <- tab0 %>% 
  select(pz, target, estimator, mu, bias, se, rmse, sc) %>% 
  mutate(bias = round(bias*100, 2), mu = round(mu, 2), se = round(se * 100, 2), 
         rmse = round(sqrt(bias^2 + se^2), 2)) %>% 
  pivot_wider(id_cols = c(sc, pz, estimator), 
              names_from = target, values_from = c(mu, bias, se, rmse)) #%>%
write.csv(tab0, file = paste0("output/simres_tab0_", today, ".csv"))

## figure 2 ----
ggthemr("fresh")
pl1 <- plotfun(sc1res)
pl2 <- plotfun(sc2res)
pl3 <- plotfun(sc3res)
pl4 <- plotfun(sc4res)
mat <- plot_grid(pl1, pl2, pl3, pl4, labels = c("A", "B", "C", "D"), nrow  = 2)
ggsave(mat, file = paste0("output/simfig_allv5b_", today, ".png"), width = 30, height = 20, units = "cm")


## figure 3 ----
mse1 <- msefun(sc1res)
mse2 <- msefun(sc2res)
mse3 <- msefun(sc3res)
mse4 <- msefun(sc4res)
mat2 <- plot_grid(mse1, mse2, mse3, mse4, labels = c("A", "B", "C", "D"), nrow  = 2)
ggsave(mat2, file = paste0("output/simfig_mseb_", today, ".png"), width = 30, height = 20, units = "cm")

# outcome dependent sampling ----

valsize = 200

true <- 0.1 #true p(y) in study sample

## run sims ----
sc0res_od <- runfig_od(0, n, valsize, J)
sc1res_od <- runfig_od(1, n, valsize, J)
sc2res_od <- runfig_od(2, n, valsize, J)
sc3res_od <- runfig_od(3, n, valsize, J)
sc4res_od <- runfig_od(4, n, valsize, J)

## table for outcome dep sampling ----
tab1 <- sc1res_od %>% mutate(sc = 1) 
tab2 <- sc2res_od %>% mutate(sc = 2)
tab3 <- sc3res_od %>% mutate(sc = 3)
tab4 <- sc4res_od %>% mutate(sc = 4)
tab <- bind_rows(tab1, tab2, tab3, tab4) %>% 
  select(val, target, estimator, mu, bias, se, rmse, sc) %>% 
  mutate(bias = round(bias*100, 2), mu = round(mu, 2), se = round(se * 100, 2), 
         rmse = round(sqrt(bias^2 + se^2), 2)) %>% 
  pivot_wider(id_cols = c(sc, val, estimator), names_from = target, values_from = c(mu, bias, se, rmse)) %>%
  filter(estimator %in% c(2,3) | val == 1)   #limit to rows i want in table


tab <- tab[, c(1,2,3,4,6,8,10,5,7,9,11)]  
write.csv(tab, file = paste0("output/simres_od_", today, ".csv"))
