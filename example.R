## example for "differential" paper ##
library(tidyverse)
library(haven)
library(resample)
source("code/example_utils.R") 

# read in data
dat <- read_sas("/Users/jkedwar/OneDrive - University of North Carolina at Chapel Hill/data/CBHIPP/deidc1june2017.sas7bdat")

esdat <- dat %>% 
  mutate(male = ifelse(c14 == 1, 1, 
                       ifelse(c14 == 2, 0, NA)), 
         edu = ifelse(c15 == 999, NA, c15), 
         drinking = ifelse(c20b ==1, 1, ifelse(c20b == 2, 0, NA)), 
         selfreporthiv = ifelse((c62 == 1 & (c63b == 1)) | 
                                  #c67 == 1 |
                                  s12 == 1, 
                                1, 0), 
         hiv = ifelse(c117b == 1, 1, 0), 
         selfreportsw = ifelse(c50 == 1 |c50 == 2, 1, 0), 
         sw = ifelse(s2 == 1, 1, 0), 
         r = ifelse(is.na(hiv), 0, 1), 
         age = c12) 

nrow(esdat)
tab1 <- table(esdat$hiv, esdat$selfreporthiv)
table(esdat$hiv)

tab1
prop.table(tab1, margin = 1)


##' a function to get se and sp
getsesp <- function(data, differential = T){
  
  if(differential == T){
    #over 40
    se1 <- mean(data[data$y == 1&   data$z==1,]$ystar, na.rm = T)
    sp1 <- 1-mean(data[data$y == 0& data$z==1,]$ystar, na.rm = T)
    
    #under 40
    se2 <- mean(data[data$y == 1&   data$z==0,]$ystar, na.rm = T)
    sp2 <- 1-mean(data[data$y == 0& data$z==0,]$ystar, na.rm = T)
    return(rbind(c(se1, sp1), c(se2, sp2)))
  }
  if(differential == F){
    se1 <- mean(data[data$y == 1,]$ystar, na.rm = T)
    sp1 <- 1-mean(data[data$y == 0,]$ystar, na.rm = T)
    return(c(se1, sp1))
  }

}

sample <- esdat %>% filter(r == 1, !is.na(selfreporthiv)) %>% mutate(old = age>30)
sample <- sample %>% 
  mutate(y = hiv, ystar = selfreporthiv, z= as.integer(old)) %>% 
  select(y, ystar, z)

table(sample$y)
table(sample$ystar)
nrow(sample)

# create study sample


getsesp(sample, differential = T)
getsesp(sample, differential = F)

table(sample$ystar, sample$y, sample$z)
nrow(sample)
table(sample$ystar)

### Table 1 -----
tab0 <- table(sample[sample$z==0,]$ystar, sample[sample$z==0,]$y)
tab1 <- table(sample[sample$z==1,]$ystar, sample[sample$z==1,]$y)
prop.table(tab0, margin = 2)
prop.table(tab1, margin = 2)
tab <- table(sample$ystar, sample$y)
prop.table(tab, margin = 2)

# draw validation data ----

set.seed(555) #432
val1 <- sample %>% 
  slice_sample(n = 5000, replace = F) %>% 
  select(y, ystar, z)
table(val1$y, val1$ystar, val1$z)
getsesp(val1, differential = T)
getsesp(val1, differential = F)

tab0 <- table(val1[val1$z==0,]$ystar, val1[val1$z==0,]$y)
tab1 <- table(val1[val1$z==1,]$ystar, val1[val1$z==1,]$y)
prop.table(tab0, margin = 2)
prop.table(tab1, margin = 2)
tab <- table(val1$ystar, val1$y)
prop.table(tab, margin = 2)


set.seed(3111) #42 33#235
val2 <- sample %>% 
  mutate(wt = ifelse(z == 1, 10, .5)) %>% #dependent on age
  slice_sample(n = 5000, replace = F, weight_by = wt) %>% 
  select(y, ystar, z)
table(val2$y, val2$ystar, val2$z)

getsesp(val2, differential = T)
getsesp(val2, differential = F)
table(val2$y, val2$ystar)

tab0 <- table(val2[val2$z==0,]$ystar, val2[val2$z==0,]$y)
tab1 <- table(val2[val2$z==1,]$ystar, val2[val2$z==1,]$y)
prop.table(tab0, margin = 2)
prop.table(tab1, margin = 2)
tab <- table(val2$ystar, val2$y)
prop.table(tab, margin = 2)

#table 2 (results) ----
tab <- table(sample$z)
prop.table(tab)

pztarget = .85

# true
true_sample <- mu0(sample, target = "sample", pztarget = NA)
ci <- ci_naive(true_sample, nrow(sample))
true_sample_ci <- paste0(round(true_sample*100, 1), " (", ci, ")")

true_target <- mu0(sample, target = "external", pztarget = pztarget)
ci <- ci_naive(true_target, nrow(sample))
true_target_ci <- paste0(round(true_target*100, 1), " (", ci, ")")

#naive
naive_sample <- mu1(sample, target = "sample", pztarget = NA)
ci <- ci_naive(naive_sample, nrow(sample))
naive_sample_ci <- paste0(round(naive_sample*100, 1), " (", ci, ")")
naive_target <- mu1(sample, target = "external", pztarget = pztarget)
ci <- ci_naive(naive_target, nrow(sample))
naive_target_ci <- paste0(round(naive_target*100, 1), " (", ci, ")")

# apply correction nd val 1
sesp <- getsesp(val1, differential = T)
nd_val1_sample <- mu2(sample, val1, target = "sample", pztarget = NA)
nd_val1_target <- mu2(sample, val1, target = "external", pztarget = pztarget)

# nd_var_val1_sample <- var_nd(val1, nd_val1_sample, naive_sample, nrow(sample))
# se1 <- sqrt(nd_var_val1_sample)
se1 <- sd(boot(500, sample, val1, mu2, target = "sample"))
lcl <- nd_val1_sample - 1.96 * se1
ucl <- nd_val1_sample + 1.96 * se1
ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
nd_val1_sample_ci <- paste0(round(nd_val1_sample*100, 1), " (", ci, ")")

se1 <- sd(boot(500, sample, val1, mu2, target = "external", pztarget = pztarget))
lcl <- nd_val1_target - 1.96 * se1
ucl <- nd_val1_target + 1.96 * se1
ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
nd_val1_target_ci <- paste0(round(nd_val1_target*100, 1), " (", ci, ")")


# apply correction nd val2
sesp <- getsesp(val2, differential = T)
nd_val2_sample <- mu2(sample, val2, target = "sample", pztarget = NA)
nd_val2_target <- mu2(sample, val2, target = "external", pztarget = pztarget)

se1 <- sd(boot(500, sample, val2, mu2, target = "sample"))
lcl <- nd_val2_sample - 1.96 * se1
ucl <- nd_val2_sample + 1.96 * se1
ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
nd_val2_sample_ci <- paste0(round(nd_val2_sample*100, 1), " (", ci, ")")

se1 <- sd(boot(500, sample, val2, mu2, target = "external", pztarget = pztarget))
lcl <- nd_val2_target - 1.96 * se1
ucl <- nd_val2_target + 1.96 * se1
ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
nd_val2_target_ci <- paste0(round(nd_val2_target*100, 1), " (", ci, ")")


#apply diff correction
sesp <- getsesp(val1, differential = T)
d_val1_sample <- mu3(sample, val1, target = "sample", pztarget = NA)
d_val1_target <- mu3(sample, val1, target = "external", pztarget = pztarget)

se1 <- sd(boot(500, sample, val1, mu3, target = "sample"))
lcl <- d_val1_sample - 1.96 * se1
ucl <- d_val1_sample + 1.96 * se1
ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
d_val1_sample_ci <- paste0(round(d_val1_sample*100, 1), " (", ci, ")")

se1 <- sd(boot(500, sample, val1, mu3, target = "external", pztarget = pztarget))
lcl <- d_val1_target - 1.96 * se1
ucl <- d_val1_target + 1.96 * se1
ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
d_val1_target_ci <- paste0(round(d_val1_target*100, 1), " (", ci, ")")

# apply correction nd val2
sesp <- getsesp(val2, differential = T)
d_val2_sample <- mu3(sample, val2, target = "sample", pztarget = NA)
d_val2_target <- mu3(sample, val2, target = "external", pztarget = pztarget)

se1 <- sd(boot(500, sample, val2, mu3, target = "sample"))
lcl <- d_val1_sample - 1.96 * se1
ucl <- d_val1_sample + 1.96 * se1
ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
d_val2_sample_ci <- paste0(round(d_val2_sample*100, 1), " (", ci, ")")

se1 <- sd(boot(500, sample, val2, mu3, target = "external", pztarget = pztarget))
lcl <- d_val1_target - 1.96 * se1
ucl <- d_val1_target + 1.96 * se1
ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
d_val2_target_ci <- paste0(round(d_val2_target*100, 1), " (", ci, ")")

r1 <- c(rep(true_sample_ci, 2), rep(true_target_ci, 2))
r2 <- c(rep(naive_sample_ci, 2), rep(naive_target_ci, 2))
r3 <- c(nd_val1_sample_ci, nd_val2_sample_ci, nd_val1_target_ci,  nd_val2_target_ci)
r4 <- c(d_val1_sample_ci, d_val2_sample_ci, d_val1_target_ci,  d_val2_target_ci)

restable <- rbind(r1, r2, r3, r4)
write.csv(restable, file = "output/exres.csv")
