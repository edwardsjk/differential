### utilities for differential misclassification paper ###

# data generating mechanisms ----

##' a function to generate data for scenarios 0-2

gen2 <- function(i, n, se1, se2, sp1, sp2, pz = 0.75){
  set.seed(i)
  y <- rbinom(n, 1, true)
  z <- rbinom(n, 1, pz)
  ystar <- ifelse(z == 1, 
                  ifelse(y == 1, rbinom(n, 1, se1), 
                         1 - rbinom(n, 1, sp1)), 
                  ifelse(y == 1, rbinom(n, 1, se2), 
                         1 - rbinom(n, 1, sp2)))
  dat1 <- tibble(y = y, ystar = ystar, z = z)
  return(dat1)
}


##' a function to generate data for scenarios 3 and 4
##' 
gen3 <- function(i,n, se1, se2, sp1, sp2, pz = 0.75){
  set.seed(i)
  z <- rbinom(n, 1, pz)
  py <- plogis(log(true/(1-true)) - log(3)*sd(z)^2/2 + log(3)*mean(z) - log(3)*z) #
  y <- rbinom(n, 1, py)
  ystar <- ifelse(z == 1, 
                  ifelse(y == 1, rbinom(n, 1, se1), 
                         1 - rbinom(n, 1, sp1)), 
                  ifelse(y == 1, rbinom(n, 1, se2), 
                         1 - rbinom(n, 1, sp2)))
  dat1 <- tibble(y = y, ystar = ystar, z = z)
  return(dat1)
}


# estimators ----


mu0 <- function(dat, target = "sample", pztarget = NA){
  if(target == "sample"){
    return(mean(dat$y))
  }
  if(target == "external"){
    return(mean(dat[dat$z == 0,]$y) * (1-pztarget) + mean(dat[dat$z==1,]$y) * pztarget)
  }
}

mu1 <- function(dat, target = "sample", pztarget = NA){
  if(target == "sample"){
    return(mean(dat$ystar))
  }
  if(target == "external"){
    return(mean(dat[dat$z == 0,]$ystar) * (1-pztarget) + mean(dat[dat$z==1,]$ystar) * pztarget)
  }
}

mu2 <- function(dat, val = "internal", target = "sample", pztarget = NA){
  if(val == "internal"){
    samp <- dat %>% sample_frac(0.2)
  }
  if(val == "target"){
    pzval <- 0.35
    prz1 <- pzval*.2/.75
    prz0 <- (0.2 - prz1*.75)/.25
    samp <- dat %>% 
      mutate(pr = ifelse(z == 0, prz0, prz1),  #gives p(z|val) = .35
             r = rbinom(nrow(dat), 1, pr)) %>% 
      filter(r == 1)#filter(z == 0) %>% sample_frac(0.536)
    # samp2 <- dat %>% filter(z == 1) %>% sample_frac(0.088)
    # samp <- bind_rows(samp1, samp2)
  }
  if(val == "other"){
    pzval <- 0.55
    prz1 <- pzval*.2/.75
    prz0 <- (0.2 - prz1*.75)/.25
    samp <- dat %>% 
      mutate(pr = ifelse(z == 0, prz0, prz1),  #gives p(z|val) = .55
             r = rbinom(nrow(dat), 1, pr)) %>% 
      filter(r == 1)

  }
  se <- mean(samp[samp$y==1,]$ystar)
  sp <- 1-mean(samp[samp$y==0,]$ystar)
  
  if(target == "sample"){
    A <- (dat$ystar - (1 - sp))/(se+sp-1)
    return(mean(A))
  }
  if(target == "external"){
    datz0 <- dat[dat$z==0,]
    datz1 <- dat[dat$z==1,]
    A1 <- (datz0$ystar - (1 - sp))/(se+sp-1)
    A2 <- (datz1$ystar -  (1 - sp))/(se+sp-1)
    return(sum(A1)/(nrow(dat[dat$z==0,])) *(1 - pztarget) + sum(A2)/(nrow(dat[dat$z==1,]))*pztarget)
  }
}

mu3 <- function(dat, val = "internal", target = "sample", pztarget = NA){
  if(val == "internal"){
    samp <- dat %>% sample_frac(0.2)
  }
  if(val == "target"){
    pzval <- 0.35
    prz1 <- pzval*.2/.75
    prz0 <- (0.2 - prz1*.75)/.25
    samp <- dat %>% 
      mutate(pr = ifelse(z == 0, prz0, prz1),  #gives p(z|val) = .35
             r = rbinom(nrow(dat), 1, pr)) %>% 
      filter(r == 1)
  }
  if(val == "other"){
    pzval <- 0.55
    prz1 <- pzval*.2/.75
    prz0 <- (0.2 - prz1*.75)/.25
    samp <- dat %>% 
      mutate(pr = ifelse(z == 0, prz0, prz1),  #gives p(z|val) = .55
             r = rbinom(nrow(dat), 1, pr)) %>% 
      filter(r == 1)
  }
  se1 <- mean(samp[samp$y==1 & samp$z==0,]$ystar)
  sp1 <- 1-mean(samp[samp$y==0& samp$z==0,]$ystar)
  se2 <- mean(samp[samp$y==1& samp$z==1,]$ystar)
  sp2 <- 1-mean(samp[samp$y==0& samp$z==1,]$ystar)
  A1 <- (dat$ystar*(dat$z==0) - (dat$z==0) * (1 - sp1))/(se1+sp1-1)
  A2 <- (dat$ystar*(dat$z==1) - (dat$z==1) * (1 - sp2))/(se2+sp2-1)
  if(target == "sample"){
    return(sum(A1+A2)/nrow(dat))
  }
  if(target == "external"){
    return(sum(A1)/(nrow(dat[dat$z==0,])) *(1 - pztarget) + sum(A2)/(nrow(dat[dat$z==1,]))*pztarget)
  }
}

# atlernate functions for plot ----


mu2_alt <- function(dat, valdata, target = "sample", pztarget = NA){
  
  se <- mean(valdata[valdata$y==1,]$ystar)
  sp <- 1-mean(valdata[valdata$y==0,]$ystar)
  
  if(target == "sample"){
    A <- (dat$ystar - (1 - sp))/(se+sp-1)
    return(mean(A))
  }
  if(target == "external"){
    datz0 <- dat[dat$z==0,]
    datz1 <- dat[dat$z==1,]
    A1 <- (datz0$ystar - (1 - sp))/(se+sp-1)
    A2 <- (datz1$ystar -  (1 - sp))/(se+sp-1)
    return(sum(A1)/(nrow(dat[dat$z==0,])) *(1 - pztarget) + sum(A2)/(nrow(dat[dat$z==1,]))*pztarget)
  }
}

##' a weighted version of mu2_alt
##' valdata must contain a variable pi with the weight

mu2_alt_w <- function(dat, valdata,  target = "sample", pztarget = NA){
  
  se <- valdata %>% filter(y == 1) %>% 
    summarize(se = sum(ystar * pi)/sum(pi)) %>% pull()
  sp <- valdata %>% filter(y == 0) %>% 
    summarize(se = sum((1 - ystar) * pi)/sum(pi)) %>% pull()
  
  if(target == "sample"){
    A <- (dat$ystar - (1 - sp))/(se+sp-1)
    return(mean(A))
  }
  if(target == "external"){
    datz0 <- dat[dat$z==0,]
    datz1 <- dat[dat$z==1,]
    A1 <- (datz0$ystar - (1 - sp))/(se+sp-1)
    A2 <- (datz1$ystar -  (1 - sp))/(se+sp-1)
    return(sum(A1)/(nrow(dat[dat$z==0,])) *(1 - pztarget) + sum(A2)/(nrow(dat[dat$z==1,]))*pztarget)
  }
}

mu3_alt <- function(dat, valdata, target = "sample", pztarget = NA){
  samp <- valdata
  se1 <- mean(samp[samp$y==1 & samp$z==0,]$ystar)
  sp1 <- 1-mean(samp[samp$y==0& samp$z==0,]$ystar)
  se2 <- mean(samp[samp$y==1& samp$z==1,]$ystar)
  sp2 <- 1-mean(samp[samp$y==0& samp$z==1,]$ystar)
  A1 <- (dat$ystar*(dat$z==0) - (dat$z==0) * (1 - sp1))/(se1+sp1-1)
  A2 <- (dat$ystar*(dat$z==1) - (dat$z==1) * (1 - sp2))/(se2+sp2-1)
  if(target == "sample"){
    return(sum(A1+A2)/nrow(dat))
  }
  if(target == "external"){
    return(sum(A1)/(nrow(dat[dat$z==0,])) *(1 - pztarget) + sum(A2)/(nrow(dat[dat$z==1,]))*pztarget)
  }
}


# simulation helper functions ----

##' second generation simulation function
##' 
runfig <- function(i, n, valn = 200, pzval = 0.5, tar = "sample", pztarget = NA, sc = 4){
  if(sc == 0){
    dat <- gen2(i, n, 1, 1, 1, 1, pz = 0.75)
    valdat <- gen2(i, valn, 1, 1, 1, 1, pz = pzval)
  }
  if(sc == 1){
    dat <- gen2(i, n, .7, .7, .9, .9, pz = 0.75)
    valdat <- gen2(i, valn, .7, .7, .9, .9, pz = pzval)
  }
  if(sc == 2){
    dat <- gen2(i, n, .7, .95, .95, .85, pz = 0.75)
    valdat <- gen2(i, valn, .7, .95, .95, .85, pz = pzval)
  }
  if(sc == 3){
    dat <- gen3(i, n, .7, .7, .9, .9, pz = 0.75)
    valdat <- gen3(i, valn, .7, .7, .9, .9, pz = pzval)
  }
  if(sc == 4){
    dat <- gen3(i, n, .7, .9, .95, .85, pz = 0.75)#ning: se2=0.9 here but 0.95 in paper, revise the paper is easier :)
    valdat <- gen3(i, valn, .7, .9, .95, .85, pz = pzval)#ning: the same issue for se2 as above
  }
  
  mu0 <- mu0(dat, target = tar, pztarget = pztarget)
  mu1 <- mu1(dat, target = tar, pztarget = pztarget)
  mu2 <- mu2_alt(dat, valdat, target = tar, pztarget = pztarget)
  mu3 <- mu3_alt(dat, valdat, target = tar, pztarget = pztarget)
  return(c(mu1, mu2, mu3, mu0))
}

##' a function to produce table of simulation results
##' 
runtab <- function(i, n, valn = 200, pzval = 0.5, tar = "sample", pztarget = NA, sc = 4){
  if(sc == 1){
    dat <- gen2(i, n, .7, .7, .9, .9, pz = 0.75)
    valdat <- gen2(i, valn, .7, .7, .9, .9, pz = pzval)
  }
  if(sc == 2){
    dat <- gen2(i, n, .7, .95, .95, .85, pz = 0.75)
    valdat <- gen2(i, valn, .7, .95, .95, .85, pz = pzval)
  }
  if(sc == 3){
    dat <- gen3(i, n, .7, .7, .9, .9, pz = 0.75)
    valdat <- gen3(i, valn, .7, .7, .9, .9, pz = pzval)
  }
  if(sc == 4){#ning: se2 issue
    dat <- gen3(i, n, .7, .9, .95, .85, pz = 0.75)
    valdat <- gen3(i, valn, .7, .9, .95, .85, pz = pzval)
  }
  mu0 <- mu0(dat, target = tar, pztarget = pztarget)#ning: revise the function from mu1 to mu0
  mu1 <- mu1(dat, target = tar, pztarget = pztarget)
  mu2 <- mu2_alt(dat, valdat, target = tar, pztarget = pztarget)
  mu3 <- mu3_alt(dat, valdat, target = tar, pztarget = pztarget)
  return(c(mu0, mu1, mu2, mu3))
}


##' function to do sims with outcome dep sampling
##' 
runfig2 <- function(i, n, valn = 200, valvar = y, valfrom = "sample", tar = "sample", 
                    pztarget = .5, sc = 4){
  if(sc == 0){
    dat <- gen2(i, n, 1, 1, 1, 1, pz = 0.75)
    valdat <- gen2(i, valn, 1, 1, 1, 1, pz = 0.75)
  }
  if(sc == 1){
    dat <- gen2(i, n, .7, .7, .9, .9, pz = 0.75)
    if(valfrom == "sample"){
      valdat1 <- dat %>% filter({{valvar}} == 1) %>% sample_n(valn/2)
      valdat0 <- dat %>% filter({{valvar}} == 0) %>% sample_n(valn/2)
    }
    if(valfrom == "target"){
      tardat <- gen2(i, n, .7, .7, .9, .9, pz = pztarget)
      valdat1 <- tardat %>% filter({{valvar}} == 1) %>% sample_n(valn/2)
      valdat0 <- tardat %>% filter({{valvar}} == 0) %>% sample_n(valn/2)
    }
    valdat <- bind_rows(valdat1, valdat0)
  }
  if(sc == 2){
    dat <- gen2(i, n, .7, .95, .95, .85, pz = 0.75)
    if(valfrom == "sample"){
      valdat1 <- dat %>% filter({{valvar}} == 1) %>% sample_n(valn/2)
      valdat0 <- dat %>% filter({{valvar}} == 0) %>% sample_n(valn/2)
    }
    if(valfrom == "target"){
      tardat <- gen2(i, n, .7, .95, .95, .85, pz = pztarget)
      valdat1 <- tardat %>% filter({{valvar}} == 1) %>% sample_n(valn/2)
      valdat0 <- tardat %>% filter({{valvar}} == 0) %>% sample_n(valn/2)
    }
    valdat <- bind_rows(valdat1, valdat0)
  }
  if(sc == 3){
    dat <- gen3(i, n, .7, .7, .9, .9, pz = 0.75)
    if(valfrom == "sample"){
      valdat1 <- dat %>% filter({{valvar}} == 1) %>% sample_n(valn/2)#ning: revise the number in sample_n() from 100 to valn/2
      valdat0 <- dat %>% filter({{valvar}} == 0) %>% sample_n(valn/2)
    }
    if(valfrom == "target"){
      tardat <- gen3(i, n, .7, .7, .9, .9, pz = pztarget)
      valdat1 <- tardat %>% filter({{valvar}} == 1) %>% sample_n(valn/2)
      valdat0 <- tardat %>% filter({{valvar}} == 0) %>% sample_n(valn/2)
    }
    valdat <- bind_rows(valdat1, valdat0)
  }
  if(sc == 4){#ning: se2 issue
    dat <- gen3(i, n, .7, .9, .95, .85, pz = 0.75)
    if(valfrom == "sample"){
      valdat1 <- dat %>% filter({{valvar}} == 1) %>% sample_n(valn/2)
      valdat0 <- dat %>% filter({{valvar}} == 0) %>% sample_n(valn/2)
    }
    if(valfrom == "target"){
      tardat <- gen3(i, n, .7, .9, .95, .85, pz = pztarget)
      valdat1 <- tardat %>% filter({{valvar}} == 1) %>% sample_n(valn/2)
      valdat0 <- tardat %>% filter({{valvar}} == 0) %>% sample_n(valn/2)
    }
    valdat <- bind_rows(valdat1, valdat0)
  }
  
  mu1 <- mu1(dat, target = tar, pztarget = pztarget)
  mu2 <- mu2_alt(dat, valdat, target = tar, pztarget = pztarget)
  mu3 <- mu3_alt(dat, valdat, target = tar, pztarget = pztarget)
  return(c(mu1, mu2, mu3))
}


##' a function to run second generation of simulations (see `runfig`)
##' 
runfigsc <- function(sc, n, valsize, J, pzvals){
  mu0 <- matrix(NA, nrow = J, ncol = length(pzvals))
  mu1 <- matrix(NA, nrow = J, ncol = length(pzvals))
  mu2 <- matrix(NA, nrow = J, ncol = length(pzvals))
  mu3 <- matrix(NA, nrow = J, ncol = length(pzvals))
  for(i in 1:length(pzvals)){
    for(j in 1:J){
      res <- runfig(j, n, valsize, pzvals[i], tar = "sample", pztarget = NA, sc = sc)
      mu2[j,i] <- res[2]
      mu3[j,i] <- res[3]
      mu1[j,i] <- res[1]
      mu0[j,i] <- res[4]
      print(c(i, j))
    }
  }
  mu0_sum <- colMeans(mu0, na.rm = T)
  mu1_sum <- colMeans(mu1, na.rm = T)
  mu2_sum <- colMeans(mu2, na.rm = T)
  mu3_sum <- colMeans(mu3, na.rm = T)
  mu0_sd <- apply(mu0, 2, FUN = sd, na.rm = T)
  mu1_sd <- apply(mu1, 2, FUN = sd, na.rm = T)
  mu2_sd <- apply(mu2, 2, FUN = sd, na.rm = T)
  mu3_sd <- apply(mu3, 2, FUN = sd, na.rm = T)
  
  resdat <- data.frame(pz = pzvals, mu1 = mu1_sum, mu2 = mu2_sum, mu3 = mu3_sum, 
                       se1 = mu1_sd, se2 = mu2_sd, se3 = mu3_sd,
                       mu0 = mu0_sum, se0 = mu0_sd,
                       target = rep("sample", length(pzvals)))
  
  # target external
  mu0_ext <- matrix(NA, nrow = J, ncol = length(pzvals))
  mu1_ext <- matrix(NA, nrow = J, ncol = length(pzvals))
  mu2_ext <- matrix(NA, nrow = J, ncol = length(pzvals))
  mu3_ext <- matrix(NA, nrow = J, ncol = length(pzvals))
  for(i in 1:length(pzvals)){
    for(j in 1:J){
      res <- runfig(j, n, valsize, pzvals[i], tar = "external", pztarget = .33, sc = sc)
      mu1_ext[j, i] <- res[1]
      mu2_ext[j, i] <- res[2]
      mu3_ext[j, i] <- res[3]
      mu0_ext[j, i] <- res[4]
      print(c(i, j))
    }
  }
  mu0_ext_sum <- colMeans(mu0_ext, na.rm = T)
  mu1_ext_sum <- colMeans(mu1_ext, na.rm = T)
  mu2_ext_sum <- colMeans(mu2_ext, na.rm = T)
  mu3_ext_sum <- colMeans(mu3_ext, na.rm = T)
  
  mu0_ext_sd <- apply(mu0_ext, 2, FUN = sd, na.rm = T)
  mu1_ext_sd <- apply(mu1_ext, 2, FUN = sd, na.rm = T)
  mu2_ext_sd <- apply(mu2_ext, 2, FUN = sd, na.rm = T)
  mu3_ext_sd <- apply(mu3_ext, 2, FUN = sd, na.rm = T)
  
  resdat2 <- data.frame(pz = pzvals, mu1 = mu1_ext_sum, mu2 = mu2_ext_sum, mu3 = mu3_ext_sum, 
                        se1 = mu1_ext_sd, se2 = mu2_ext_sd, se3 = mu3_ext_sd, 
                        mu0 = mu0_ext_sum, se0 = mu0_ext_sd, 
                        target = rep("external", length(pzvals)))
  
  # define truth
  if(sc %in% c(0,1)) dat_true <- gen2(i, 1e6, .7, .7, .9, .9)
  if(sc == 2) dat_true <- gen2(i, 1e6, .7, .95, .95, .85)
  if(sc == 3) dat_true <- gen3(i, 1e6, .7, .7, .9, .9)
  if(sc == 4) dat_true <- gen3(i, 1e6, .7, .9, .95, .85)#ning: se2 issue
  
  true_samp <- mu0(dat_true, target = "sample", pztarget = NA)
  true_ext <- mu0(dat_true, target = "external", pztarget = .33)
  
  allres <- bind_rows(resdat, resdat2) %>% 
    pivot_longer(cols = c(mu0, mu1, mu2, mu3, se0, se1, se2, se3), names_to = c(".value", "estimator"), 
                 names_pattern = "(.*)(.)") %>% 
    mutate(true = ifelse(target == "sample", true_samp, true_ext), 
           bias = mu - true, rmse = sqrt(bias^2 + se^2))
  

  return(allres)
  #return(list(allres, mu3))
}


##' a function to run outcome dependent sims (see `runfig2`)
##' 
runfig_od <- function(sc, n, valsize, J){
  mu1_mat <- matrix(NA, nrow = J, ncol = 4)
  mu2_mat <- matrix(NA, nrow = J, ncol = 4)
  mu3_mat <- matrix(NA, nrow = J, ncol = 4)

  for(j in 1:J){
      res <- runfig2(j, n, valsize, valvar = y, valfrom = "sample", 
                     tar = "sample", pztarget = 0.33, sc = sc)
      mu2_mat[j, 1] <- res[2]
      mu3_mat[j, 1] <- res[3]
      mu1_mat[j, 1] <- res[1]
      
      res <- runfig2(j, n, valsize, valvar = y, valfrom = "target", 
                     tar = "sample", pztarget = 0.33, sc = sc)
      mu2_mat[j, 2] <- res[2]
      mu3_mat[j, 2] <- res[3]
      mu1_mat[j, 2] <- res[1]
      
      res <- runfig2(j, n, valsize, valvar = ystar, valfrom = "sample", tar = "sample", 
                     pztarget = 0.33, sc = sc)
      mu2_mat[j, 3] <- res[2]
      mu3_mat[j, 3] <- res[3]
      mu1_mat[j, 3] <- res[1]
      
      res <- runfig2(j, n, valsize, valvar = ystar, valfrom = "target", tar = "sample", 
                     pztarget = 0.33, sc = sc)
      mu2_mat[j, 4] <- res[2]
      mu3_mat[j, 4] <- res[3]
      mu1_mat[j, 4] <- res[1]
      
      print(c(j))
  }
  
  mu1_sum <- colMeans(mu1_mat, na.rm = T)
  mu2_sum <- colMeans(mu2_mat, na.rm = T)
  mu3_sum <- colMeans(mu3_mat, na.rm = T)
  mu1_sd <- apply(mu1_mat, 2, FUN = sd, na.rm = T)
  mu2_sd <- apply(mu2_mat, 2, FUN = sd, na.rm = T)
  mu3_sd <- apply(mu3_mat, 2, FUN = sd, na.rm = T)
  
  resdat <- data.frame(val=c(1,2,3,4), mu1 = mu1_sum, mu2 = mu2_sum, mu3 = mu3_sum, 
                       se1 = mu1_sd, se2 = mu2_sd, se3 = mu3_sd,
                       target = rep("sample", 4))
  
  # target external
  mu1_ext <- matrix(NA, nrow = J, ncol = 4)
  mu2_ext <- matrix(NA, nrow = J, ncol = 4)
  mu3_ext <- matrix(NA, nrow = J, ncol = 4)
  
    for(j in 1:J){
      res <- runfig2(j, n, valsize, valvar = y, valfrom = "sample", tar = "external", pztarget = .33, sc = sc)
      mu1_ext[j, 1] <- res[1]
      mu2_ext[j, 1] <- res[2]
      mu3_ext[j, 1] <- res[3]
      
      res <- runfig2(j, n, valsize, valvar = y, valfrom = "target", tar = "external", pztarget = .33, sc = sc)
      mu1_ext[j, 2] <- res[1]
      mu2_ext[j, 2] <- res[2]
      mu3_ext[j, 2] <- res[3]
      
      res <- runfig2(j, n, valsize, valvar = ystar, valfrom = "sample", tar = "external", pztarget = .33, sc = sc)
      mu1_ext[j, 3] <- res[1]
      mu2_ext[j, 3] <- res[2]
      mu3_ext[j, 3] <- res[3]
      
      res <- runfig2(j, n, valsize, valvar = ystar, valfrom = "target", tar = "external", pztarget = .33, sc = sc)
      mu1_ext[j, 4] <- res[1]
      mu2_ext[j, 4] <- res[2]
      mu3_ext[j, 4] <- res[3]
      
      print(c(j))
    }
  
  mu1_ext_sum <- colMeans(mu1_ext, na.rm = T)
  mu2_ext_sum <- colMeans(mu2_ext, na.rm = T)
  mu3_ext_sum <- colMeans(mu3_ext, na.rm = T)
  
  mu1_ext_sd <- apply(mu1_ext, 2, FUN = sd, na.rm = T)
  mu2_ext_sd <- apply(mu2_ext, 2, FUN = sd, na.rm = T)
  mu3_ext_sd <- apply(mu3_ext, 2, FUN = sd, na.rm = T)
  
  resdat2 <- data.frame(val=c(1,2,3,4), mu1 = mu1_ext_sum, mu2 = mu2_ext_sum, mu3 = mu3_ext_sum, 
                        se1 = mu1_ext_sd, se2 = mu2_ext_sd, se3 = mu3_ext_sd, 
                        target = rep("external", 4))
  
  # define truth
  if(sc %in% c(0,1)) dat_true <- gen2(111, 1e6, .7, .7, .9, .9) #111 is arbitrary seed
  if(sc == 2) dat_true <- gen2(111, 1e6, .7, .95, .95, .85)
  if(sc == 3) dat_true <- gen3(111, 1e6, .7, .7, .9, .9)
  if(sc == 4) dat_true <- gen3(111, 1e6, .7, .9, .95, .85)#ning: se2 issue
  
  true_samp <- mu0(dat_true, target = "sample", pztarget = NA)
  true_ext <- mu0(dat_true, target = "external", pztarget = .33)
  
  allres <- bind_rows(resdat, resdat2) %>% 
    pivot_longer(cols = c(mu1, mu2, mu3, se1, se2, se3), names_to = c(".value", "estimator"), 
                 names_pattern = "(.*)(.)") %>% 
    mutate(true = ifelse(target == "sample", true_samp, true_ext), 
           bias = mu - true, rmse = sqrt(bias^2 + se^2))
  
  return(allres)
  
}

##' a function to plot bias and CIs from simulations
##' 
ggthemr("flat")
plotfun <- function(data){
  data$targetf = factor(data$target, levels = c("sample", "external"))
  data <- data %>% 
    filter(estimator > 0)
  pl1 <- ggplot()+
    geom_line( aes(x = pz, y = bias*100, color = estimator, group = estimator), 
               data = data) +
    geom_point(aes(x = pz, y = bias*100, group = estimator, color = estimator, shape = estimator), 
               data = data ) +
    geom_ribbon(aes(x = pz, 
                    ymin = bias*100 - 1.96*se*100, 
                    ymax = bias*100 + 1.96*se*100, 
                    fill = estimator), alpha = 0.2, data = data)+
    scale_color_manual("key", values = c("lightgrey", "steelblue", "black"),
                       labels = c("Naive",
                                  "Assuming nondifferential",
                                  "Assuming differential by Z"),
                       name = "Assumption about\nmisclassification")+
    # scale_color_viridis_d(labels = c("Naive", option = "D",
    #                               "Assuming nondifferential",
    #                               "Assuming differential by Z"),
    #                    name = "Assumption about\nmisclassification")+
    scale_fill_manual(guide = "none", values = c("lightgrey", "steelblue", "grey"))+
    scale_shape_discrete(labels = c("Naive",
                                    "Assuming nondifferential",
                                    "Assuming differential by Z"),
                         name = "Assumption about\nmisclassification")+
    #scale_fill_viridis_d(guide = "none", option = "D")+
    facet_wrap(~targetf, labeller = labeller(targetf = c("external" = "External Target\nP(Z = 1) = 0.33", 
                                                         "sample" = "Main study sample\nP(Z = 1) = 0.75")))+
    xlab("P(Z=1) in validation sample") +
    labs(color  = "Guide name",  shape = "Guide name") +
    scale_x_continuous(breaks = c(0, .2, .4, .6, .8, 1), limits = c(0, 1))+
    coord_cartesian(expand = F, ylim = c(-20, 20))+
    ylab("Bias") + #theme_classic()+
    theme(panel.border = element_blank(), strip.background = element_blank(), 
          panel.spacing = unit(2, "lines"), legend.position = c(0.8, 0.145), 
          legend.title = element_blank(), 
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA))
  return(pl1)
}


##' a function to plot mse for simulations
##' 
msefun <- function(data){
  data$targetf = factor(data$target, levels = c("sample", "external"))
  data <- data %>% 
    mutate(rmse = sqrt((bias*100)^2 + (se*100)^2)) %>% 
    filter(estimator > 0)
  pl1 <- ggplot()+
    geom_line(aes(x = pz, y = rmse, color = estimator), 
              data = data ) +
    geom_point(aes(x = pz, y = rmse, color = estimator, shape = estimator), 
               data = data ) +
    scale_shape_discrete(labels = c("Naive", 
                                    "Assuming nondifferential",
                                    "Assuming differential by Z"),
                         name = "Assumption about\nmisclassification")+
    scale_color_manual(values =  c("lightgrey", "steelblue", "black"), 
                       labels = c("Naive", 
                                  "Assuming nondifferential",
                                  "Assuming differential by Z"),
                       name = "Assumption about\nmisclassification")+
    scale_fill_manual(guide = "none", values = c("lightgrey", "steelblue", "black"))+
    facet_wrap(~targetf, labeller = labeller(targetf = c("external" = "External Target\nP(Z = 1) = 0.33", 
                                                         "sample" = "Main study sample\nP(Z = 1) = 0.75")))+
    xlab("P(Z=1) in validation sample") +
    scale_x_continuous(breaks = c(0, .2, .4, .6, .8, 1), limits = c(0, 1))+
    coord_cartesian(expand = F, ylim=c(0, 15))+
    ylab("RMSE") + #theme_classic()+
    theme(panel.border = element_blank(), strip.background = element_blank(), 
          panel.spacing = unit(2, "lines"), legend.position = c(0.8, 0.9), 
          legend.title = element_blank(), 
          legend.background = element_rect(fill = "transparent", color = NA), 
          legend.key.size = unit(.5, "cm"))
  return(pl1)
}
