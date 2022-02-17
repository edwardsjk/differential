### utilities for differential misclassification paper EXAMPLE ###


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

mu2 <- function(dat, valdata, target = "sample", pztarget = NA){
  samp <- valdata
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

mu3 <- function(dat, valdata, target = "sample", pztarget = NA){
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

ci_naive <- function(p, n){
  var <- p*(1-p)/n
  se <- sqrt(var)
  lcl = p - 1.96*se
  ucl <- p + 1.96 * se
  ci <- paste0(round(lcl*100, 1), ", ", round(ucl*100, 1))
  return(ci)
}


boot <- function(B, data, valdat, estimator, ...){
  if(B == 0) return(0)
  if(B>0){
    risks <- matrix(0, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(data), B)
    vdatbi<-samp.bootstrap(nrow(valdat), B)
    #vdatbi_y0<-samp.bootstrap(nrow(valdat[y==0,]), B)
    for(i in 1:B){
      dati <- data[datbi[,i],]
      vdati <- valdat[vdatbi[,i],]
      risks[i,] <- estimator(dati, vdati, ...)
    }
  }
  return(risks)
}
