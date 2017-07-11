#This runs the original CMSY as published
runOrig<-function()
{
  #------------------------------------------------------------------
  # Uniform sampling of the r-k space
  #------------------------------------------------------------------
  # get random set of r and k from log space distribution 
  ri = exp(runif(n, log(start.r[1]), log(start.r[2])))  
  ki = exp(runif(n, log(start.k[1]), log(start.k[2])))  
  
  #---------------------------------------------------------------------
  # 1 - Call CMSY-SchaeferMC function to preliminary explore the r-k space
  #---------------------------------------------------------------------
  cat("First Monte Carlo filtering of r-k space with ",n," points...\n")
  MCA <-  SchaeferMC(ri=ri, ki=ki, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                     pt=F, duncert=duncert, startbins=10, ni=ni)
  
  mdat.all <- MCA[[1]]
  rv.all   <- mdat.all[,1]
  kv.all   <- mdat.all[,2]
  btv.all  <- mdat.all[,3:(2+nyr+1)]
  btc.all  <- mdat.all[,(nyr+2)]-mdat.all[,3] 
  # count viable trajectories and r-k pairs 
  n.viable.b   <- length(mdat.all[,1])
  n.viable.pt <- length(unique(mdat.all[,1]))
  cat("Found ",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
  
  #----------------------------------------------------------------------- 
  # 2 - if the lower bound of k is too high, reduce it by half and rerun
  #-----------------------------------------------------------------------
  if(length(kv.all[kv.all < 1.1*start.k[1] & rv.all < mean(start.r)]) > 10) {
    cat("Reducing lower bound of k, resampling area with",n,"additional points...\n")
    start.k <- c(0.5*start.k[1],start.k[2])
    ri1 = exp(runif(n, log(start.r[1]), log(start.r[2])))  
    ki1 = exp(runif(n, log(start.k[1]), log(start.k[2])))  
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                       pt=F, duncert=dataUncert, startbins=10, ni=ni)
    mdat.all <- rbind(mdat.all,MCA[[1]])
    rv.all   <- mdat.all[,1]
    kv.all   <- mdat.all[,2]
    btv.all  <- mdat.all[,3:(2+nyr+1)]
    n.viable.b   <- length(mdat.all[,1])
    n.viable.pt <- length(unique(mdat.all[,1]))
    cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
  }
  
  #-------------------------------------------------------------------
  # 3 - if few points were found then resample and shrink the log k space
  #-------------------------------------------------------------------
  if (n.viable.b <= 1000){
    log.start.k.new  <- log(start.k) 
    max.attempts     <- 3
    current.attempts <- 1
    startbins        <- 10  
    while (n.viable.b <= 1000 && current.attempts <= max.attempts){
      if(n.viable.pt > 0) {
        log.start.k.new[1] <- mean(c(log(start.k[1]), min(log(kv.all))))
        log.start.k.new[2] <- mean(c(log.start.k.new[2], max(log(kv.all)))) }
      n.new <- n*current.attempts #add more points
      ri1 = exp(runif(n.new, log(start.r[1]), log(start.r[2])))  
      ki1 = exp(runif(n.new, log.start.k.new[1], log.start.k.new[2]))
      cat("Shrinking k space: repeating Monte Carlo in the interval [",exp(log.start.k.new[1]),",",exp(log.start.k.new[2]),"]\n")
      cat("Attempt ",current.attempts," of ",max.attempts," with ",n.new," additional points...","\n")
      if(current.attempts==2 & n.viable.b < 50){
        duncert   <- 2*dataUncert
        sigR      <- 2*sigmaR
        startbins <- 20
        cat("Doubling startbins, catch and process error, and number of variability patterns \n")   
      }
      MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                         pt=F, duncert=duncert, startbins=startbins, ni=2*ni)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
      current.attempts=current.attempts+1 #increment the number of attempts
    }
    if(n.viable.b < 5) {
      cat("Only",n.viable.pt,"viable r-k pairs found, check data and settings \n")
    }  
  }
  if(length(rv.all)>0)
  {
    #------------------------------------------------------------------
    # 4 - if tip of viable r-k pairs is 'thin', do extra sampling there
    #------------------------------------------------------------------
    if(length(rv.all[rv.all > 0.9*start.r[2]]) < 5) { 
      l.sample.r        <- quantile(rv.all,0.6)
      add.points        <- ifelse(is.na(current.attempts)==T,n,ifelse(current.attempts==2,2*n,ifelse(length(rv.all)>500,3*n,6*n)))
      cat("Final sampling in the tip area above r =",l.sample.r,"with",add.points,"additional points...\n")
      log.start.k.new <- c(log(0.8*min(kv.all)),log(max(kv.all[rv.all > l.sample.r])))
      
      ri1 = exp(runif(add.points, log(l.sample.r), log(start.r[2])))  
      ki1 = exp(runif(add.points, log.start.k.new[1], log.start.k.new[2]))
      MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                         pt=F, duncert=duncert, startbins=10, ni=ni)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
    }
  }
  return(mdat.all)
}

#This code is identical to runOrig except it calls ShaeferMC2 which is a vectorized version of ShaeferMC
runVect<-function()
{
  #------------------------------------------------------------------
  # Uniform sampling of the r-k space
  #------------------------------------------------------------------
  # get random set of r and k from log space distribution 
  ri1 = exp(runif(n, log(start.r[1]), log(start.r[2])))  
  ki1 = exp(runif(n, log(start.k[1]), log(start.k[2])))  
  
  #---------------------------------------------------------------------
  # 1 - Call CMSY-SchaeferMC function to preliminary explore the r-k space
  #---------------------------------------------------------------------
  cat("First Monte Carlo filtering of r-k space with ",n," points...\n")
  MCA <-  SchaeferMC2(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                     pt=F, duncert=dataUncert, startbins=10, ni=ni)
  
  mdat.all <- MCA[[1]]
  rv.all   <- mdat.all[,1]
  kv.all   <- mdat.all[,2]
  btv.all  <- mdat.all[,3:(2+nyr+1)]
  btc.all  <- mdat.all[,(nyr+2)]-mdat.all[,3] 
  # count viable trajectories and r-k pairs 
  n.viable.b   <- length(mdat.all[,1])
  n.viable.pt <- length(unique(mdat.all[,1]))
  cat("Found ",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
  
  #----------------------------------------------------------------------- 
  # 2 - if the lower bound of k is too high, reduce it by half and rerun
  #-----------------------------------------------------------------------
  if(length(kv.all[kv.all < 1.1*start.k[1] & rv.all < mean(start.r)]) > 10) {
    cat("Reducing lower bound of k, resampling area with",n,"additional points...\n")
    start.k <- c(0.5*start.k[1],start.k[2])
    ri1 = exp(runif(n, log(start.r[1]), log(start.r[2])))  
    ki1 = exp(runif(n, log(start.k[1]), log(start.k[2])))  
    MCA <-  SchaeferMC2(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                       pt=F, duncert=dataUncert, startbins=10, ni=ni)
    mdat.all <- rbind(mdat.all,MCA[[1]])
    rv.all   <- mdat.all[,1]
    kv.all   <- mdat.all[,2]
    btv.all  <- mdat.all[,3:(2+nyr+1)]
    n.viable.b   <- length(mdat.all[,1])
    n.viable.pt <- length(unique(mdat.all[,1]))
    cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
  }
  
  #-------------------------------------------------------------------
  # 3 - if few points were found then resample and shrink the log k space
  #-------------------------------------------------------------------
  if (n.viable.b <= 1000){
    log.start.k.new  <- log(start.k) 
    max.attempts     <- 3
    current.attempts <- 1
    startbins        <- 10  
    while (n.viable.b <= 1000 && current.attempts <= max.attempts){
      if(n.viable.pt > 0) {
        log.start.k.new[1] <- mean(c(log(start.k[1]), min(log(kv.all))))
        log.start.k.new[2] <- mean(c(log.start.k.new[2], max(log(kv.all)))) }
      n.new <- n*current.attempts #add more points
      ri1 = exp(runif(n.new, log(start.r[1]), log(start.r[2])))  
      ki1 = exp(runif(n.new, log.start.k.new[1], log.start.k.new[2]))
      cat("Shrinking k space: repeating Monte Carlo in the interval [",exp(log.start.k.new[1]),",",exp(log.start.k.new[2]),"]\n")
      cat("Attempt ",current.attempts," of ",max.attempts," with ",n.new," additional points...","\n")
      if(current.attempts==2 & n.viable.b < 50){
        duncert   <- 2*dataUncert
        sigR      <- 2*sigmaR
        startbins <- 20
        cat("Doubling startbins, catch and process error, and number of variability patterns \n")   
      }
      MCA <-  SchaeferMC2(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                         pt=F, duncert=duncert, startbins=startbins, ni=2*ni)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
      current.attempts=current.attempts+1 #increment the number of attempts
    }
    if(n.viable.b < 5) {
      cat("Only",n.viable.pt,"viable r-k pairs found, check data and settings \n")
    }  
  }
  
  if(length(rv.all)>0)
  {
    #------------------------------------------------------------------
    # 4 - if tip of viable r-k pairs is 'thin', do extra sampling there
    #------------------------------------------------------------------
    if(length(rv.all[rv.all > 0.9*start.r[2]]) < 5) { 
      l.sample.r        <- quantile(rv.all,0.6)
      add.points        <- ifelse(is.na(current.attempts)==T,n,ifelse(current.attempts==2,2*n,ifelse(length(rv.all)>500,3*n,6*n)))
      cat("Final sampling in the tip area above r =",l.sample.r,"with",add.points,"additional points...\n")
      log.start.k.new <- c(log(0.8*min(kv.all)),log(max(kv.all[rv.all > l.sample.r])))
      
      ri1 = exp(runif(add.points, log(l.sample.r), log(start.r[2])))  
      ki1 = exp(runif(add.points, log.start.k.new[1], log.start.k.new[2]))
      MCA <-  SchaeferMC2(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                         pt=F, duncert=duncert, startbins=10, ni=ni)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
    }
  }
  return(mdat.all)
}


getResults<-function(mdat.all,start.r,ct.raw)
{
  #------------------------------------
  # get results from CMSY
  #------------------------------------
  # get estimate of most probable r as 75th percentile of mid log.r-classes
  # get unique combinations of r-k
  unique.rk         <- unique(mdat.all[,1:2])
  # get remaining viable log.r and log.k 
  log.rs           <- log(unique.rk[,1])
  log.ks           <- log(unique.rk[,2])
  # get vectors with numbers of r and mid values in classes
  # determine number of classes as a function of r-width
  r.width         <- (max(unique.rk[,1])-start.r[1])/(start.r[2]-start.r[1])
  classes         <- ifelse(r.width>0.8,100,ifelse(r.width>0.5,50,ifelse(r.width>0.3,25,12)))
  hist.log.r      <- hist(x=log.rs, breaks=classes, plot=T)
  log.r.counts    <- hist.log.r$counts
  log.r.mids      <- hist.log.r$mids
  # get most probable log.r as 75th percentile of mids with counts > 0
  log.r.est       <- as.numeric(quantile(log.r.mids[which(log.r.counts > 0)],0.75))
  median.log.r    <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.50))
  lcl.log.r       <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.5125))
  ucl.log.r       <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.9875))
  sd.log.r.est    <- (ucl.log.r - log.r.est) / 1.96 
  r.est           <- exp(log.r.est)
  lcl.r.est       <- exp(log.r.est-1.96*sd.log.r.est)
  ucl.r.est       <- exp(log.r.est+1.96*sd.log.r.est)
  
  # get r-k pairs above median of mids
  rem            <- which(unique.rk[,1] > exp(median.log.r))
  rem.log.r      <- log(unique.rk[,1][rem])
  rem.log.k      <- log(unique.rk[,2][rem])
  # do linear regression of log k ~ log r with slope fixed to -1 (from Schaefer)
  reg            <- lm(rem.log.k ~ 1 + offset(-1*rem.log.r))
  int.reg        <- as.numeric(reg[1])
  sd.reg      <- sd(resid(reg))
  # get estimate of log(k) from y where x = log.r.est
  log.k.est      <- int.reg + (-1) * log.r.est
  # get estimates of ucl of log.k.est from y + SD where x = ucl.log.r  
  ucl.log.k     <- int.reg + (-1) * lcl.log.r + sd.reg
  # get estimates of sd.log.k.est from upper confidence limit of log.k.est  
  sd.log.k.est   <- (ucl.log.k - log.k.est) / 1.96 
  lcl.log.k      <- log.k.est - 1.96*sd.log.k.est
  ucl.log.k      <- log.k.est + 1.96*sd.log.k.est
  
  k.est       <- exp(log.k.est)
  lcl.k.est   <- exp(lcl.log.k)
  ucl.k.est   <- exp(ucl.log.k)
  
  # get MSY from remaining log r-k pairs
  log.MSY.est     <- mean(rem.log.r + rem.log.k - log(4))
  sd.log.MSY.est  <- sd(rem.log.r + rem.log.k - log(4))
  lcl.log.MSY.est <- log.MSY.est - 1.96*sd.log.MSY.est 
  ucl.log.MSY.est <- log.MSY.est + 1.96*sd.log.MSY.est
  MSY.est         <- exp(log.MSY.est)
  lcl.MSY.est     <- exp(lcl.log.MSY.est)
  ucl.MSY.est     <- exp(ucl.log.MSY.est)
  
  # get predicted biomass vectors as median and quantiles 
  # only use biomass trajectories from r-k pairs within the confidence limits
  rem.btv.all <- mdat.all[which(mdat.all[,1] > lcl.r.est & mdat.all[,1] < ucl.r.est 
                                & mdat.all[,2] > lcl.k.est & mdat.all[,2] < ucl.k.est),3:(2+nyr+1)]
  median.btv <- apply(rem.btv.all,2, median)
  median.btv.lastyr  <- median.btv[length(median.btv)-1]
  nextyr.bt  <- median.btv[length(median.btv)]
  lcl.btv    <- apply(rem.btv.all,2, quantile, probs=0.025)
  q.btv      <- apply(rem.btv.all,2, quantile, probs=0.25)
  ucl.btv    <- apply(rem.btv.all,2, quantile, probs=0.975)
  lcl.median.btv.lastyr <- lcl.btv[length(lcl.btv)-1]
  ucl.median.btv.lastyr <- ucl.btv[length(lcl.btv)-1]  
  lcl.nextyr.bt <- lcl.btv[length(lcl.btv)]
  ucl.nextyr.bt <- ucl.btv[length(lcl.btv)]
  
  # get F derived from predicted CMSY biomass
  F.CMSY      <- ct.raw/(median.btv[1:nyr]*k.est)
  Fmsy.CMSY  <- r.est/2 # Fmsy from CMSY 
  
  # --------------------------------------------
  # Get results for management
  # --------------------------------------------
  MSY   <-MSY.est; lcl.MSY<-lcl.MSY.est; ucl.MSY<-ucl.MSY.est 
  Bmsy  <-k.est/2; lcl.Bmsy<-lcl.k.est/2; ucl.Bmsy<-ucl.k.est/2
  Fmsy  <-r.est/2; lcl.Fmsy<-lcl.r.est/2; ucl.Fmsy<-ucl.r.est/2
  B.Bmsy<-2*median.btv[1:nyr];lcl.B.Bmsy<-2*lcl.btv[1:nyr];ucl.B.Bmsy<-2*ucl.btv[1:nyr]
  if(is.na(sel.yr)==F){B.Bmsy.sel<-2*median.btv[yr==sel.yr]}
  
  B          <-B.Bmsy*Bmsy;lcl.B<-lcl.B.Bmsy*Bmsy;ucl.B<-ucl.B.Bmsy*Bmsy
  B.last     <-B[nyr];lcl.B.last<-lcl.B[nyr];ucl.B.last<-ucl.B[nyr]
  B.Bmsy.last<-B.Bmsy[nyr];lcl.B.Bmsy.last<-lcl.B.Bmsy[nyr];ucl.B.Bmsy.last<-ucl.B.Bmsy[nyr]
  
  Fm           <- ct.raw/B;lcl.F<-ct.raw/ucl.B;ucl.F<-ct.raw/lcl.B
  Fmsy.vec     <- ifelse(B.Bmsy>0.5,Fmsy,Fmsy*2*B.Bmsy)
  lcl.Fmsy.vec <- ifelse(B.Bmsy>0.5,lcl.Fmsy,lcl.Fmsy*2*B.Bmsy)
  ucl.Fmsy.vec <- ifelse(B.Bmsy>0.5,ucl.Fmsy,ucl.Fmsy*2*B.Bmsy)
  F.Fmsy       <- Fm/Fmsy.vec; lcl.F.Fmsy<-lcl.F/Fmsy.vec; ucl.F.Fmsy<-ucl.F/Fmsy.vec
  
  F.last     <-Fm[nyr];lcl.F.last<-lcl.F[nyr];ucl.F.last<-ucl.F[nyr]
  Fmsy.last  <-Fmsy.vec[nyr];lcl.Fmsy.last<-lcl.Fmsy.vec[nyr];ucl.Fmsy.last<-ucl.Fmsy.vec[nyr]
  F.Fmsy.last<-F.Fmsy[nyr];lcl.F.Fmsy.last<-lcl.F.Fmsy[nyr];ucl.F.Fmsy.last<-ucl.F.Fmsy[nyr]
  
  if(is.na(sel.yr)==F){
    B.sel<-B.Bmsy.sel*Bmsy
    F.sel<-ct.raw[yr==sel.yr]/B.sel
    F.Fmsy.sel<-F.sel/Fmsy.vec[yr==sel.yr]
  }
  
}
