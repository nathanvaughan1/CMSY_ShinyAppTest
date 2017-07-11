
#This is the original CMSY code identical to that published on GitHub
SchaeferParallelSearch<-function(ni, nyr,sigR,duncert,ct,int.yr,intbio, startbt, ki,i, ri,int.yr.i, nstartbt, yr, end.yr, endbio, npoints, pt){
  ptm<-proc.time()
  # create vectors for viable r, k and bt
  inmemorytable <- vector()
  # parallelised for the points in the r-k space
  inmemorytable <- foreach (i = 1 : npoints, .combine='rbind', .packages='foreach', .inorder=TRUE) %dopar%{
    nsbt = length(startbt)
    VP   <- FALSE
    for(nj in 1:nsbt) {  
      # create empty vector for annual biomasses
      bt <- vector()
      j<-startbt[nj]
      # set initial biomass, including 0.1 process error to stay within bounds
      bt[1]=j*ki[i]*exp(rnorm(1,0, 0.1*sigR))  ## set biomass in first year
      # repeat test of r-k-startbt combination to allow for different random error
      for(re in 1:ni)   {
        #loop through years in catch time series
        for (t in 1:nyr)  {  # for all years in the time series
          xt=rnorm(1,0, sigR) # set new process error for every year  
          zlog.sd = sqrt(log(1+(duncert)^2))
          zt=rlnorm(1,meanlog = 0, sdlog = zlog.sd) # model the catch error as a log normal distribution.
          # calculate biomass as function of previous year's biomass plus surplus production minus catch
          bt[t+1] <- ifelse(bt[t]/ki[i] >= 0.25,
                            bt[t]+ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt,
                            bt[t]+(4*bt[t]/ki[i])*ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt) # assuming reduced r at B/k < 0.25
          
          # if biomass < 0.01 k, discard r-k-startbt combination
          if(bt[t+1] < 0.01*ki[i]) { 
            break
          } # stop looping through years, go to next upper level
          # intermediate year check
          if ((t+1)==int.yr.i && (bt[t+1]>(intbio[2]*ki[i]) || bt[t+1]<(intbio[1]*ki[i]))) { 
            break 
          }  
        } # end of loop of years
        # if loop was broken or last biomass falls outside of expected ranges 
        # do not store results, go directly to next startbt
        if(t < nyr || bt[yr==end.yr] > (endbio[2]*ki[i]) || bt[yr==end.yr] < (endbio[1]*ki[i]) ) { 
          next 
        } else {
          #each vector will be finally appended to the others found by the threads - this is done by the .combine='rbind' option
          inmemorytablerow<-c(i,j,ri[i],ki[i],bt[1:(nyr+1)]/ki[i])
          if (length(inmemorytablerow)==(4+nyr+1)){
            if (VP==FALSE)
            {
              inmemorytable<-inmemorytablerow
            }
            else
            {
              inmemorytable<-rbind(inmemorytable,inmemorytablerow)
            }
            VP<-TRUE
          }
        }
      } # end of repetition for random error
    } # end of j-loop of initial biomasses 
    # instruction necessary to make the foreach loop see the variable:
    if (length(inmemorytable)==0)
    {inmemorytable<-vector(length=4+nyr+1)*NA}
    else
    {inmemorytable}
  }#end loop on points
  
  #create the output matrix
  mdat        <- matrix(data=NA, nrow = npoints*nstartbt, ncol = 2+nyr+1) 
  npointsmem = dim(inmemorytable)[1]
  npointscols = dim(inmemorytable)[2]
  #reconstruction of the processing matrix after the parallel search
  if (npointsmem>0 && npointscols>0){
    for (idxr in 1:npointsmem){
      i = inmemorytable[idxr,1]
      if (!is.na(i)){
        j = inmemorytable[idxr,2]
        mdatindex<-((i-1)*nstartbt)+which(startbt==j)  
        mdat[mdatindex,1]           <- inmemorytable[idxr,3]
        mdat[mdatindex,2]           <- inmemorytable[idxr,4]
        mdat[mdatindex,3:(2+nyr+1)] <- inmemorytable[idxr,5:(4+nyr+1)]
        if(pt==T) points(x=ri[i], y=ki[i], pch=".", cex=4, col="gray")
      }
    }
  }
  ptm<-proc.time()-ptm
  mdat <- na.omit(mdat)
  return(mdat)
}

SchaeferMC <- function(ri, ki, startbio, int.yr, intbio, endbio, sigR, pt, duncert, startbins, ni) {
  
  # create vector for initial biomasses
  startbt     <- seq(from =startbio[1], to=startbio[2], by = (startbio[2]-startbio[1])/startbins)
  nstartbt    <- length(startbt)
  npoints     <- length(ri)
  # get index of intermediate year
  int.yr.i     <- which(yr==int.yr) 
  
  #loop through r-k pairs with parallel search
  mdat<-SchaeferParallelSearch(ni, nyr,sigR,duncert,ct,int.yr,intbio, startbt, ki, i, ri, int.yr.i, nstartbt, yr, end.yr, endbio, npoints,pt)
  
  cat("\n")
  return(list(mdat))
} # end of SchaeferMC function

#This code follows the same logic as the original CMSY but is vectorized and chuncks the data before sending to foreach
SchaeferParallelSearch2<-function(ri, ki, ct, startbt, intbio, endbio, int.yr.i, ni, sigR, duncert, pt){
  
  nyr<-length(ct)        
  nsbt = length(startbt)
  inmemorytable <- matrix(NA,nrow=0,ncol=(nyr+3))
  
  #set uncertianty in catch 
  zlog.sd = sqrt(log(1+(duncert)^2))
  #seperate r and k values into one chunk for each parrallel core
  rb1<-matrix(ri,nrow=ncores_for_computation)
  kb1<-matrix(ki,nrow=ncores_for_computation)
  #print points if graphing
  if(pt){ points(x=ri, y=ki, pch=".", cex=4, col="gray")}
  
  inmemorytable <- foreach (i = 1 : ncores_for_computation, .combine='rbind', .packages='foreach', .inorder=TRUE) %dopar%{
    
    biomass<-matrix(nrow=(length(rb1[i,])*nsbt*ni),ncol=(3*nyr+3))
    
    #Fill columns with r, K, and (B[1]/K)
    biomass[,1]<-c(outer(c(outer(rep(1,ni),rep(1,length(startbt)),"*")),rb1[i,],"*"))
    biomass[,2]<-c(outer(c(outer(rep(1,ni),rep(1,length(startbt)),"*")),kb1[i,],"*"))
    biomass[,(2*nyr+3)]<-c(outer(c(outer(exp(rep(rnorm(1,0,0.1*sigR),ni)),startbt,"*")),rep(1,length(kb1[i,])),"*"))
    
    #Loop over years to calculate biomass series
    for (t in 1:nyr)  { 
      biomass[,(t+2)]=rnorm(length(biomass[,1]),0, sigR) # set new process error for every year  
      biomass[,(t+nyr+2)]=rlnorm(length(biomass[,1]),meanlog = 0, sdlog = zlog.sd) # model the catch error as a log normal distribution.
      biomass[,(t+2*nyr+3)]<-ifelse(biomass[,(t+2*nyr+2)] >= 0.25,
                                    biomass[,(t+2*nyr+2)]+biomass[,1]*biomass[,(t+2*nyr+2)]*(1-biomass[,(t+2*nyr+2)])*exp(biomass[,(t+2)])-(ct[t]/biomass[,2])*biomass[,(t+nyr+2)],
                                    biomass[,(t+2*nyr+2)]+(4*biomass[,(t+2*nyr+2)])*biomass[,1]*biomass[,(t+2*nyr+2)]*(1-biomass[,(t+2*nyr+2)])*exp(biomass[,(t+2)])-(ct[t]/biomass[,2])*biomass[,(t+nyr+2)]) # assuming reduced r at B/k < 0.25
      
      #remove biomass series where stock collapses
      biomass<-biomass[(!is.na(biomass[,(t+2*nyr+3)])),,drop=FALSE]
      biomass<-biomass[biomass[,(t+2*nyr+3)]>=0.01,,drop=FALSE]
      #remove biomass series where intermittent biomass is outside bounds
      if((t+1)==int.yr.i){
        biomass<-biomass[biomass[,(t+2*nyr+3)]>=intbio[1],,drop=FALSE]
        biomass<-biomass[biomass[,(t+2*nyr+3)]<=intbio[2],,drop=FALSE]
      }
    }
    #remove biomass series where final biomass is outside bounds
    biomass<-biomass[biomass[,(3*nyr+2)]>=endbio[1],,drop=FALSE]
    biomass<-biomass[biomass[,(3*nyr+2)]<=endbio[2],,drop=FALSE]
    #remove duplicate r,K,B[1] combinations (only retain one catch jitter)
    if(length(biomass[,1])>1){biomass<-biomass[!duplicated(biomass[,c(1:2,(2*nyr+3))]),,drop=FALSE]}
    # instruction necessary to make the foreach loop see the variable:
    inmemorytable<-rbind(inmemorytable,biomass[,c(1,2,(2*nyr+3):(3*nyr+3))])
  }#end parallelization
  
  return(inmemorytable)
}

SchaeferMC2 <- function(ri, ki, startbio, int.yr, intbio, endbio, sigR, pt, duncert, startbins, ni) {
  
  # create vector for initial biomasses
  startbt     <- seq(from =startbio[1], to=startbio[2], by = (startbio[2]-startbio[1])/startbins)
  # get index of intermediate year
  int.yr.i     <- which(yr==int.yr) 
  #loop through r-k pairs with parallel search
  mdat<-SchaeferParallelSearch2(ri, ki, ct, startbt, intbio, endbio, int.yr.i, ni, sigR, duncert, pt)
  
  cat("\n")
  return(list(mdat))
} # end of SchaeferMC function
