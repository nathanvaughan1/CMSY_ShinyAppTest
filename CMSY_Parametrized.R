


runCMSY <- function(catch_file_ = NULL, id_file_ = null, stocks_ = NA, dataUncert_ = 0.1, sigmaR_ = 0.1, n_ = 1000, ni_ = 3, nab_ = 5) {
        catch_file <<- catch_file_
        id_file <<- id_file_
        stocks <<- stocks_
        dataUncert <<- dataUncert_
        sigmaR <<- sigmaR_
        n <<- n_
        ni <<- ni_
        nab <<- nab_
  
        cat("Running with parameters:\n")
        cat("Catch File: ", catch_file, "\n")
        cat("ID File: ", id_file, "\n")
        cat("Stocks: ", stocks, "\n")
        cat("Data uncertainty: ", dataUncert, "\n")
        cat("Overall process error for CMSY: ", sigmaR, "\n")
        cat("Initial number of r-k pairs: ", n, "\n")
        cat("Iterations for r-k-startbiomass combinations: ", ni, "\n")
        cat("Minimum number of years with abundance data to run BSM: ", nab, "\n\n\n")
  
        library("parallel")
        library("foreach")
        library("doParallel")
        library("gplots")
        #-----------------------------------------
        # Some general settings
        #-----------------------------------------
        # set.seed(999) # use for comparing results between runs
        #rm(list=ls(all=TRUE)) # clear previous variables etc
        options(digits=3) # displays all numbers with three significant digits as default
        graphics.off() # close graphics windows from previous sessions
        FullSchaefer <- F    # initialize variable; automatically set to TRUE if enough abundance data are available
        ncores_for_computation<<-(detectCores()) # cores to be used for parallel processing of CMSY
        cl           <- makeCluster(ncores_for_computation)
        registerDoParallel(cl, cores = ncores_for_computation)
        
        #-----------------------------------------
        # Required settings, File names
        #-----------------------------------------
        #catch_file<-"/home/enrico/Work/BlueBridge/CMSY/CMSY_with_UserGuide_and_examples/O_Stocks_Catch_14_Med.csv"
        #id_file<-"/home/enrico/Work/BlueBridge/CMSY/CMSY_with_UserGuide_and_examples/O_Stocks_ID_17_Med.csv"
        
        #----------------------------------------
        # Select stock to be analyzed
        #----------------------------------------
        #stocks      <-NA
        #stocks <- "MERLMER_SA"
        # If the input files contain more than one stock, specify below the stock to be analyzed
        # If the line below is commented out (#), all stocks in the input file will be analyzed
        #stocks <-  "SARDPIL_SA_EXAMPLE" # c("SEPIOFF_CY","MICRPOU_IS","EPINGUA_IS","CHAMGAL_SA","CORYHIP_SA","ILLECOI_SA")
        
        #-----------------------------------------
        # General settings for the analysis
        #-----------------------------------------
        #dataUncert   <- 0.1  # set observation error as uncertainty in catch - default is SD=0.1
        #sigmaR       <- 0.1 # overall process error for CMSY; SD=0.1 is the default
        #n            <- 10000 # initial number of r-k pairs
        #n            <- 1000 # initial number of r-k pairs
        n.new        <- n # initialize n.new
        #ni           <- 3 # iterations for r-k-startbiomass combinations, to test different variability patterns; no improvement seen above 3
        #nab          <- 5 # default=5; minimum number of years with abundance data to run BSM
        mgraphs      <- T # set to TRUE to produce additional graphs for management
        save.plots   <- T # set to TRUE to save graphs to JPEG files
        close.plots  <- T # set to TRUE to close on-screen plots after they are saved, to avoid "too many open devices" error in batch-processing
        write.output <- T # set to TRUE if table with results in output file is wanted; expects years 2004-2010 to be available
        force.cmsy   <- T # set to TRUE if CMSY results are to be preferred over BSM results
        select.yr    <- NA # option to display F, B, F/Fmsy and B/Bmsy for a certain year; default NA
        
        #----------------------------------------------
        #  FUNCTIONS
        #----------------------------------------------
        # Schaefer biomass calculation functions
        #----------------------------------------------
        
        source(paste0(getwd(),"/SchaeferMCfuncs.R"))
        
        #----------------------------------------------
        #r-K-b[1] parameter space sampling strategie functions
        #----------------------------------------------
        
        source(paste0(getwd(),"/AltCMSYs.R"))
        
        #-----------------------------------------------
        # Function for moving average
        #-----------------------------------------------
        ma    <- function(x){
          cat(x, "\n")
          cat(rep(1/3,3), "\n")
          cat(typeof(x), "\n")
          cat(typeof(rep(1/3,3)), "\n")
          x.1    <-   stats::filter(x,rep(1/3,3),sides=1)
          x.1[1] <- x[1]
          x.1[2] <- (x[1]+x[2])/2
          return(x.1)
        }
        
        #-----------------------------------------------
        #Function for truncated normal random sample
        #-----------------------------------------------
        rtnorm<-function(n,mean=0,sd=1,min=-3,max=3)
        {
          vals<-mean+(qnorm(runif(n,pnorm((min-mean)/sd),pnorm((max-mean)/sd))))*sd
        }
        
        #---------------------------------------------
        # END OF FUNCTIONS
        #---------------------------------------------
        
        #------------------------------------------
        # Read data and assign to vectors
        #------------------------------------------
        
        #catch_file<-"Example_Catch1.csv"
        #id_file<-"Example_ID1.csv"
        # Read data
        cdat         <- read.csv(catch_file, header=T, dec=".", stringsAsFactors = FALSE)
        cinfo        <- read.csv(id_file, header=T, dec=".", stringsAsFactors = FALSE)
        
        #---------------------------------
        # Analyze stock(s)
        #---------------------------------
        if(is.na(stocks[1])==TRUE){
          stocks         <- as.character(cinfo$Stock) # Analyze stocks in sequence of ID file
          # stocks         <- sort(as.character(cinfo$Stock)) # Analyze stocks in alphabetic order
          # stocks         <- as.character(cinfo$Stock[cinfo$Subregion=="Sardinia"]) # Analyze stocks in Region
        }
        
        resAll<-list()
        resAll$Orig<-list()
        resAll$Vect<-list()
        timeResults<-matrix(NA,nrow=(length(stocks)),ncol=3)
        countResults<-matrix(NA,nrow=(length(stocks)),ncol=3)
        # analyze one stock after the other
        i=0
        for(stock in stocks) {
          i=i+1
          cat("Starting species",i," with name ", stock," \n")
          # assign data from cinfo to vectors
          res          <<- as.character(cinfo$Resilience[cinfo$Stock==stock])
          start.yr     <<- as.numeric(cinfo$StartYear[cinfo$Stock==stock])
          end.yr       <<- as.numeric(cinfo$EndYear[cinfo$Stock==stock])
          r.low        <<- as.numeric(cinfo$r.low[cinfo$Stock==stock])
          r.hi         <<- as.numeric(cinfo$r.hi[cinfo$Stock==stock])
          user.log.r   <<- ifelse(is.na(r.low)==F & is.na(r.hi)==F,TRUE,FALSE)     
          stb.low      <<- as.numeric(cinfo$stb.low[cinfo$Stock==stock])
          stb.hi       <<- as.numeric(cinfo$stb.hi[cinfo$Stock==stock])
          int.yr       <<- as.numeric(cinfo$int.yr[cinfo$Stock==stock])
          intb.low     <<- as.numeric(cinfo$intb.low[cinfo$Stock==stock])
          intb.hi      <<- as.numeric(cinfo$intb.hi[cinfo$Stock==stock])
          endb.low     <<- as.numeric(cinfo$endb.low[cinfo$Stock==stock])
          endb.hi      <<- as.numeric(cinfo$endb.hi[cinfo$Stock==stock])
          btype        <<- as.character(cinfo$btype[cinfo$Stock==stock])
          force.cmsy   <<- ifelse(force.cmsy==T,T,cinfo$force.cmsy[cinfo$Stock==stock])
          comment      <- as.character(cinfo$Comment[cinfo$Stock==stock])
          # set global defaults for uncertainty
          duncert      <<- dataUncert
          sigR         <<- sigmaR
          
          # check for common errors
          if (length(btype)==0){
            cat("ERROR: Could not find the stock in the ID input file - check that the stock names match in ID and Catch files and that commas are used (not semi-colon)")
            return (NA) }
          if(start.yr < cdat$yr[cdat$Stock==stock][1]){
            cat("ERROR: start year in ID file before first year in catch file\n")
            return (NA)}
          
          # extract data on stock
          yr           <<- as.numeric(cdat$yr[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])
          
          if (length(yr)==0){
            cat("ERROR: Could not find the stock in the Catch input files - Please check that the code is written correctly")
            return (NA)
          }
          
          ct.raw           <- as.numeric(cdat$ct[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])/1000  ## assumes that catch is given in tonnes, transforms to '000 tonnes
          if(btype=="biomass" | btype=="CPUE" ) {
            bt <<- as.numeric(cdat$bt[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])/1000  ## assumes that biomass is in tonnes, transforms to '000 tonnes
          } else {bt <- NA}
          
          if(is.na(mean(ct.raw))){
            cat("ERROR: Missing value in Catch data; fill or interpolate\n")  
          }
          nyr          <<- length(yr) # number of years in the time series
          
          # change catch to 3 years moving average where value is average of past 3 years 
          ct              <<- ma(ct.raw)
          
          # initialize vectors for viable r, k, bt, and all in a matrix
          mdat.all    <<- matrix(data=vector(),ncol=2+nyr+1)
          
          # initialize other vectors anew for each stock
          current.attempts <- NA
          
          # use start.yr if larger than select year
          if(is.na(select.yr)==F) {
            sel.yr <<- ifelse(start.yr > select.yr,start.yr,select.yr)
          } else sel.yr <<- NA
          
          #----------------------------------------------------
          # Determine initial ranges for parameters and biomass
          #----------------------------------------------------
          # initial range of r from input file
          if(is.na(r.low)==F & is.na(r.hi)==F) {
            start.r <<- c(r.low,r.hi)
          } else {
            # initial range of r based on resilience
            if(res == "High") {
              start.r <<- c(0.6,1.5)} else if(res == "Medium") {
                start.r <<- c(0.2,0.8)}    else if(res == "Low") {
                  start.r <<- c(0.05,0.5)}  else { # i.e. res== "Very low"
                    start.r <<- c(0.015,0.1)} 
          }
          
          # get index of years with lowest and highest catch between start+3 and end-3 years
          min.yr.i     <<- which.min(ct[4:(length(ct)-3)])+3
          max.yr.i     <<- which.max(ct[4:(length(ct)-3)])+3
          min.ct       <<- ct[min.yr.i]
          max.ct       <<- ct[max.yr.i]
          
          # use initial biomass range from input file if stated
          if(is.na(stb.low)==F & is.na(stb.hi)==F) {
            startbio <<- c(stb.low,stb.hi)
          } else {
            # if start year < 1960 assume high biomass
            if(start.yr < 1960) {startbio <- c(0.5,0.9)} else {
              # else use medium prior biomass range
              startbio <<- c(0.2,0.6)} }
          
          # use year and biomass range for intermediate biomass from input file
          if(is.na(intb.low)==F & is.na(intb.hi)==F) {
            int.yr   <<- int.yr
            intbio   <<- c(intb.low,intb.hi)
            
            # if contrast in catch is low, use initial range again in mid-year
          } else if(min(ct)/max(ct) > 0.6) {
            int.yr    <<- as.integer(mean(c(start.yr, end.yr)))
            intbio    <<- startbio 
            
            # else if year of minimum catch is after max catch then use min catch
          } else if(min.yr.i > max.yr.i) {
            int.yr    <- yr[min.yr.i-1]
            if(startbio[1]>=0.5 &  (int.yr-start.yr) < (end.yr-int.yr) & 
               (min.ct/max.ct) > 0.3) intbio <<- c(0.2,0.6) else intbio <- c(0.01,0.4)
               
               # else use max catch  
          } else {
            # assume that biomass range in year before maximum catch was high or medium
            int.yr    <<- yr[max.yr.i-1]
            intbio    <<- if((startbio[1]>=0.5 & (int.yr-start.yr) < (end.yr-int.yr))| # if initial biomass is high, assume same for intermediate
                            # ((min.ct/max.ct < 0.3 & (max.yr.i - min.yr.i) < 25))) c(0.5,0.9) else c(0.2,0.6) }
                            (((max.ct-min.ct)/max.ct)/(max.yr.i-min.yr.i) > 0.04)) c(0.5,0.9) else c(0.2,0.6) } # if incease is steep, assume high, else medium
          # end of intbio setting
          
          # final biomass range from input file
          if(is.na(endb.low)==F & is.na(endb.hi)==F) {
            endbio   <<- c(endb.low,endb.hi)
          } else {
            # else use mean final catch/max catch to estimate final biomass
            rawct.ratio=ct.raw[nyr]/max(ct)
            endbio  <<- if(ct[nyr]/max(ct) > 0.8) {c(0.4,0.8)} else if(rawct.ratio < 0.5) {c(0.01,0.4)} else {c(0.2,0.6)}
            
            # if default endbio is low (0.01-0.4), check whether the upper bound should be lower than 0.4 for depleted stocks
            if(endbio[2]==0.4){
              if(rawct.ratio< 0.05) {endbio[2] <<- 0.1} else
                if(rawct.ratio< 0.15) {endbio[2] <<- 0.2} else
                  if(rawct.ratio< 0.35) {endbio[2] <<- 0.3} else {endbio[2] <<- 0.4}
            }
          } # end of final biomass setting
          
          # initial prior range of k values, assuming min k will be larger than max catch / prior for r 
          if(mean(endbio) <= 0.5) {
            start.k <<- c(max(ct)/start.r[2],4*max(ct)/start.r[1])} else {
              start.k <<- c(2*max(ct)/start.r[2],12*max(ct)/start.r[1])} 
          cat("startbio=",startbio,ifelse(is.na(stb.low)==T,"default","expert"),
              ", intbio=",int.yr,intbio,ifelse(is.na(intb.low)==T,"default","expert"),
              ", endbio=",endbio,ifelse(is.na(endb.low)==T,"default","expert"),"\n")
          
          
          #Sample acceptable r-K-b combos useing original CMSY code
          stOrig<-proc.time()[3]
          resOrig<-runOrig()
          endOrig<-proc.time()[3]
          tOrig<-endOrig-stOrig
          
          #Sample acceptable r-K-b combos useing vectorized CMSY code  
          stVect<-proc.time()[3]
          resVect<-runVect()
          endVect<-proc.time()[3]
          tVect<-endVect-stVect
          
          #Store run time results for each species in a matrix
          timeResults[i,1]<-tOrig
          timeResults[i,2]<-tVect
          
          #Store number of accepted results for each species in a matrix
          countResults[i,1]<-length(resOrig[,1])
          countResults[i,2]<-length(resVect[,1])
          
          #Store full results for each species in a list (These could be fead back in to the getResults function to calculate MSY etc)
          resAll$Orig[[i]]<-resOrig
          resAll$Vect[[i]]<-resVect
          resAll$comments[[i]] <- "Some comments can be added here..."
          #Not currently running all the calculations after getting accepted r-K-b combos as this code is just to test speed.
          
        } # end of stocks loop
        
        #stop parallel processing clusters
        stopCluster(cl)
        stopImplicitCluster()
        cat("Returning object....\n")
        return (resAll);
}
