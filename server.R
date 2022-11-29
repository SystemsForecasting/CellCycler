require(deSolve)
require(ggplot2)
require(reshape2)


# function to find area under curve using Simpson's rule with step = 0.1
simpson <- function(fun, rng, h = 0.1) {
  # rng is time range, n supposed to be an even integer
  a <- rng[1]
  b <- rng[2]
  n <- (b-a)/h
  x <- seq(a, b, by=h)
  if (n == 2) {
    s <- fun(x[1]) + 4*fun(x[2]) +fun(x[3])
  } else {
    s <- fun(x[1]) + fun(x[n+1]) + 2*sum(fun(x[seq(2,n,by=2)])) + 4 *sum(fun(x[seq(3,n-1, by=2)]))
  }
  s <- s*h/3
  return(s)
}

# function to find area under curve using Simpson's rule from time series
simpsonvec <- function(t,y) {
  nt <- length(t)
  n <- 2*floor(nt/2) # round off so even number
  rnd <- n-nt # can be 0 or 0.1
  a <- t[1]
  b <- t[2] - rnd
  h <- 0.1
  t2 <- seq(a, b, h)
  if (n == 2) {
    s <- y[1] + 4*y[2] +y[3]
  } else {
    s <- y[1] + y[n+1] + 2*sum(y[seq(2,n,by=2)]) + 4 *sum(y[seq(3,n-1, by=2)])
  }
  s <- s*h/3
  return(s)
}
# 
# function to define 1-compartment pk model
# input is parameter list pkpar
# creates a list with first element pk as a function of time
createonecmptpkmodel <- function(dose,schedule,ka=0.2,v=25,cl=10) {
  toff <- schedule[1]
  drugfunc <- function(t) {
    t <- t - toff
    t <- pmax(0,t)
    dose*ka/(v*ka-cl)*(exp(-cl/v*t) - exp(-ka*t))
  }
  pkpar <- list(dose=dose,ka=ka,v=v,cl=cl)
  pkmodellist <- list(drugfunc=drugfunc,pkpar=pkpar)
}

# general look-up model, finds ypk corresponding to time tpk for use in ode solver
createlookupmodel <- function(tpk,ypk) {
  drugfunc <- approxfun(c(tpk,1e10),c(ypk,0)) # add point at (infinity,0)
  return(drugfunc)
}

# runs an ode-based pk model, returns a function which looks up value for use in ode
# parms is a list with dose, schedule, ...
# step is -1 for infusion, toggles dose on and off
createlookupfromode <- function(func,tmax,xini,parms,iout=1,step=1) {
  dose <- parms$dose
  schedule <- parms$schedule
  schedule <- schedule[schedule < tmax]  # truncate to max time
  nsch <- length(schedule)
  tsch1 <- schedule[1]
  out <- list(xini)
  nx <- length(xini)
  indy <- seq(2,nx+1) # index for results, time is first column
  if (tsch1 > 0) {  # no drug for initial segment
    tsim1 <- seq(0,tsch1,by=0.1)
    out <- ode(y=xini,times=tsim1,func=func,parms=parms) # run ode solver from t=0
    np <- length(tsim1)
    xini <- out[np,indy]  # update xini for next segment
  }
  for (n in 1:nsch) {
    tschn <- schedule[n]  # start point
    if (n < nsch) {
      tschnp1 <- schedule[n+1]  # end point is start of next schedue
    }
    else {
      tschnp1 <- tmax # end point is final time
    }
    tsim1 <- seq(tschn,tschnp1,by=0.1)
    np <- length(tsim1)
    xini[1] <- xini[1] + (step^(n+1))*dose  # add dose, toggle if step=-1 for infusion
    outn <- ode(y=xini,times=tsim1,func=func,parms=parms) # run for this segment
    xini <- outn[np,indy] # update xini for next segment
    nout <- dim(out)[1]
    out <- rbind(out,outn)  # append with duplicated time point
  }
  out <- as.data.frame(out)
  names(out)[c(1,(iout+1))] <- c('time','y') # output index given by iout
  createlookupmodel(unlist(out$time),unlist(out$y)) # lookup model based on time series
}

decaymodel <- function(t,x,parms) {
  ke <- parms$ke  # elimination rate
  dx <- -ke*x  # output
  list(dx)
}

stepmodel <- function(t,x,parms) {
#   dose <- parms$dose
#   ke <- parms$ke  # elimination rate
  dx <- 0  # output, steady state equals dose
  list(dx)
}

onecmptodemodel <- function(t,x,parms) {
  # central cmpt C, dC/dt = ka*A - k10*C
  # absorption cpmt A, dA/dt = -ka*A
  # here ke = k10
  # dose is per unit volume, normalised to bioavailability F=1
  # Cl = V*ke, typical values are Cl=10, V=25, so ke=0.4
  ke <- parms$ke  # elimination rate
  ka <- parms$ka  # absorption rate
  dx1 <- -ka*x[1]
  dx2 <- ka*x[1] - ke*x[2]  # output
  list(c(dx1,dx2))
}

twocmptodemodel <- function(t,x,parms) {
  # dose is per unit volume, normalised to bioavailability F=1
  # central cmpt C, dC/dt = k21*P - (k12+k10)*C = ka*P - (kp+ke)*C
  # peripheral cpmt P, dP/dt = k12*C - k21*P = kp*C - ka*P
  # here ke = k10, ka = k21, kp = k12
  ke <- parms$ke  # elimination rate k10
  ka <- parms$ka  # absorption rate from peripheral k21
  kp <- parms$kp # rate to peripheral k12
  dx1 <- kp*x[2] - ka*x[1]
  dx2 <- ka*x[1] - (kp + ke)*x[2]  # output
  list(c(dx1,dx2))
}

# function to assign phase factor to each compartment of cell cycle
# based on number of cmpts and fractions in g1, s, and m phase
getphasefactor <- function(frac,ncomp) {
  ind_frac <- floor(frac*ncomp)
  g2frac <- ncomp - sum(ind_frac)  # compute G2 fraction from others
  ind_frac <- c(g1=ind_frac[1], s=ind_frac[2],g2=g2frac, m=ind_frac[3]) # order G1, S, G2, M
  ind_phase <- factor(rep(c('G1','S','G2','M'),ind_frac))
  return(ind_phase)
}


# function to define cell model including pk model, drug effect, tumour dynamics
# output is list for ode solver containing func, xini and parameter list
# parameter list includes pk function which is defined separately
# also includes drugphase which specifies action of drug on cells
createcellmodel <- function(ncomp,tdoub,phasefrac,pkmodel,pkmodel2,drugcoef,drugcoef2) {
  ind_phase <- getphasefactor(phasefrac,ncomp)
  drugphase1 <- drugcoef$drugphase
  drugphase2 <- drugcoef2$drugphase
  nspecies <- 2*ncomp + 6
  xini <- rep(0,nspecies)
  xini[1:ncomp] <- rep(1, ncomp)/ncomp # equal volumes in each cmpt
  xini[(ncomp+1):(2*ncomp)] <- 0 # zero volume in dmg
  xini[2*ncomp+1] <- 1 # proliferating volume
  xini[2*ncomp+2] <- 0 # tumour radius starts at 0
  species <- paste('c',1:ncomp,sep='')
  species <- paste(species,ind_phase,sep='_')
  speciesdmg <- paste(species,'dmg',sep='_')
  species <- c(species,speciesdmg,'vp','rt','ap','ra','ap2','ra2') # names carry through to plot variables
  names(xini) <- species
  # pval is parameter list passed to ode solver
  pval <- list(ncomp=ncomp,
               tdoub=tdoub,
               drugfunc=pkmodel$drugfunc,
               drugfunc2=pkmodel2$drugfunc,
               drugphase1=drugphase1,
               drugphase2=drugphase2,
               ind_phase=ind_phase,
               drugcoef=drugcoef,
               drugcoef2=drugcoef2) # drug coefs
  list(func=odecellcycle,pval=pval,xini=xini)
}

# function to run the cell model
# calls drugfunc from pk model which is defined separately
odecellcycle <- function(t, x, parameters) {
  td <- parameters$tdoub
  nc <- parameters$ncomp
  ind_phase <- parameters$ind_phase
  drugphase1 <- parameters$drugphase1
  drugphase2 <- parameters$drugphase2
  drugfunc <- parameters$drugfunc
  drugfunc2 <- parameters$drugfunc2
  kapcoef <- parameters$drugcoef$kapcoef
  kdmgcoef <- parameters$drugcoef$kdmgcoef
  krepcoef <- parameters$drugcoef$krepcoef
  kapcoef2 <- parameters$drugcoef2$kapcoef
  kdmgcoef2 <- parameters$drugcoef2$kdmgcoef
  krepcoef2 <- parameters$drugcoef2$krepcoef

  if (drugphase1 == 'All') ind_drug <- !logical(length(ind_phase))
  else ind_drug <- ind_phase == drugphase1
  if (drugphase2 == 'All') ind_drug2 <- !logical(length(ind_phase))
  else ind_drug2 <- ind_phase == drugphase2
  
  nspecies <- 2*nc + 6
  dx <- rep(0,nspecies)
  grfac <- 2^(1/nc) # growth due to shift from previous compartment
  kfac <- log(2)/td/(grfac-1) 
  k <- kfac*rep(1,nc) # rates through phases
  # growth doubles after completion of cell cycle
  # kfac scales growth rate to give desired log(2)/tdoub
  kap <- rep(0,nc) # apoptosis, same as pa parameters
  kap[ind_drug] <- kapcoef*drugfunc(t)
  kdmg <- rep(0,nc) # rate of damage
  kdmg[ind_drug] <- kdmgcoef*drugfunc(t)
  
  kap2 <- rep(0,nc) # apoptosis, same as pa parameters
  kap2[ind_drug2] <- kapcoef2*drugfunc2(t)
  kdmg2 <- rep(0,nc) # rate of damage
  kdmg2[ind_drug2] <- kdmgcoef2*drugfunc2(t)

    # x tracks volume per cmpt, for number of cells scale by lamda
  # lamda[n] <- 1/nc/(grfac-1)/2^(1-n/nc)  volume per cell in cmpt i
  for (n in 1:nc) {
    nm1 <- ((n-2) %% nc) + 1
    vdmg <- (kdmg[n]+kdmg2[n])*x[n] - krepcoef*x[n+nc]  # damage, same repair rate
    dx[n] <- grfac*k[nm1]*x[nm1] - k[n]*x[n] - (kap[n] + kap2[n])*x[n] - vdmg
    dx[n+nc] <- vdmg # rate for damaged cells
    dx[2*nc+3] <- dx[2*nc+3] + kap[n]*x[n] # volume of cells killed drug1
    dx[2*nc+5] <- dx[2*nc+5] + kap2[n]*x[n] # volume of cells killed drug2
  }
  dx[2*nc+1] <- sum(dx[1:(2*nc)])  # total growth, includes damaged cells
  # dp = x(nspecies-2) # volume rate
  # p = sum(x(1:(nspecies-2))) # proliferating volume
  # dr/dt = dgr*(dp/dt)/p
  # max rate when dp/dt = log(2)/td*p is dr/dt= dgr*log(2)/td
  # normalise to dgr=1
  dx[2*nc+2] <- dx[2*nc+1]/x[2*nc+1]
  dx[2*nc+4] <- dx[2*nc+3]/x[2*nc+1] # missing growth due to cell death drug1
  dx[2*nc+6] <- dx[2*nc+5]/x[2*nc+1] # missing growth due to cell death drug2
  list(dx)
}

runcellmodel <- function(cellmodel,tmax) {
  tsim <- seq(0,tmax,by=0.1)
  out <- ode(y=cellmodel$xini,times=tsim,func=cellmodel$func,parms=cellmodel$pval)
  out <- as.data.frame(out)
  ind_phase <- cellmodel$pval$ind_phase
  ncomp <- cellmodel$pval$ncomp
  ycmpt <- out[,(1:(ncomp+1))] # growing cells
  ydmg <- out[,c(1,(ncomp+2):(2*ncomp+1))] # damaged cells 
  # calculate sums over each phase for plotting
  ycmptsum <- melt(ycmpt,id='time')
  ycmptsum$phase <- rep(ind_phase,each=length(ycmpt$time)) # add phase column
  ycmptsum <- aggregate(value~time + phase,ycmptsum, sum)
  ydmgsum <- melt(ydmg,id='time')
  ydmgsum$phase <- rep(ind_phase,each=length(ydmg$time)) # add phase column
  ydmgsum <- aggregate(value~time + phase,ydmgsum, sum)
  # calculate derived quantities
  ypro <- out[,c(1,2*ncomp+2)] # proliferating volume
  yrad <- out[,c(1,2*ncomp+3)] # tumour radius gain
  yapo <- out[,c(1,2*ncomp+4)] # volume of cells killed by drug
  yradapo <- out[,c(1,2*ncomp+5)] # radius loss due to cell death by drug
  yapo2 <- out[,c(1,2*ncomp+6)] # volume of cells killed by drug2
  yradapo2 <- out[,c(1,2*ncomp+7)] # radius loss due to cell death by drug2
  
  list(ycmptsum=ycmptsum,ydmgsum=ydmgsum,ypro=ypro,yrad=yrad,yapo=yapo,yradapo=yradapo,
       yapo2=yapo2,yradapo2=yradapo2,ind_phase=ind_phase)
}

# function to read file and return matrix of radius data
# rows are times, columns are experiments
loaddata <- function(filename, treat='untreated') {
  if (is.null(filename)) {
    filename <- 'www/ACCX16_TRT.csv'  # default
  }
  raddat <- read.table(filename, header = FALSE)
  if (dim(raddat)[2] > 1) {   # format has times in first row
    raddat <- as.data.frame(t(raddat))  # transpose
    # experiments start in column 2, time is in column 1
    nexp = ncol(raddat) - 1
    nt = nrow(raddat)
    treatnames <- list(treat)
  } else {    # format is same as ACCX16_TRT.csv
    diamdata <- read.csv(filename)
    diamdata$VOL<-as.numeric(as.character(diamdata$VOL))  # convert from factor
    # Dosing started from DAY 0 so take data from then
    diamdata <- diamdata[diamdata$DAYS >= 0, ]
    diamdata$TRT <- gsub('Control', 'untreated', diamdata$TRT)
    DIAM <- (6/pi*diamdata$VOL)^(1/3)  # compute diameters
    diamdata <- cbind(diamdata,DIAM=DIAM)  # add column diameter
    treatnames <- as.character(unique(diamdata$TRT))
    curtreatind <- which(diamdata$TRT == treat) # indices of diamdata for this treatment
    curidind <- as.character(unique(diamdata$ID[curtreatind])) # IDs with this treatment
    timevec <- 24*diamdata[diamdata$ID == curidind[1],]$DAYS  # time points associated with first ID
    diamdata <- diamdata[curtreatind,]
    nt <- length(timevec)
    ntot <- length(curtreatind)  # total number of rows
    nexp <- ntot/nt  # number of experiments
    expmat <- matrix(diamdata$DIAM,nrow=nt,ncol=nexp) # convert vector to matrix
    raddat <- data.frame(timevec, expmat)
  }
  colnames(raddat) <- c("hr", paste("exp",1:nexp, sep = ""))
  rownames(raddat) <- paste("t",1:nt,sep = "")
  expmat <- as.matrix(raddat[,2:ncol(raddat)])
  fit <- lm(expmat ~ raddat$hr) # fit over all experiments
  lincoefs <- rowMeans(fit$coefficients)
  names(lincoefs) <- c('intercept','slope')
  list(radiusdata=raddat,lincoefs=lincoefs,treatnames=treatnames)
}

# function to calculate standard error of vector x
std.err <- function(x){
  result <- sqrt(var(x)/length(x))
  return(result)
}

ploterrbars <- function(raddat) {
  # add error bars
  raddat.t <- raddat$hr
  raddat.mn <- apply(raddat,1,mean)
  raddat.se <- apply(raddat,1,std.err) # standard errors over rows of matrix
  errdat <- data.frame(t=raddat.t,mn=raddat.mn,se=raddat.se)
  g <- g + geom_errorbar(data=errdat,aes(x=t,y=mn,ymin=mn-2*se, ymax=mn+2*se, color=NULL), size=1,width=2)
  
}

# function to plot pk
plotpk <- function(drugfunc,type,tmax) {
  tpk <- seq(0,tmax,by=0.1)
  ypk <- drugfunc(tpk)
  pkres <- data.frame(time=tpk,y=ypk)
  if (type==1) {
    g <- geom_line(data=pkres,aes(x=time,y=y),color='red',size=2)
  }
  if (type==2) {
    g <- geom_line(data=pkres,aes(x=time,y=y),color='blue',size=2)
  }
  return(g)
}

plotresults <- function(results,type,phaseflag) {
  varnames <- c('ycmptsum','ydmgsum')
  y <- results[[varnames[as.numeric(type)]]]
  type <- as.numeric(type)
  g <- ggplot()
  if (1 %in% type) {
    g <- g + geom_line(data=results$ycmptsum,aes(x=time,y=value,color=phase),size=1,alpha=0.8)
  }
  if (2 %in% type) {
    g <- g + geom_line(data=results$ydmgsum,aes(x=time,y=value,color=phase),size=1,alpha=0.8,
                       linetype='dashed')
  }
  g <- g + geom_line(data=results$ypro,aes(x=time,y=vp),size=2,alpha=0.8)
  g <- g + scale_color_manual(values=c("blue","green","yellow","red")) #g1,g2,m,s
  g <- g + labs(x="time (hrs)", y="volume (relative)")
  
}


shinyServer(function(input, output, session) {
  
  raddatlist <- loaddata('www/ACCX16_TRT.csv', 'untreated')  # load default data for overlay
  updateSelectInput(session, "treatment", choices = raddatlist$treatnames)
  
  # set values for parameter input from stored pkparlist
  # when model is toggled, set input value to value stored in list
  
  updatesliders <- reactive({
    tmax <- input$tmax
    updateSliderInput(session,'rangeAxis',max=tmax,value=c(0,tmax)) # set tmax on sliders
    updateSliderInput(session,'pkrangeAxis',max=tmax,value=c(0,tmax)) # does not update slider when not in view
    updateSliderInput(session,'pkrangeAxis2',max=tmax,value=c(0,tmax)) # does not update slider when not in view
    updateSliderInput(session,'radrangeAxis',max=tmax,value=c(0,tmax)) 
  })
  
  readpkmodel <- function(pknum=1,pkFile=NULL) {
    tdoub <- input$tdoub
    tmax <- input$tmax
    updatesliders() 
    if (pknum == 1) {
      schedule <- input$sch1
      wks <- input$wks1
      dosein <- input$dose1
      pktype <- input$pktype1
      ka <- input$pkka1
      ke <- input$pkke1
      kp <- input$pkkp1
    }
    else if (pknum == 2) {
      schedule <- input$sch2
      wks <- input$wks2
      dosein <- input$dose2
      pktype <- input$pktype2
      ka <- input$pkka2
      ke <- input$pkke2
      kp <- input$pkkp2
    }
    schedule <- as.numeric(unlist(strsplit(schedule,","))) # string to vector
    wkshr <- seq(0,wks-1)*24*7  # start of each week in hours
    schedule <- as.vector(sapply(wkshr,function(x) x+schedule))
    tpkoff <- schedule[1]  # just take first number for now
    if (pktype=='K-PD') {  # simple decay
      pkpar <- list(dose=dosein, schedule=schedule,ke=ke)
      pkmodelfunc <- createlookupfromode(decaymodel,tmax,xini=0,pkpar,iout=1) 
      pkmodellist <- list(drugfunc=pkmodelfunc,pkpar=pkpar)
    }
    else if (pktype=='Step') {  # constant infusion
      pkpar <- list(dose=dosein, schedule=schedule)
      pkmodelfunc <- createlookupfromode(stepmodel,tmax,xini=0,pkpar,iout=1,step=-1) 
      pkmodellist <- list(drugfunc=pkmodelfunc,pkpar=pkpar)
    }
    else if (pktype=='1-cmpt') {
      pkpar <- list(dose=dosein, schedule=schedule,ke=ke,ka=ka)
      pkmodelfunc <- createlookupfromode(onecmptodemodel,tmax,xini=c(0,0),pkpar,iout=2)
      pkmodellist <- list(drugfunc=pkmodelfunc,pkpar=pkpar)
    }
    else if (pktype=='2-cmpt') {
      pkpar <- list(dose=dosein, schedule=schedule,ke=ke,ka=ka,kp=kp)
      pkmodelfunc <- createlookupfromode(twocmptodemodel,tmax,xini=c(0,0),pkpar,iout=2)
      pkmodellist <- list(drugfunc=pkmodelfunc,pkpar=pkpar)
    }
    else if (!is.null(pkFile)) {
      pkmodelcode <- source(pkFile$datapath)  # code containing model
      pkmodellist <- pkmodelcode$value()  # list with function and params
    }
    return(pkmodellist)
  }
  
  runpkmodel1 <- eventReactive(input$runpkButton1, {
    # create lookup model by running the model
    tmax <- input$tmax
    tsim <- seq(0, tmax, by = 0.1) # for lookup function
    pkFile <- input$sourcepkfile
    pkmodellist <- readpkmodel(pknum=1,pkFile) # list with function and params
    ysim <- pkmodellist$drugfunc(tsim)  # perform simulation
    auc <- simpsonvec(tsim,ysim)
    pkmodelfunc <- createlookupmodel(tsim,ysim)
    return(list(pkmodelfunc=pkmodelfunc,auc=auc))
  })
  
  runpkmodel2 <- eventReactive(input$runpkButton2, {
    # create lookup model by running the model
    tmax <- input$tmax
    tsim <- seq(0, tmax, by = 0.1) # for lookup function
    pkFile <- input$sourcepkfile2
    pkmodellist <- readpkmodel(pknum=2,pkFile) # list with function and params
    ysim <- pkmodellist$drugfunc(tsim)  # perform simulation
    auc <- simpsonvec(tsim,ysim)
    pkmodelfunc <- createlookupmodel(tsim,ysim)
    return(list(pkmodelfunc=pkmodelfunc,auc=auc))
  })
  
  runmodel <- eventReactive(input$runButton, {
    tdoub <- input$tdoub
    tmax <- input$tmax
    updatesliders() 
    ncomp <- as.numeric(input$ncompstr) 
    G2phase <- 1 - input$G1phase - input$Sphase - input$Mphase
    phasefrac <- c(input$G1phase,input$Sphase,input$Mphase) # G1, S, M
    tsim <- seq(0, tmax, by = 0.1) # for lookup function
    pkmodel <- readpkmodel(input$sourcepkfile,pknum=1) # runs without pk run button
    pkmodel2 <- readpkmodel(input$sourcepkfile2,pknum=2) 
    drugcoef <- list(kapcoef=input$kap, # apoptosis rate
                     kdmgcoef=input$kdmg, #damage
                     krepcoef=input$krep,   # repair
                     drugphase1=input$drugphase1)
    drugcoef2 <- list(kapcoef=input$kap2, # apoptosis rate
                      kdmgcoef=input$kdmg2, #damage
                      krepcoef=input$krep2,   # repair
                      drugphase2=input$drugphase2)
    cellmodel <- createcellmodel(ncomp = ncomp,   # number of compartments
                                 tdoub = tdoub,    # cell doubling time
                                 phasefrac = phasefrac,
                                 pkmodel = pkmodel,
                                 pkmodel2 = pkmodel2,
                                 drugcoef = drugcoef,
                                 drugcoef2 = drugcoef2)
    withProgress(message = 'Cycling ...', value = 0.1, {
                 results <- runcellmodel(cellmodel,tmax)
    })
  })
  
  output$volPlot <- renderPlot({
    results <- runmodel()
    g <- plotresults(results,input$plotVariables,input$phaseSum)
    g <- g + coord_cartesian(xlim = input$rangeAxis) 
    if (input$logAxis==TRUE) {g <- g + scale_y_log10()}
    print(g)
  })
  
  output$pkPlot <- renderPlot({
    tdoub <- input$tdoub
    tmax <- input$tmax
    tsim <- seq(0, tmax, by = 0.1) # for lookup function
    pkmodelfunc <- runpkmodel1()$pkmodelfunc #readpkmodel(input$sourcepkfile,pknum=1) 
    #pkmodel2 <- runpkmodel2() #readpkmodel(input$sourcepkfile,pknum=2) 
    plotsel <- as.numeric(input$pkplotVariables) # model 1 and/or 2
    g <- ggplot()
    g <- g + plotpk(pkmodelfunc,type=1,tmax)
    g <- g + coord_cartesian(xlim = input$pkrangeAxis)
    if (input$pklogAxis==TRUE) {g <- g + scale_y_log10()}
    g <- g + labs(x="time (hrs)", y="concentration")
    print(g)
  })
  
  output$pkPlot2 <- renderPlot({
    tdoub <- input$tdoub
    tmax <- input$tmax
    tsim <- seq(0, tmax, by = 0.1) # for lookup function
    pkmodelfunc2 <- runpkmodel2()$pkmodelfunc #readpkmodel(input$sourcepkfile,pknum=2) 
    plotsel2 <- as.numeric(input$pkplotVariables2) # model 1 and/or 2
    g <- ggplot()
    g <- g + plotpk(pkmodelfunc2,type=2,tmax)
    g <- g + coord_cartesian(xlim = input$pkrangeAxis2)
    if (input$pklogAxis2==TRUE) {g <- g + scale_y_log10()}
    g <- g + labs(x="time (hrs)", y="concentration")
    print(g)
  })
  
  updatetreatments <- observeEvent(input$dataFile, {
    inFile <- input$dataFile
    treat <- input$treatment
    raddatlist <- loaddata(inFile$datapath, treat)
    updateSelectInput(session, "treatment", choices = raddatlist$treatnames)
    updateCheckboxGroupInput(session, "tumplotOverlay", selected = 'overlay')
  })
  
  output$radPlot <- renderPlot({
    rini <- input$rini
    dgr <- input$dgr
    results <- runmodel()
    g <- ggplot() + theme(legend.position="none")
    g <- g + labs(x="time (hrs)", y="radius (mm)")
    if ('overlay' %in% input$tumplotOverlay) {  #(input$showOverlay==TRUE) {
      inFile <- input$dataFile
      treat <- input$treatment
      raddatlist <- loaddata(inFile$datapath, treat)
      raddat <- raddatlist$radiusdata
      raddat_melt <- melt(raddat,id='hr')
      g <- g + geom_point(data=raddat_melt,aes(x=hr,y=value,color=variable),size=2,alpha=0.5)
      g <- g + geom_line(data=raddat_melt,aes(x=hr,y=value,color=variable),size=1,alpha=0.2)

      if ('linear fit' %in% input$tumplotOverlay) { # if (input$showLinearFit==TRUE) {
        g <- g + geom_smooth(data=raddat_melt,aes(x=hr,y=value,color=NULL), method = "lm",
                             alpha=0.8,size=1) # fits entire data set, doesn't account for grouped
      }
    }
    time <- results$yrad$time
    control <- rini + dgr*log(2)/input$tdoub*time # control
    radgr <- rini + dgr*results$yrad$rt  # scale by dgr and add initial radius
    radapo <- dgr*results$yradapo$ra # scale by dgr
    radapo2 <- dgr*results$yradapo2$ra2 # scale by dgr
    raddmg <- control - radgr - radapo - radapo2 # radius lost to damage
    if (1 %in% input$tumplotVariables) {        # show control
      radapo <- control - radapo
      radapo2 <- radapo - radapo2
      raddmg <- radapo2 - raddmg
      g <- g + geom_line(aes(x=time,y=control),color='grey',size=1,alpha=1)
      if (2 %in% input$tumplotVariables) {      # add ribbons for drug effect
        g <- g + geom_ribbon(aes(ymin=radapo, ymax=control, x=time), 
                             fill='red',color='red',alpha=0.2)
        g <- g + geom_ribbon(aes(ymin=radapo2, ymax=radapo, x=time), 
                             fill='blue',color='blue',alpha=0.2)
        g <- g + geom_ribbon(aes(ymin=radgr, ymax=radapo2, x=time), 
                             fill='green',color='green',alpha=0.2)
      }
    } 
    if (2 %in% input$tumplotVariables) {          # show drug effect
      g <- g + geom_line(aes(x=time,y=radapo),color='red',size=2,alpha=0.5)
      g <- g + geom_line(aes(x=time,y=radapo2),color='blue',size=2,alpha=0.5)
      g <- g + geom_line(aes(x=time,y=raddmg),color='green',size=2,alpha=0.5)
    }
    g <- g + geom_line(aes(x=time,y=radgr),color='black',size=2)
    g <- g + labs(x="time (hrs)", y="radius (mm)")
    g <- g + coord_cartesian(xlim = input$radrangeAxis)
    print(g)
  })
 
  output$resultsText <- renderTable({ 
    rini <- input$rini
    dgr <- input$dgr
    results <- runmodel()
    tsim <- results$yrad[,1]
    tf <- tail(tsim,1)
    radgain <- dgr*tail(results$yrad[2],1)  # final minus initial radius
    raddeath1 <- dgr*tail(results$yradapo[2],1) # loss due to cell death drug1, scale by dgr
    raddeath2 <- dgr*tail(results$yradapo2[2],1) # loss due to cell death drug2
    raddeathtot <- raddeath1 + raddeath2
    radmax <- dgr*log(2)/input$tdoub*tf # max possible gain
    delrad <- radmax-radgain
    raddam <- delrad - raddeathtot
    sumvar <- c('Radius gain','Control gain','Radius loss','Death PK1','% of total',
                'Death PK2','% of total',
                'Damage','% of total')
    sumval <- as.numeric(c(radgain,radmax,delrad,
                           raddeath1,100*(raddeath1/delrad),
                           raddeath2,100*(raddeath2/delrad),
                           raddam,100*(raddam/delrad)))
    sumunit <- c(' mm',' mm', ' mm',' mm', 'percent', ' mm', 'percent', ' mm', 'percent') 
    tablesummary <- data.frame(quantity=sumvar,value=sumval,unit=sumunit)
    tablesummary <- format(tablesummary,digits=3,width=8,format='f')
  },include.rownames=FALSE,include.colnames=FALSE)

  output$resultsTextDrug1 <- renderTable({ 
    results <- runmodel()
    auc1 <- runpkmodel1()$auc
    raddeath1 <- tail(results$yradapo[2],1) # loss due to cell death drug1
    radperdose1 <- raddeath1/auc1
    sumvar <- c('PK1 AUC','Death/AUC')
    sumval <- as.numeric(c(auc1,radperdose1))
    sumunit <- c('dose-hr','mm/dose/hr')
    tablesummary <- data.frame(quantity=sumvar,value=sumval,unit=sumunit)
    tablesummary <- format(tablesummary,digits=3,width=8,format='f')
  },include.rownames=FALSE,include.colnames=FALSE)
  
  output$resultsTextDrug2 <- renderTable({ 
    results <- runmodel()
    auc2 <- runpkmodel2()$auc
    raddeath2 <- tail(results$yradapo2[2],1) # loss due to cell death drug2
    radperdose2 <- raddeath2/auc2
    sumvar <- c('PK2 AUC','Death/AUC')
    sumval <- as.numeric(c(auc2,radperdose2))
    sumunit <- c('dose-hr','mm/dose/hr')
    tablesummary <- data.frame(quantity=sumvar,value=sumval,unit=sumunit)
    tablesummary <- format(tablesummary,digits=3,width=8,format='f')
  },include.rownames=FALSE,include.colnames=FALSE)
  
  output$fitcoefsText <- renderTable({ 
    if (!('overlay' %in% input$tumplotOverlay)) {
      return(NULL)  # only show table if overlay is selected
    }
    raddatlist <- loaddata(input$dataFile$datapath,input$treatment)
    intercept <- raddatlist$lincoefs[[1]] # estimate for r0 if control
    slope <- raddatlist$lincoefs[[2]]
    tdoub <-input$tdoub  # model doubling time
    dgrest <- slope*tdoub/log(2)
    rini <- input$rini  # model r0
    delr <- rini - intercept # estimate for drug effect if control slope
    sumvar <- c('Slope data','Doubling time model','Growing layer data',
                'Initial radius data','Delta radius')
    sumval <- as.numeric(c(slope,tdoub,dgrest,intercept,delr))
    sumunit <- c('mm/hr','hr','mm','mm','mm')
    tablesummary <- data.frame(quantity=sumvar,value=sumval,unit=sumunit)
    tablesummary <- format(tablesummary,digits=3,width=8,format='f')
  },include.rownames=FALSE,include.colnames=FALSE)
  
  
  output$G2phase <- renderText({
    g1phase <- as.character(1 - input$G1phase - input$Sphase - input$Mphase)
    c('G2 phase ____________',g1phase)
    })
  
  output$repairrate <- renderText({
    krep <- as.character(input$krep)
    c('Repair rate _________',krep, ' (as drug 1)')
  })
  
    output$tdoubstdev <- renderText({
      ncomp <- as.numeric(input$ncompstr) 
      tdoubstdev <- 1/sqrt(ncomp)
      tdoubstdev <- format(tdoubstdev,digits = 2)
      c('Proportional standard deviation of doubling times: 1/sqrt(N) = ',tdoubstdev)
    })
  
    output$discret <- renderText({
      ncomp <- as.numeric(input$ncompstr) 
      discret <- 1/ncomp
      discret <- format(discret,digits = 3)
      c('Discretisation of cell cycle: 1/N = ',discret)
    })
    
  output$downloadData <- downloadHandler(filename = function() { 
      paste('results', '.csv', sep='') 
    },
    content = function(file) {
      rini <- input$rini
      dgr <- input$dgr
      results <- runmodel()
      results$yrad$rt <- rini + dgr*results$yrad$rt # scale for initial and growing layer thickness
      write.csv(results$yrad, file)
    }
  )
  
  # Whenever a field is filled, aggregate all form data
  formData <- reactive({
    fields <- names(input)
    data <- sapply(fields, function(x) input[[x]])
  })
  
  output$savefile <- downloadHandler(
    filename = function() { 
      'modelsettings.csv'
    },
    content = function(file) {
      data <- as.data.frame(t(formData()))
      write.csv(t(data),file)
    }
  )
  
  observeEvent(input$readfile, {
    loadInputData()
  })
  
  as.numeric.factor <- function(x) {as.numeric(gsub(",", "", levels(x)[x]))}

  loadInputData <- function() {  # read input data from file
    data <- read.csv(input$readfile$datapath)
    inputnames <- as.character(data[,1])
    inputvalues <- as.numeric.factor(data[,2])
    np <- length(inputvalues)
    for (n in 1:np) {
      if (!is.na(inputvalues[n]) & is.numeric(inputvalues[n])) {
        updateNumericInput(session, inputnames[n], value = inputvalues[n])
      }
    }
    # Change values for selected input$drugphase, which is not numeric
    iphase <- grep('drugphase1',data[,1])
    updateSelectInput(session, "drugphase1",
                      selected = as.character(data[iphase,2])
    )
    iphase2 <- grep('drugphase2',data[,1])
    updateSelectInput(session, "drugphase2",
                      selected = as.character(data[iphase2,2])
    )
    
  }
  
  
})