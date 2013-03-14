library(shiny)

if (try(library(ggplot2)) == "try-error") {install.packages("ggplot2")}
library(ggplot2)

if (try(library(devtools)) == "try-error") {install.packages("devtools")}
library(devtools)

if (try(library(earlywarnings)) == "try-error") {install_github(repo = "earlywarnings-R", username = "earlywarningtoolbox", subdir = "earlywarnings", ref = "master")}
library(earlywarnings)


# Simulated data
simulateddata <- read.csv("fold_simulated_data.csv")

# Real data
climatedata <- read.csv("climate_data.csv")

#userdata <- climatedata
 
shinyServer(function(input, output) {
 
    datasetInput <- reactive({
        switch(input$timeseries,
        "simulated - overharvested resource" = simulateddata,
        "real-world - climate data" = climatedata)
    })
    
  output$plot <- reactivePlot(function() {

    qda_RShiny(datasetInput(), winsize = input$winsize, detrending = input$detrending, logtransform = input$logtransform, interpolate = input$interpolate, analysis = input$analysis, cutoff = input$cutoff, detection.threshold = input$detection.threshold, grid.size = input$grid.size)

  }, height=500)
})

qda_RShiny <- function(timeseries, param = NULL, winsize = 50, detrending=c("no","gaussian","linear","first-diff"), bandwidth=NULL, boots = 100, s_level = 0.05, cutoff=0.05, detection.threshold = 0.002, grid.size = 50, logtransform=FALSE, interpolate=FALSE, analysis = c("Indicator trend analysis", "Trend significance analysis","Potential analysis")){
  
  timeseries <- data.matrix(timeseries)
  
  if (analysis == "Indicator trend analysis") { 
    
    message("Indicator trend analysis")
    generic_RShiny(timeseries,winsize,detrending,bandwidth,logtransform,interpolate,AR_n=FALSE,powerspectrum=FALSE)
    
  } else if (analysis == "Trend significance analysis") {
    
    message("Trend significance analysis")
    s <- surrogates_RShiny(timeseries,winsize,detrending,bandwidth,boots,s_level,logtransform,interpolate)
    print(s)
    
  } else if (analysis == "Potential analysis") { 
    
    if (ncol(timeseries) == 1) {
      observations <- as.vector(timeseries[,1])
      timepoints <- 1:nrow(timeseries)
    } else if (ncol(timeseries) > 1) {
      observations <- as.vector(timeseries[,2])
      timepoints <- as.vector(timeseries[,1])
    }
    
    message("Potential analysis")
    p <- movpotential_ews(observations, timepoints, detection.threshold = detection.threshold, grid.size = grid.size, plot.cutoff = cutoff)
    print(p)
    
  }
  #   print("Sensitivity of trends")
  #   s <- sensitivity_RShiny(timeseries,winsizerange=c(25,75),incrwinsize,detrending=detrending, bandwidthrange=c(5,100),incrbandwidth,logtransform=FALSE,interpolate=FALSE)
}

# generic_Rshiny for estimating only AR1 and Variance in moving windows with various options for pretreating the data
# 26 Feb 2013

generic_RShiny<-function(timeseries, winsize = 50, detrending = c("no", "gaussian", "linear", "first-diff"), bandwidth=NULL, logtransform=FALSE, interpolate=FALSE, AR_n = FALSE, powerspectrum = FALSE){  
  
  #   require(lmtest)
  #   require(nortest)
  #   require(stats)
  #   require(som)
  #   require(Kendall)
  #   require(KernSmooth)
  #   require(moments)
  
  #timeseries<-ts(timeseries)
  timeseries<-data.matrix(timeseries) #strict data-types the input data as tseries object for use in later steps
  if (ncol(timeseries) == 1){
    Y=timeseries
    timeindex=1:dim(timeseries)[1]
  }else if(dim(timeseries)[2]==2){
    Y<-timeseries[,2]
    timeindex<-timeseries[,1]
  }else{
    warning("not right format of timeseries input")
  }
  #return(timeindex)
  # Interpolation
  if (interpolate){
    YY<-approx(timeindex,Y,n=length(Y),method="linear")
    Y<-YY$y
  }else{
    Y<-Y}
  
  # Log-transformation
  if (logtransform){
    Y<-log(Y+1)}
  
  # Detrending  
  detrending<-match.arg(detrending)	
  if (detrending=="gaussian"){
    if (is.null(bandwidth)){
      bw<-round(bw.nrd0(timeindex))}else{
        bw<-round(length(Y)*bandwidth/100)}
    smYY<-ksmooth(timeindex,Y,kernel="normal",bandwidth=bw,range.x=range(timeindex),x.points=timeindex)
    if (timeindex[1]>timeindex[length(timeindex)]){
      nsmY<-Y-rev(smYY$y)
      smY<-rev(smYY$y) 
    }else{  	nsmY<-Y-smYY$y
             smY<-smYY$y}
  }else if(detrending=="linear"){
    nsmY<-resid(lm(Y~timeindex))
    smY<-fitted(lm(Y~timeindex))
  }else if(detrending=="first-diff"){
    nsmY<-diff(Y)
    timeindexdiff<-timeindex[1:(length(timeindex)-1)]
  }else if(detrending=="no"){
    smY<-Y
    nsmY<-Y
  }
  
  
  # Rearrange data for indicator calculation
  mw<-round(length(Y)*winsize/100)
  omw<-length(nsmY)-mw+1 ##number of moving windows
  low<-6
  high<-omw
  nMR<-matrix(data=NA,nrow=mw,ncol=omw)
  x1<-1:mw
  for (i in 1:omw){ 	 
    Ytw<-nsmY[i:(i+mw-1)]
    nMR[,i]<-Ytw}
  
  # Calculate indicators
  nARR<-numeric()
  nSD<-numeric()
  
  nSD<-apply(nMR, 2, sd, na.rm = TRUE)
  for (i in 1:ncol(nMR)){
    nYR<-ar.ols(nMR[,i],aic= FALSE, order.max=1, dmean=FALSE, 		intercept=FALSE)
    nARR[i]<-nYR$ar
  }
  
  nVAR=sqrt(nSD)
  
  # Estimate Kendall trend statistic for indicators
  timevec<-seq(1,length(nARR))
  KtAR<-cor.test(timevec,nARR,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
  KtVAR<-cor.test(timevec,nVAR,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
  
  # Plotting
  # Generic Early-Warnings
  #dev.new()
  par(mar=(c(1,2,0.5,2)+0),oma=c(2,2,2,2),mfrow=c(4,1))  
  plot(timeindex,Y,type="l",ylab="",xlab="",xaxt="n",lwd=2,las=1,xlim=c(timeindex[1],timeindex[length(timeindex)]))
  legend("bottomleft","data",,bty = "n")
  if(detrending=="gaussian"){
    lines(timeindex,smY,type="l",ylab="",xlab="",xaxt="n",lwd=2,col=2,las=1,xlim=c(timeindex[1],timeindex[length(timeindex)]))
  }
  if(detrending=="no"){
    plot(c(0,1),c(0,1),ylab="",xlab="",yaxt="n",xaxt="n",type="n",las=1)
    text(0.5,0.5,"no detrending - no residuals")
  }else if (detrending=="first-diff"){
    limit<-max(c(max(abs(nsmY))))
    plot(timeindexdiff,nsmY,ylab="",xlab="",type="l",xaxt="n",lwd=2,las=1,ylim=c(-	limit,limit),xlim=c(timeindexdiff[1],timeindexdiff[length(timeindexdiff)]))
    legend("bottomleft","first-differenced",bty = "n")		}else{
      limit<-max(c(max(abs(nsmY))))
      plot(timeindex,nsmY,ylab="",xlab="",type="h",xaxt="n",las=1,lwd=2,ylim=c(-	limit,limit),xlim=c(timeindex[1],timeindex[length(timeindex)]))
      legend("bottomleft","residuals",bty = "n")}
  plot(timeindex[mw:length(nsmY)],nARR,ylab="",xlab="",type="l",xaxt="n",col="green",lwd=2,las=1,xlim=c(timeindex[1],timeindex[length(timeindex)])) #3
  legend("bottomright",paste("trend ",round(KtAR$estimate,digits=3)),bty = "n")
  legend("bottomleft","autocorrelation",bty = "n")
  plot(timeindex[mw:length(nsmY)],nVAR,ylab="",xlab="",type="l",col="blue", lwd=2, las=1,xlim=c(timeindex[1],timeindex[length(timeindex)]))
  legend("bottomright",paste("trend ",round(KtVAR$estimate,digits=3)),bty = "n")
  legend("bottomleft","variance",bty = "n")
  mtext("time",side=1,line=2,cex=0.8)
  mtext("Generic Early-Warnings: Autocorrelation - Variance",side=3,line=0.2, outer=TRUE)#outer=TRUE print on the outer margin
  
  # Output
  out<-data.frame(timeindex[mw:length(nsmY)],nARR,nSD)
  colnames(out)<-c("timeindex","ar1","sd")
  #return(out)
  
}


# surrogates_Rshiny for estimating significance of trends for variance and autocorrelation
# 6 March 2013

surrogates_RShiny<-function(timeseries,winsize=50,detrending=c("no","gaussian","linear","first-diff"),bandwidth=NULL,boots=50,s_level=0.05,logtransform=FALSE,interpolate=FALSE){
  
  timeseries<-data.matrix(timeseries)
  if (dim(timeseries)[2]==1){
    Y=timeseries
    timeindex=1:dim(timeseries)[1]
  }else if(dim(timeseries)[2]==2){
    Y<-timeseries[,2]
    timeindex<-timeseries[,1]
  }else{
    warning("not right format of timeseries input")
  }
  
  # Interpolation
  if (interpolate){
    YY<-approx(timeindex,Y,n=length(Y),method="linear")
    Y<-YY$y
  }else{
    Y<-Y}
  
  # Log-transformation
  if (logtransform){
    Y<-log(Y+1)}
  
  # Detrending	
  detrending<-match.arg(detrending)	
  if (detrending=="gaussian"){
    if (is.null(bandwidth)){
      bw<-round(bw.nrd0(timeindex))}else{
        bw<-round(length(Y)*bandwidth)/100}
    smYY<-ksmooth(timeindex,Y,kernel=c("normal"), bandwidth=bw, range.x=range(timeindex),n.points=length(timeindex))
    nsmY<-Y-smYY$y
    smY<-smYY$y
  }else if(detrending=="linear"){
    nsmY<-resid(lm(Y~timeindex))
    smY<-fitted(lm(Y~timeindex))
  }else if(detrending=="first-diff"){
    nsmY<-diff(Y)
    timeindexdiff<-timeindex[1:(length(timeindex)-1)]
  }else if(detrending=="no"){
    smY<-Y
    nsmY<-Y
  }
  
  
  # Rearrange data for indicator calculation
  mw<-round(length(Y)*winsize)/100
  omw<-length(nsmY)-mw+1
  low<-6
  high<-omw
  nMR<-matrix(data=NA,nrow=mw,ncol=omw)
  for (i in 1:omw){
    Ytw<-nsmY[i:(i+mw-1)]
    nMR[,i]<-Ytw}
  # Estimate indicator
  
  indic_ar1<-apply(nMR,2,function(x){nAR1<-ar.ols(x,aic= FALSE, order.max=1,dmean=FALSE,intercept=FALSE)
                                     nAR1$ar})
  
  indic_var<-apply(nMR,2,var)
  
  # Calculate trend statistics
  timevec<-seq(1,length(indic_ar1))
  Kt_ar1<-cor.test(timevec,indic_ar1,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
  Ktauestind_ar1orig<-Kt_ar1$estimate
  
  Kt_var<-cor.test(timevec,indic_var,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
  Ktauestind_varorig<-Kt_var$estimate
  
  # Fit ARMA model based on AIC
  #   arma=matrix(,4,5)
  #   for (ij in 1:4){
  #     for (jj in 0:4){
  ARMA<-arima(nsmY, order = c(1,0,0),include.mean = FALSE)
  #       arma[ij,jj+1]=ARMA$aic		
  # 		print(paste("AR","MA", "AIC"),quote=FALSE)
  # 		print(paste(ij,jj,ARMA$aic),zero.print=".",quote=FALSE)
# }
# }

# Simulate ARMA(p,q) model fitted on residuals
#   ind=which(arma==min(arma),arr.ind=TRUE)
#   ARMA<-arima(nsmY, order = c(ind[1],0,ind[2]-1),include.mean = FALSE)

Ktauestind_ar1 <- numeric()
Ktauestind_var <- numeric()

for (jjj in 1:boots){
  x=arima.sim(n = length(nsmY), list(ar = ARMA$coef[1], ma = 0, sd=sqrt(ARMA$sigma2)))
  
  ## Rearrange data for indicator calculation
  nMR1<-matrix(data=NA,nrow=mw,ncol=omw)
  for (i in 1:omw){   
    Ytw<-x[i:(i+mw-1)]
    nMR1[,i]<-Ytw}
  
  # Estimate indicator
  
  indic_ar1<-apply(nMR1,2,function(x){nAR1<-ar.ols(x,aic= FALSE, order.max=1,dmean=FALSE,intercept=FALSE)
                                      nAR1$ar})
  
  indic_var<-apply(nMR1,2,var)
  
  # Calculate trend statistics
  timevec<-seq(1,length(indic_ar1))
  Kt_ar1<-cor.test(timevec,indic_ar1,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
  Ktauestind_ar1[jjj]<-Kt_ar1$estimate
  
  Kt_var<-cor.test(timevec,indic_var,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
  Ktauestind_var[jjj]<-Kt_var$estimate
  
}

# Estimate probability of false positive
q_ar1<-sort(Ktauestind_ar1,na.last=NA)
Kpos_ar1<-max(which(Ktauestind_ar1orig>q_ar1),na.rm=TRUE)
p<-(boots+1-Kpos_ar1)/boots
print(paste('significance autocorrelation p = ',p,' estimated from ',boots,' surrogate ARMA timeseries'))

q_var<-sort(Ktauestind_var,na.last=NA)
Kpos_var<-max(which(Ktauestind_varorig>q_var),na.rm=TRUE)
p<-(boots+1-Kpos_var)/boots
print(paste('significance variance p = ',p,' estimated from ',boots,' surrogate ARMA timeseries'))

# Plotting
layout(matrix(1:2,1,2))
par(font.main=10,mar=(c(4.6,3.5,0.5,2)+0.2),mgp=c(2,1,0),oma=c(0.5,0.5,2,0),cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
hist(Ktauestind_ar1,freq=TRUE,nclass=20,xlim=c(-1,1),col="green",main=NULL,xlab="Surrogate trend estimates",ylab="occurrence")#,ylim=c(0,boots))
abline(v=q_ar1[s_level*boots],col="red",lwd=2)
abline(v=q_ar1[(1-s_level)*boots],col="red",lwd=2)
points(Ktauestind_ar1orig,0,pch=21, bg="black", col = "black", cex=4)
title("Autocorrelation",cex.main=1.3)

hist(Ktauestind_var,freq=TRUE,nclass=20,xlim=c(-1,1),col="blue",main=NULL,xlab="Surrogate trend estimates",ylab="occurrence")#,ylim=c(0,boots))
abline(v=q_var[s_level*boots],col="red",lwd=2)
abline(v=q_var[(1-s_level)*boots],col="red",lwd=2)
points(Ktauestind_varorig,0,pch=21, bg="black", col = "black", cex=4)
title("Variance" ,cex.main=1.3)
}

