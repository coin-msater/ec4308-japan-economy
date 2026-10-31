library(tidyverse)
library(lubridate)


# Read the four horizon datasets written by the preprocessing script
df_h1  <- read_csv("data/processed/h1.csv")
df_h3  <- read_csv("data/processed/h3.csv")
df_h6  <- read_csv("data/processed/h6.csv")
df_h12 <- read_csv("data/processed/h12.csv")



make_series=function(df_h,lag){
  cbind(c(df_h$x1_lag0, tail(df_h$target,lag))) #one-column matrix with x1 only
}


#Inputs for the function:

#1) Y: one-column matrix with the x1 series

#2) lag - the forecast horizon

#3) max_p - the largest lag order BIC may choose (candidates are AR(1) ... AR(max_p))

runAR=function(Y,lag,max_p=12,p_fixed=NULL){
  
  Y=Y[1:(nrow(Y)-lag+1),,drop=FALSE] #drop the last h-1 rows so that the window ends at the forecast origin (nothing is dropped when lag=1)
  
  Y2=cbind(Y[,1]) #Y variable (x1)
  aux=embed(Y2,max_p+lag) #create max_p lags + forecast horizon shift (=lag option)
  y=aux[,1] #  Y variable aligned/adjusted for missing data due do lags
  X=aux[,-c(1:(ncol(Y2)*lag))] # lags of Y (predictors) corresponding to forecast horizon   
  
  X.out=tail(aux,1)[1:ncol(X)] #retrieve the last max_p observations of the window - these are the regressors for the h-step forecast (same construction for every horizon)
  
  candidates=if(is.null(p_fixed)) seq(1,ncol(X),1) else p_fixed
  bb=Inf #initialize the "best BIC" at a huge number
  for(i in candidates){ 
    m=lm(y~X[,1:i]) #estimate AR(i) by OLS
    crit=BIC(m) #retrieve BIC
    if(crit<bb){ #if BIC improved
      bb=crit #updated the new best BIC
      model=m #save the model object
      ar.coef=coef(model) #save coefficients
      p=i #save the selected lag order
    }
  }
  coef=rep(0,ncol(X)+1) #blank vector for coefficients
  coef[1:length(ar.coef)]=ar.coef #fill in the model coefficients (zeros may remain if model less than order max_p)
  
  pred=c(1,X.out)%*%coef #make a forecast using the last few observations: a direct h-step forecast
  
  return(list("model"=model,"pred"=pred,"coef"=coef,"p"=p)) #save estimated AR regression, prediction, estimated coefficients and selected order
}


#This function will repeatedly call the previous function in the rolling window h-step forecasting

#Inputs for the function:

#1) Y: one-column matrix with the x1 series

#2) nprev_covid - number of out-of-sample observations (at the end of the sample): 352 = Sep 1992 - Dec 2021

#3) nprev_nocovid - number of those that are pre-COVID (the first 328 = Sep 1992 - Dec 2019)

#4) lag - the forecast horizon

#5) max_p - the largest lag order BIC may choose


ar.rolling.window=function(Y,nprev_covid,nprev_nocovid,lag=1,max_p=12,p_fixed=NULL){
  
  save.coef=matrix(NA,nprev_covid,max_p+1) #blank matrix for coefficients at each iteration
  save.pred=matrix(NA,nprev_covid,1) #blank for forecasts
  save.p=numeric(nprev_covid) #blank for selected lag orders
  for(i in nprev_covid:1){  #NB: backwards FOR loop: going from 352 down to 1
    
    Y.window=Y[(1+nprev_covid-i):(nrow(Y)-i),,drop=FALSE] #define the estimation window (first one: 1 to nrow-nprev, then 2 to nrow-nprev+1 etc.)
    fact=runAR(Y.window,lag,max_p,p_fixed) #call the function to fit the AR(p) selected on BIC and generate h-step forecast
    save.coef[(1+nprev_covid-i),]=fact$coef #save estimated coefficients
    save.pred[(1+nprev_covid-i),]=fact$pred #save the forecast
    save.p[1+nprev_covid-i]=fact$p #save the selected lag order
  }
  
  #Some helpful stuff:
  real=Y[,1] #get actual values
  plot(real,type="l")
  lines(c(rep(NA,length(real)-nprev_covid),save.pred),col="red") #padded with NA for blanks, plot predictions vs. actual
  
  actual=tail(real,nprev_covid) #realized values over the whole test window
  
  mse=mean((actual-save.pred)^2) #MSE, COVID-inclusive (all nprev_covid forecasts)
  rmse=sqrt(mse) #RMSE, COVID-inclusive
  mae=mean(abs(actual-save.pred)) #MAE, COVID-inclusive
  mse.pre=mean((actual[1:nprev_nocovid]-save.pred[1:nprev_nocovid])^2) #MSE, pre-COVID (first nprev_nocovid forecasts)
  rmse.pre=sqrt(mse.pre) #RMSE, pre-COVID
  mae.pre=mean(abs(actual[1:nprev_nocovid]-save.pred[1:nprev_nocovid])) #MAE, pre-COVID
  errors=c("mse"=mse,"rmse"=rmse,"mae"=mae,"mse.pre"=mse.pre,"rmse.pre"=rmse.pre,"mae.pre"=mae.pre) #stack errors in a vector
  
  return(list("pred"=save.pred,"coef"=save.coef,"p"=save.p,"errors"=errors)) #return forecasts, history of estimated coefficients, selected orders, and errors for the period.
}



# Sample split (by forecast origin, same as the tuned models)
# Training set: 120 + 100 = 220 obs, since the AR has nothing to tune.
# Test set: 
#   328 obs (Sep 1992 - Dec 2019) for pre-covid; 
#   352 obs (Sep 1992 - Dec 2021) for COVID-inclusive


nprev_nocovid <- 328  # pre-COVID test window
nprev_covid <- 352   # COVID-inclusive test window



select_p=function(Y,lag,nprev,max_p=12){
  ntrain=nrow(Y)-nprev                        # first 220 obs, before the test window
  runAR(Y[1:ntrain,,drop=FALSE],lag,max_p)$p  # BIC on that window only
}

run_h=function(df_h,h){
  Y=make_series(df_h,h)
  p=select_p(Y,h,nprev_covid)
  res=ar.rolling.window(Y,nprev_covid,nprev_nocovid,lag=h,p_fixed=p)
  res$p_selected=p
  res
}


ar_h1_bic  <- run_h(df_h1,1)
ar_h3_bic  <- run_h(df_h3,3)
ar_h6_bic  <- run_h(df_h6,6)
ar_h12_bic <- run_h(df_h12,12)

c(h1=ar_h1_bic$p_selected, h3=ar_h3_bic$p_selected,
  h6=ar_h6_bic$p_selected, h12=ar_h12_bic$p_selected)   # selected p per horizon



err_stats=function(res,df_h,h){
  actual=tail(make_series(df_h,h)[,1],nprev_covid)  #the 352 out-of-sample actuals
  e=as.numeric(actual-res$pred)
  
  idx=list(pre=1:nprev_nocovid,                     #328 forecasts, Sep 1992 - Dec 2019
           covid_incl=1:nprev_covid)                #352 forecasts, Sep 1992 - Dec 2021
  
  map_dfr(names(idx),function(nm){
    ee=e[idx[[nm]]]
    tibble(horizon=h,period=nm,n=length(ee),
           mse=mean(ee^2),rmse=sqrt(mean(ee^2)),mae=mean(abs(ee)))
  })
}

oos_tbl=bind_rows(
  err_stats(ar_h1_bic, df_h1, 1),
  err_stats(ar_h3_bic, df_h3, 3),
  err_stats(ar_h6_bic, df_h6, 6),
  err_stats(ar_h12_bic,df_h12,12)
)

oos_tbl