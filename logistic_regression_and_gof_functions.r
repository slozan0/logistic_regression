######
# makes the sucess array, that is if there are 5 dead out 10 at 0.1 mg it creates
# 1 1 1 1 1 0 0 0 0 0 
# 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
# lozano.saul@gmail.com
# parameter "df" is a dataframe object  
######
MakeBinomialPairs <- function(df){
  retX <- list()
  retY <- list()
  for(i in 1:nrow(df)){
    x <- rep(log(df$dose[i]), df$batch[i])
    y <- c(rep(1, df$dead[i]), rep(0, df$batch[i]-df$dead[i]))
    retX <- append(retX, x)
    retY <- append(retY, y)
    retX <- unlist(retX, recursive=TRUE)
    retY <- unlist(retY, recursive=TRUE)
  }  
  return(myDf <- data.frame(logDose=retX, dead=retY))
}

####
# A function to do the Hosmer-Lemeshow test in R.
# Original by www.r-bloggers.com/author/denishaine/
# modified by lozano.saul@gmail.com
# param y are the observed and yhat the adjusted by the model and g?
####
Hosmerlem <- function (y, yhat, g) 
{
  try
  #todo what to do if cut is not unique maybe add a try loop.
  cutyhat <- tryCatch({
    cut(yhat, breaks = quantile(yhat, probs = seq(0, 1, 1/g)), include.lowest = T)
    },
    error=function(cond) {
      message("(╯°□°）╯ ︵ ┻━┻ The model has many ties in its predicted probabilities")
      message("(ಠ_ಠ) Hosmer-Lemeshow test is not applicable in this case")
      message("(ಠ_ಠ) Too few covariate values?")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      stop()
    },
  warning=function(cond) {
    message("Here's the original warning message:")
    message(cond)
    # Choose a return value in case of warning
  },
  finally={
  }
  ) 
  
  obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
  expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq <- sum((obs - expect)^2/expect)
  P <- 1 - pchisq(chisq, g - 2)
  c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
}
######
#reverse prediction of x given a logistic regression 
# lozano.saul@gmail.com
# parameter "y" is the observed value, intercept is the beta 0 and the sloep is beta 1 
###### 
PredictX <- function(y, intercept, slope){
  x <- (log(y/(1-y)) - intercept) / slope
  return(x)  
}

#####################################################################
# Summary: calculate the CI for beta1 "slope/dose" parameter using the delta method
# based on code by Dason Kurkiewicz http://dasonk.com/r/2013/02/09/Using-the-delta-method/
# modified by lozano.saul@gmail.com Nov 25th, 2014
# parameters 
################################################
ConfInt <- function(estimatedDose, p, glmModel){
  standerr <- deltamethod(~ (log(p/(1-p)) - x1)/x2, coef(glmModel), vcov(glmModel))
  ci <- estimatedDose + c(-1, 1) * 1.96 * standerr
  return(ci)
}

######
# Checks if a library is present. If present it loads. 
# If the library is not present tries to install
# and load the requested library.
# lozano.saul@gmail.com 25th Nov, 2014
# parameter(s) "libname" is the name of the library to check in string
######
Check_n_Load <- function(libName){
  if(require(libName, character.only=TRUE)){
    print(paste("@(^_^)@ ", libName," is loaded correctly"))
  } 
  else {
    print(paste("@(◔ ◡ ◔)@, trying to install", libName ))
    install.packages(libName)
    if(require(libName)){
      print(paste(libName, "@(^_^)@ installed and loaded"))
    } else {
      stop(paste("(╯°□°）╯ ︵ ┻━┻ could not install", libname))
    }
  }
}

Check_n_Load("msm")# needed for deltamethod function
Check_n_Load("rms")# needed le Cessie – van Houwelingen – Copas – Hosmer unweighted sum of squares test 
                   # for global goodness of fit 