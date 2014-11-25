# code for the statistical analysis of dose respose bottle assay
# lozano.saul@gmail.com 25th nov, 2014

source("logistic_regression_and_gof_functions.r")

bioAssay <- read.csv(file="bioassayIrma.txt", header=TRUE)
dfBioAssay <- MakeBinomialPairs(bioAssay)

glmModel <- glm(dead~logDose, data=dfBioAssay, family=binomial)
summary(glmModel)

intercept=glmModel$coef[[1]]
slope=glmModel$coef[[2]]

ld50 <- PredictX(.5,intercept, slope)
ld90 <- PredictX(.9,intercept, slope)
ld99 <- PredictX(.99,intercept, slope)

# Doing the Hosmer-Lemeshow test
Hosmerlem(dfBioAssay$dead, fitted(glmModel),9)

#le Cessie – van Houwelingen – Copas – Hosmer unweighted sum of squares test 
#for global goodness of fit 
lrmModel <- lrm(dead ~ logDose, data=dfBioAssay, linear.predictors=TRUE,
         method="lrm.fit", x=TRUE, y=TRUE, model = TRUE)
residuals(lrmModel, type = "gof")


# lets do Confidence intervals for ld50
p <- 0.5
ci <- ConfInt(ld50, p, glmModel)
c(estimatedx = ld50, ci=ci)

# lets do Confidence intervals for ld90
p <- 0.9
ci <- ConfInt(ld90, p, glmModel)
c(estimatedx = ld90, ci=ci)

# lets do Confidence intervals for ld99
p <- 0.99
ci <- ConfInt(ld99, p, glmModel)
c(estimatedx = ld99, ci=ci)