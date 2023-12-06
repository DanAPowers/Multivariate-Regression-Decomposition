# mvdecompExample.R

setwd('~/documents/R_progs/mvdecomp')
source('~/documents/R_progs/mvdecomp/svydecomp_functions.R')

# libraries
require(survey)    # all models
require(sjstats)   # negative binomial

#  raw data
nhanes <- read.csv(file = 'nhanes2.csv')
#  set svy design
#  full design
nhsvy <- svydesign(ids=~psuid, weights=~leadwt, strata=~stratid, nest=TRUE, data=nhanes)
# no design (ignore warnings)
nhsvy <- svydesign(ids=~1, nest=TRUE, data=nhanes)
# note: robust std. errors are estimated for all models
# group difference decomposition: define groups
# define groups: black/non black
Asub <- subset(nhsvy, black == 1)
Bsub <- subset(nhsvy, black == 0)
################################################################################
# call: decomp.model(formula, Asub, Bsub, scale=1, printit=FALSE, reverse=FALSE)
# Asub and Bsub are lists that act as data frames for the models
################################################################################
#
# poisson/negative binomial models
#

m1a <- decomp.poisson(lead ~ age + rural + female + offset(logBMI), Asub, Bsub, 
                    scale=1, printit=TRUE, 
                    reverse=FALSE)

m1b <- decomp.poisson(lead ~ age + rural + female, Asub, Bsub, 
                    scale=1, printit=TRUE, 
                    reverse=FALSE)

m2a <- decomp.negbin(lead ~ age + rural + female + offset(logBMI), Asub, Bsub, 
                   scale=1, printit=TRUE, 
                   reverse=FALSE)

m2b <- decomp.negbin(lead ~ age + rural + female, Asub, Bsub, 
                   scale=1, printit=TRUE, 
                   reverse=FALSE)
#
# data for binary and linear models
#
nhanes <- read.csv(file='nhanes1.csv')
# set svy design
nhsvy <- svydesign(ids=~psuid, weights=~finalwgt, strata=~stratid, nest=TRUE, data=nhanes)
# no design (ignore warnings)
nhsvy <- svydesign(ids=~1, nest=TRUE, data=nhanes)
# define groups (female/male)
Asub <- subset(nhsvy, female == 1)
Bsub <- subset(nhsvy, female == 0)


# logit
#
m3a <- decomp.logit(highbp ~ age + BMI + rural + black, Asub, Bsub, 
                   scale=10, printit=TRUE)

# define different groups
Asub <- subset(nhsvy, black == 1)
Bsub <- subset(nhsvy, black == 0)

m3b <- decomp.logit(highbp ~ age + BMI + rural + female, Asub, Bsub, 
                   scale=10, printit=TRUE)

# probit
#
m4 <- decomp.probit(highbp ~ age + BMI + rural + female, Asub, Bsub, 
                   scale=10, printit=TRUE)

# complementary log log
#
m5 <- decomp.cloglog(highbp ~ age + BMI + rural + female, Asub, Bsub, 
                    scale=10, printit=TRUE)

#
# linear
# make subsets etc. (female/male)
Asub <- subset(nhsvy, female == 1)
Bsub <- subset(nhsvy, female == 0)

m6 <- decomp.linear(tcresult ~ age + rural + black, Asub, Bsub, printit=TRUE)


# averaging two decompistions: 
require(tidyverse)

# average and graph
m1 <- decomp.linear(tcresult ~ age + rural + black, Asub, Bsub)
m2 <- decomp.linear(tcresult ~ age + rural + black, Asub, Bsub, reverse=TRUE)

E.ave <- ave.decompE(m1,m2,1)
C.ave <- ave.decompC(m1,m2,1)
E.tot <- tot.decompE(m1,m2,1)
C.tot <- tot.decompC(m1,m2,1)
# simple tables
list(Detailed=E.ave[, -c(1)],Total=E.tot)
list(Detailed=C.ave[, -c(1)],Total=C.tot)

require(tidyverse)
# coef plots
pE <- ggplot(E.ave, aes(b, term)) +
  geom_point(color="seagreen") +
  geom_errorbarh(aes(xmin=b.lower, xmax=b.upper), 
                 linewidth=.3,
                 height=.3,
                 linetype=1,
                 color="turquoise4") +
  geom_vline(xintercept = 0, lty = 2, color="green") +
  theme_gray() +
  labs(
    x = "E-coefs",
    y = NULL,
    title = "Plot of Characteristic Effects (averaged)"
  ) 

pC <- ggplot(C.ave, aes(b, term)) +
  geom_point(color="seagreen") +
  geom_errorbarh(aes(xmin=b.lower, xmax=b.upper), 
                 linewidth=.3,
                 height=.3,
                 linetype=1,
                 color="turquoise4") +
  geom_vline(xintercept = 0, lty = 2, color="green") +
  theme_gray() +
  labs(
    x = "C-coefs",
    y = NULL,
    title = "Plot of Coefficient Effects (averaged)"
  ) 

gridExtra::grid.arrange(pE,pC, nrow=1, ncol=2)


# source this first

ave.decompE <- function(D1, D2, scale) {
  s  <- scale
  vn <- D1$vnames
  lv <- length(vn)
  b1 <- D1$bE[2:lv]
  b2 <- D2$bE[2:lv]
  v1 <- D1$varbE[2:lv,2:lv]
  v2 <- D2$varbE[2:lv,2:lv]
  b <- s * (b1 - b2)/2
  v <- s^2 * (v1 + v2)/4
  se.b <- sqrt(diag(v))
  #colnames(b, vn)
  return(data.frame(term=vn[2:lv], b=b, se=se.b, z=b/se.b, 
                    b.lower = b-1.96*se.b,
                    b.upper = b+1.96*se.b))
}

ave.decompC <- function(D1, D2, scale) {
  s  <- scale
  vn <- D1$vnames
  b1 <- D1$bC
  b2 <- D2$bC
  v1 <- D1$varbC
  v2 <- D2$varbC
  b <- s * (b1 - b2)/2
  v <- s^2 * (v1 + v2)/4
  se.b <- sqrt(diag(v))
  #colnames(b,vn)
  return(data.frame(term = vn, b=b, se=se.b, z=b/se.b, 
                    b.lower = b-1.96*se.b,
                    b.upper = b+1.96*se.b))
}


tot.decompE <- function(D1, D2, scale) {
  vn <- "E"
  s  <- scale
  E1 <- D1$E
  E2 <- D2$E
  v1 <- D1$varE
  v2 <- D2$varE
  E <- s * (E1 - E2)/2
  V <- (v1^2 + v2^2)/4
  seE  <- sqrt(V)*s
  return(data.frame(term=vn, b=E,se=seE, z=E/seE, 
                    b.lower=E - 1.96*seE,
                    b.upper=E + 1.96*seE))
}

tot.decompC <- function(D1,D2,scale) {
  s  <- scale
  vn <- "C"
  C1 <- D1$C
  C2 <- D2$C
  v1 <- D1$varC
  v2 <- D2$varC
   C <- s * (C1 - C2)/2
   V <- (v1^2 + v2^2)/4
   c
seC  <- sqrt(V)*s
return(data.frame(term=vn, b=C,se=seC, z=C/seC,
                  b.lower=C - 1.96*seC,
                  b.upper=C + 1.96*seC))
}


