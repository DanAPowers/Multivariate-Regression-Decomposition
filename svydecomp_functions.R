#########################################################################
# Author: Dan Powers (11/22/23)
# functions for multivariate regression decomposition
# of difference in first moments--E[F(XB)]--between group a and group b: 
#
#  D = E[F(XaBa)] - E[F(XbBb)] 
#    = {E[F(XaBa)]  - E[F(XbBa)]} + {E[F(XbBa)] - E[F(XbBb)]|
# 
#  X <- characteristics (observed)
#  B <- effects         (estimated)
# 
# use svyglm for linear and nonlinear regression models
# compare to Stata -mvdcmp-
########################################################################


##################################
# model_setup_svy.R
##################################
require(survey)
###############
mod.out <- function(m, sub) {
  b   <- coef(m)
  x   <- m$x
  v   <- m$cov.unscaled
  sub <- sub
  w   <- weights(sub)
  mget(ls())
} 

mod.outp <- function(m,sub) {
  b   <- coef(m)
  x   <- m$x
  v   <- m$cov.unscaled
  off <- m$offset
  if (is.null(off)) off <- rep(0,nrow(x))
  sub <- sub
  mget(ls())
} 

mod.outnb <- function(m,m2,sub) {
  lb  <- length(coef(m)) 
  b   <- coef(m)[2:lb]
  x   <- m2$x
  v   <- m$sandwich[2:lb, 2:lb]
  off <- m2$offset
  if (is.null(off)) off <- rep(0,nrow(x))
  sub <- sub
  mget(ls())
} 

decomp.logit  <- function(form, designA, designB, scale=1, printit=FALSE, reverse=FALSE) {
  A <- substitute(svyglm(terms(form), family=quasibinomial(link=logit), design = designA,                
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  B <- substitute(svyglm(terms(form), family=quasibinomial(link=logit), design = designB,
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  if (reverse) outA <- mod.out(eval(B),Bsub)
  else outA <- mod.out(eval(A),Asub)
  if (reverse) outB <- mod.out(eval(A),Asub)
  else outB <- mod.out(eval(B), Bsub)
  Alist <- list(b=outA$b,x=outA$x,v=outA$v, sub=outA$sub)
  Blist <- list(b=outB$b,x=outB$x,v=outB$v, sub=outB$sub)
  decomp_logit(Alist, Blist, scale, printit)
}



decomp.probit  <- function(form, designA, designB, scale=NULL, printit=FALSE, reverse=FALSE) {
  A <- substitute(svyglm(terms(form), family=quasibinomial(link=probit), design = designA,                
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  B <- substitute(svyglm(terms(form), family=quasibinomial(link=probit), design = designB,
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  if (reverse) outA <- mod.out(eval(B),Bsub)
  else outA <- mod.out(eval(A),Asub)
  if (reverse) outB <- mod.out(eval(A),Asub)
  else outB <- mod.out(eval(B), Bsub)
  Alist <- list(b=outA$b,x=outA$x,v=outA$v, sub=outA$sub)
  Blist <- list(b=outB$b,x=outB$x,v=outB$v, sub=outB$sub)
  decomp_probit(Alist, Blist, scale, printit)
}


decomp.cloglog  <- function(form, designA, designB, scale=NULL, printit=FALSE, reverse=FALSE) {
  A <- substitute(svyglm(terms(form), family=quasibinomial(link=cloglog), design = designA,                
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  B <- substitute(svyglm(terms(form), family=quasibinomial(link=cloglog), design = designB,
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  if (reverse) outA <- mod.out(eval(B),Bsub)
  else outA <- mod.out(eval(A),Asub)
  if (reverse) outB <- mod.out(eval(A),Asub)
  else outB <- mod.out(eval(B), Bsub)
  Alist <- list(b=outA$b,x=outA$x,v=outA$v, sub=outA$sub)
  Blist <- list(b=outB$b,x=outB$x,v=outB$v, sub=outB$sub)
  decomp_cloglog(Alist, Blist, scale, printit)
}


decomp.linear  <- function(form, designA, designB, scale=NULL, printit=FALSE, reverse=FALSE) {
  A <- substitute(svyglm(terms(form), family=gaussian, design = designA,                
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  B <- substitute(svyglm(terms(form), family=gaussian, design = designB,
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  if (reverse) outA <- mod.out(eval(B),Bsub)
  else outA <- mod.out(eval(A),Asub)
  if (reverse) outB <- mod.out(eval(A),Asub)
  else outB <- mod.out(eval(B), Bsub)
  Alist <- list(b=outA$b,x=outA$x,v=outA$v, w = outA$w, sub=outA$sub)
  Blist <- list(b=outB$b,x=outB$x,v=outB$v, w = outB$w, sub=outB$sub)
  decomp_linear(Alist, Blist, scale, printit)
}

decomp.poisson  <- function(form, designA, designB, scale=NULL, printit=FALSE, reverse=FALSE) {
  A <- substitute(svyglm(terms(form), family=quasipoisson(), design = designA, 
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  B <- substitute(svyglm(terms(form), family=quasipoisson(), design = designB,
                         epsilon = 1e-9, maxit = 100, x=TRUE, rescale=TRUE))
  # invoke out.modp
  if (reverse) outA <- mod.outp(eval(B),Bsub)
  else outA <- mod.outp(eval(A),Asub)
  if (reverse) outB <- mod.outp(eval(A),Asub)
  else outB <- mod.outp(eval(B), Bsub)
  Alist <- list(b=outA$b,x=outA$x,v=outA$v, sub=outA$sub, off=outA$off)
  Blist <- list(b=outB$b,x=outB$x,v=outB$v, sub=outB$sub, off=outB$off)
  
  decomp_poisson(Alist, Blist, scale, printit)
}

decomp.negbin  <- function(form, designA, designB, scale=NULL, printit=FALSE, reverse=FALSE) {
  A  <- substitute(svyglm.nb(form, design = designA))
  A1 <- substitute(svyglm(terms(form), family=quasipoisson(), design = designA, x=TRUE))              
  B  <- substitute(svyglm.nb(form, design = designB))
  B1 <- substitute(svyglm(terms(form), family=quasipoisson(), design = designB, x=TRUE))   
  
  # invoke out.modp
  if (reverse) outA <- mod.outp(eval(B),eval(B1), Bsub)
  else outA <- mod.outnb(eval(A),eval(A1), Asub)
  if (reverse) outB <- mod.outp(eval(A), eval(A1), Asub)
  else outB <- mod.outnb(eval(B), eval(B1), Bsub)
  Alist <- list(b=outA$b,x=outA$x,v=outA$v, sub=outA$sub, off=outA$off)
  Blist <- list(b=outB$b,x=outB$x,v=outB$v, sub=outB$sub, off=outB$off)
  
  decomp_poisson(Alist, Blist, scale, printit)
}
##################################
# svy logistic-specific functions
##################################

decomp_logit <- function(A,B, scale=NULL, printit)  {
  
  if (is.null(scale)) scale <- 1  
  
  # translate args
  Asub  <- A$sub
  Bsub  <- B$sub
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  varbA <- A$v
  varbB <- B$v
  
  # get svy weighted x-means here
  
  m     <- as.data.frame(svymean(~xA, Asub, na=TRUE))
  mA    <- m$mean
  m     <- as.data.frame(svymean(~xB, Bsub, na=TRUE))
  mB    <- m$mean
  
  # logistic CDF -- F(x)
  F <- function(b,x) {
    xb <- x%*%b
    F  <- exp(xb)/(1 + exp(xb)) 
    return(F)
  }
  
  # logistic pdf -- dF(x)
  f <- function(b,x) {
    xb <- x%*%b
    f <-  exp(xb)/(1 + exp(xb))^2
    return(f)
  }
  
  decomp_svybinary(bA, bB, xA, xB, varbA, varbB, 
                    mA, mB, Asub, Bsub, printit, scale, F, f)
}

################################
# svy probit-specific functions
################################

decomp_probit <- function(A,B,scale=NULL, printit)  {
  if (is.null(scale)) scale <- 1  
  
  # translate args
  Asub  <- A$sub
  Bsub  <- B$sub
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  varbA <- A$v
  varbB <- B$v
  
  # get svy weighted x-means here
  
  m     <- as.data.frame(svymean(~xA, Asub, na=TRUE))
  mA    <- m$mean
  m     <- as.data.frame(svymean(~xB, Bsub, na=TRUE))
  mB    <- m$mean
  
  #  normal CDF -- F(x)
  F <- function(b,x) {
    xb <- x%*%b
    F  <- pnorm(xb) 
    return(F)
  }
  
  # normal pdf -- dF(x)
  f <- function(b,x) {
    xb <- x%*%b
    f <-  dnorm(xb)
    return(f)
  }
  
  decomp_svybinary(bA, bB, xA, xB, varbA, varbB, 
                   mA, mB, Asub, Bsub, printit, scale, 
                   F, f)
}

################################
# svy cloglog-specific functions
################################

decomp_cloglog <- function(A,B, scale=NULL, printit)  {
  
  if (is.null(scale)) scale <- 1  
  
  # translate args
  
  Asub  <- A$sub
  Bsub  <- B$sub
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  varbA <- A$v
  varbB <- B$v
  
  # get svy weighted x-means here
  
  m     <- as.data.frame(svymean(~xA, Asub, na=TRUE))
  mA    <- m$mean
  m     <- as.data.frame(svymean(~xB, Bsub, na=TRUE))
  mB    <- m$mean
  
  # extreme-value CDF -- F(x)
  F <- function(b,x) {
    xb <- x%*%b
    F <- 1 - exp(-exp(xb)) 
    return(F)
  }
  # extreme-value pdf -- dF(x)
  f <- function(b,x) {
    xb <- x%*%b
    f <- exp(-exp(xb))
    return(f)
  }
  
  decomp_svybinary(bA, bB, xA, xB, varbA, varbB, 
                   mA, mB, Asub, Bsub, printit, scale, 
                   F, f)
}

###############################
# functions for svy count model
###############################
########################
# model_setup_svycount.R
########################

decomp_poisson <- function(A, B, scale=NULL, printit=FALSE)  {
  
  if (is.null(scale)) scale <- 1 
  
  # translate args
  
  Asub  <- A$sub
  Bsub  <- B$sub
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  varbA <- A$v
  varbB <- B$v
  offA  <- A$off
  offB  <- B$off
  xpA   <- exp(offA) # define exposure for for denominators of rates.
  xpB   <- exp(offB)
  mxpA  <- svymean(~xpA, Asub, na=TRUE) # apply svy weighting
  mxpB  <- svymean(~xpB, Bsub, na=TRUE)
  # get svy weighted x-means here
  m     <- as.data.frame(svymean(~xA, Asub, na=TRUE))
  mA    <- m$mean
  m     <- as.data.frame(svymean(~xB, Bsub, na=TRUE))
  mB    <- m$mean
  
  #############################
  # poisson-specific functions
  #############################
  
  # E(x)
  F <- function(b,x,off) {
    xb <- x%*%b
    F <- exp(xb + off)  
    return(F)
  }
  # dE(x)
  f <- function(b,x,off) {
    xb <- x%*%b
    f <- exp(xb + off) 
    return(f)
  }
  
  decomp_svycount(bA, bB, xA, xB, varbA, varbB, offA, offB, mxpA, mxpB,
                  mA, mB, Asub, Bsub, printit, scale, F, f)
}


#################################
# functions for svy linear model
#################################

decomp_linear <- function(A,B,scale=NULL, printit)  {
  
  if (is.null(scale)) scale <- 1  
  
  # translate args
  Asub  <- A$sub
  Bsub  <- B$sub
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  varbA <- A$v
  varbB <- B$v
  wA    <- weights(Asub)
  wB    <- weights(Bsub)
  
  # get svy weighted x-means here
  
  m     <- as.data.frame(svymean(~xA, Asub, na=TRUE))
  mA    <- m$mean
  m     <- as.data.frame(svymean(~xB, Bsub, na=TRUE))
  mB    <- m$mean
  
  #  linear model -- E(x)
  F <- function(b,x,w) {
    xb <- x%*%b
    F  <- xb * w
    return(F)
  }
  
  #  linear model -- dE(x)
  
  f <- function(w) {
    f  <- w
    return(f)
  }
  
  decomp_svylinear(bA, bB, xA, xB, varbA, varbB, 
                   mA, mB, Asub, Bsub, wA, wB, printit, scale, 
                   F, f)
}


###########################################################
#                      MAIN                               #
###########################################################
# functions specific to complex survey binary models
###########################################################
decomp_svybinary <- function(bA, bB, xA, xB, varbA, varbB, 
                   mA, mB, Asub, Bsub, printit, scale, F, f) {
  
Wdx.F  <- function(b,x1,x2) {
  # Yun weight function (E composition) 
  # b      = coef
  # x1, x2 = weighted means
  #
  A <- (x1-x2)%*%b
  Wdx <- NULL
  for (i in 1:length(b)){
    Wdx[i] <- (x1[i] - x2[i])*b[i] / A
  }
  return(Wdx)
}

Wdb.F  <- function(b1,b2,x) {
  # Yun weight function (C coefficients)
  # b1, b2 = coef
  # x      = weighted mean
  #
  A <- x%*%(b1-b2)
  Wdb <- NULL
  for (i in 1:length(x)){
    Wdb[i] <- (x[i]*(b1[i] - b2[i])) / A
  }
  return(Wdb)
}

dW.F <- function(b,x1,x2) {
  # derivative of Wdx  (K x K result) 
  # b      = coef
  # x1, x2 = weighted mean
  #
  dW <- array(NA,c(length(b), length(b)))
  A <- (x1-x2)%*%b   
  for (k in 1:length(b)) {
    for (l in 1:length(b)) {
      dW[k,l]  <-  as.numeric(k==l) * ((x1[k] - x2[k])/A) - 
        (b[k] * ( (x1[k] - x2[k])*(x1[l] - x2[l]) )/A^2 )
    }
  }
  return(dW)
}     

dwA.F <- function(b1,b2,x2) {
  # derivative of Wdb part A --  K x K  
  # b1, b2 = coef
  # x2      = weighted mean
  #
  dwA1 <- array(NA, c(length(b1), length(b2)))
  A <- x2%*%(b1-b2)
  for (k in 1:length(b1)){
    for (l in 1:length(b2)) {
      dwA1[k,l] <- as.numeric(k==l) * (x2[k]/A) - (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 
    } 
  }
  return(dwA1)
}

dwB.F <- function(b1,b2,x2) { 
  # derivative of Wdb part B
  # b1, b2 = coef
  # x2      = weighted mean
  #  
  dwB1 <- array(NA, c(length(b1), length(b2)))
  A <- x2%*%(b1-b2)
  for (k in 1:length(b1)){
    for (l in 1:length(b2)) {
      dwB1[k,l] <- (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 - as.numeric(k==l) * (x2[k]/A) 
    } 
  }
  return(dwB1)
}

wbA  <-  dwA.F(bA,bB,mB)
wbB  <-  dwB.F(bA,bB,mB)
#############
# Yun weights
#############
Wdx  <-  Wdx.F(bA, mA, mB)
Wdb  <-  Wdb.F(bA, bB, mB)
############################################
# Convention: 1st moment A - 1st moment B
# group A is treated as reference (standard)
# group B is treated as comparison
############################################
# decomp total: E & C
E   <- svymean(~F(bA, xA), Asub) -  svymean(~F(bA, xB), Bsub)
# Coefficients + Unexplained: 
C   <- svymean(~F(bA, xB), Bsub) -  svymean(~F(bB, xB), Bsub)

# get dEdb for variance estimator
dWx  <- dW.F(bA, mA, mB)
dEdb <- array(NA, c(length(bA), length(bA)))
fxA  <- NULL
fxB  <- NULL
for (k in 1:length(bA)) {
  for (l in 1:length(bA)) {
    # 
    df <- svymean(~as.vector(f(bA,xA)*xA[,l]), Asub) - 
          svymean(~as.vector(f(bA,xB)*xB[,l]), Bsub)
    dEdb[k,l] <- Wdx[k] * df + dWx[k,l] * E
  }
}

##   Variances
#    Composition (E)
Var.E.k <- dEdb%*%varbA%*%t(dEdb)
seWdx   <- sqrt(diag(Var.E.k))
dCdbA   <- array(NA, c(length(bA), length(bA)))
dCdbB   <- array(NA, c(length(bB), length(bB)))
for (k in 1:length(bA)){
  for (l in 1:length(bB)) { 
    dCdbA[k,l] <- Wdb[k] * svymean(~as.vector(f(bA,xB)*xB[,l]), Bsub) + wbA[k,l]*C
    dCdbB[k,l] <- wbB[k,l] * C - Wdb[k]*svymean(~as.vector(f(bB,xB)*xB[,l]),Bsub)
  }
} 

##  Variances
#   Coefficients (C)
Var.C.k <- dCdbA%*%varbA%*%t(dCdbA) + dCdbB%*%varbB%*%t(dCdbB)
seWdb   <- sqrt(diag(Var.C.k))

# get asymptotic variance of E_k and C_k
# E_k
dEdb.0 <- NULL
for (k in 1:length(bA)){
  dEdb.0[k] <- svymean(~as.vector(f(bA,xA)*xA[,k]), Asub) - 
               svymean(~as.vector(f(bA,xB)*xB[,k]), Bsub)
}

#C_k
dCdb.0A <- NULL
for (k in 1:length(bA)){
  dCdb.0A[k] <- svymean(~as.vector(f(bA,xB)*xB[,k]),Bsub)
}
dCdb.0B <- NULL
for (k in 1:length(bB)){
  dCdb.0B[k] <- svymean(~as.vector(f(bB,xB)*xB[,k]),Bsub)
}

#         1 x k   k x k     k x 1
var.E0 <- t(dEdb.0)%*%varbA%*%dEdb.0
s.E0   <- sqrt(var.E0)

var.C0 <- t(dCdb.0A)%*%varbA%*%dCdb.0A + t(dCdb.0B)%*%varbB%*%dCdb.0B 
s.C0   <- sqrt(var.C0)

R <- E + C
v.names <- names(bA)
if (is.null(v.names)) {
  v.names <- paste("Var", 1:length(bA), sep = "") # if names not available
}
# 
ymA <- svymean(~F(bA, xA), Asub) * scale
ymB <- svymean(~F(bB, xB), Bsub) * scale


if (printit==TRUE) {
  
  cat("","\n")
  cat("Binary Regression Decomposition","\n")
  cat("","\n")
  cat("Mean Group A           ", formatC(ymA, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Mean Group B           ", formatC(ymB, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("__________________________________________________________________", "\n")
  
  cat("Total Difference (E + C)", formatC((E+C)*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Raw Contrib......(E,  C)", formatC(E*scale, flag=" ", dig=4, wid=10, format="f"), 
      formatC(C*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Std. Errors......(E,  C)", formatC(s.E0*scale, flag=" ", dig=4, wid=10, format="f"),
      formatC(s.C0*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Z-values.........(E,  C)", formatC(E/s.E0, flag=" ", dig=4, wid=10, format="f"),
      formatC(C/s.C0, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Pct Contrib......(E,  C)", formatC(100*E/(E+C), flag=" ", dig=4, wid=10, format="f"),
      formatC(100*C/(E+C), flag=" ", dig=4, wid=10, format="f"), "\n" )
  cat("","\n")
  cat("SUMMARY MEASURES","\n")
  cat("        Comparison Group", "\n")
  cat("        ----------------", "\n")
  cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
  cat("__________________________________________________________________", "\n")
  for (i in 1:length(mA)){
    cat(formatC(v.names[i],format="s", wid=20),
        formatC(mA[i],flag=" ",dig=4,wid=10, format="f"),
        formatC(bA[i],wid=10,dig=4, format="f"),
        formatC(sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), 
        formatC(bA[i]/sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), "\n")
  }
  cat("","\n")
  cat("        Reference Group", "\n")
  cat("        ---------------", "\n")
  cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
  cat("__________________________________________________________________", "\n")
  for (i in 1:length(mB)){
    cat(formatC(v.names[i],format="s", wid=20),
        formatC(mB[i],flag=" ",   dig=4,wid=10, format="f"),
        formatC(bB[i],wid=10,     dig=4, format="f"),
        formatC(sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), 
        formatC(bB[i]/sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), "\n")
  }
  cat("","\n")
  cat("DUE TO DIFFERENCE IN COMPOSITION (E)", "\n")
  cat("","\n")
  cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "        PCT",  "\n")
  cat("___________________________________________________________________", "\n")
  
  idx <- grepl("(Intercept)", v.names)
  
  for (i in 1:length(bA)){
    if (!idx[i]) {
    cat(formatC(v.names[i],format="s", wid=20),
        formatC(E*Wdx[i]*scale,       flag=" ", dig=4, wid=10, format="f"),
        formatC(seWdx[i]*scale,       flag=" ", dig=4, wid=10, format="f"),
        formatC(E*Wdx[i]/seWdx[i],    flag=" ", dig=4, wid=10, format="f"),
        formatC(100*(E*Wdx[i]/(E+C)), flag="+", dig=1, wid=10, format="f"), "\n")
    }
  }
  cat("","\n")
  cat("DUE TO DIFFERENCE IN COEFFICIENTS (C)", "\n")
  cat("","\n")
  cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "         PCT"  ,"\n")
  cat("___________________________________________________________________", "\n")
  for (i in 1:length(bA)){
    cat(formatC(v.names[i],format="s", wid=20),
        formatC(C*Wdb[i]*scale,       flag=" ",  dig=4, wid=10, format="f"),
        formatC(seWdb[i]*scale,       flag=" ",  dig=4, wid=10, format="f"),
        formatC(C*Wdb[i]/seWdb[i],    flag=" ",  dig=4, wid=10, format="f"), 
        formatC(100*(C*Wdb[i]/(E+C)), flag="+",  dig=1, wid=10, format="f"),"\n")
    
  }
} # end printit

# collect things to return
K        <- length(bA)
bE       <- matrix(E*Wdx,K,1) 
varbE    <- matrix(Var.E.k,K,K)
bC       <- matrix(C*Wdb,K,1)
varbC    <- matrix(Var.C.k,K,K)
E        <- E[[1]]
varE     <- var.E0
C        <- C[[1]]
varC     <- var.C0
rownames(bE) <- v.names
rownames(bC) <- v.names
colnames(varbE) <- rownames(varbE) <- v.names
colnames(varbC) <- rownames(varbC) <- v.names

# returns to model object
return(list( vnames=v.names, 
             scale=scale,
             coefE=bE, 
             vcovE=varbE, 
             coefC=bC, 
             vcovC=varbC,
             E=E,
             C=C,    
             varE=varE,
             varC=varC )
       )


}  # end decomp_svybinary

#####################################################
#   functions specific to complex survey count models
#####################################################
#
# lineage: svydecomp_functions_count.R
#
decomp_svycount <- function(bA, bB, xA, xB, varbA, varbB, offA, offB, mxpA, mxpB,
                            mA, mB, Asub, Bsub, printit, scale, F, f) {
  
  Wdx.F  <- function(b,x1,x2) {
    # Yun weight function (E composition) 
    # b      = coef
    # x1, x2 = weighted means
    #
    A <- (x1-x2)%*%b
    Wdx <- NULL
    for (i in 1:length(b)){
      Wdx[i] <- (x1[i] - x2[i])*b[i] / A
    }
    return(Wdx)
  }
  
  Wdb.F  <- function(b1,b2,x) {
    # Yun weight function (C coefficients)
    # b1, b2 = coef
    # x      = weighted mean
    #
    A <- x%*%(b1-b2)
    Wdb <- NULL
    for (i in 1:length(x)){
      Wdb[i] <- (x[i]*(b1[i] - b2[i])) / A
    }
    return(Wdb)
  }
  
  dW.F <- function(b,x1,x2) {
    # derivative of Wdx  (K x K result) 
    # b      = coef
    # x1, x2 = weighted mean
    #
    dW <- array(NA,c(length(b), length(b)))
    A <- (x1-x2)%*%b   
    for (k in 1:length(b)) {
      for (l in 1:length(b)) {
        dW[k,l]  <-  as.numeric(k==l) * ((x1[k] - x2[k])/A) - 
          (b[k] * ( (x1[k] - x2[k])*(x1[l] - x2[l]) )/A^2 )
      }
    }
    return(dW)
  }     
  
  dwA.F <- function(b1,b2,x2) {
    # derivative of Wdb part A --  K x K  
    # b1, b2 = coef
    # x2      = weighted mean
    #
    dwA1 <- array(NA, c(length(b1), length(b2)))
    A <- x2%*%(b1-b2)
    for (k in 1:length(b1)){
      for (l in 1:length(b2)) {
        dwA1[k,l] <- as.numeric(k==l) * (x2[k]/A) - (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 
      } 
    }
    return(dwA1)
  }
  
  dwB.F <- function(b1,b2,x2) { 
    # derivative of Wdb part B
    # b1, b2 = coef
    # x2     = weighted mean
    #  
    dwB1 <- array(NA, c(length(b1), length(b2)))
    A <- x2%*%(b1-b2)
    for (k in 1:length(b1)){
      for (l in 1:length(b2)) {
        dwB1[k,l] <- (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 - as.numeric(k==l) * (x2[k]/A) 
      } 
    }
    #          cat(dwB1,"\n")
    return(dwB1)
  }
  
  wbA  <-  dwA.F(bA,bB,mB)
  wbB  <-  dwB.F(bA,bB,mB)
  
  #############
  # Yun weights
  #############
  Wdx  <-  Wdx.F(bA, mA, mB)
  Wdb  <-  Wdb.F(bA, bB, mB)
  ############################################
  # Convention: 1st moment A - 1st moment B
  # group A is treated as reference (standard)
  # group B is treated as comparison
  ############################################
  # decomp total: E & C
  #####################
  E   <- svymean(~F(bA, xA, offA), Asub)/mxpA -  
         svymean(~F(bA, xB, offB), Bsub)/mxpB
  # Coefficients + Unexplained: 
  C   <- svymean(~F(bA, xB, offB), Bsub)/mxpB -  
         svymean(~F(bB, xB, offB), Bsub)/mxpB
  
  # get dEdb for variance estimator
  dWx  <- dW.F(bA, mA, mB)
  dEdb <- array(NA, c(length(bA), length(bA)))
  for (k in 1:length(bA)) {
    for (l in 1:length(bA)) {
      # 
      df <- svymean(~as.vector(f(bA,xA,offA)*xA[,l]), Asub)/mxpA - 
            svymean(~as.vector(f(bA,xB,offB)*xB[,l]), Bsub)/mxpB
      dEdb[k,l] <- Wdx[k] * df + dWx[k,l] * E
    }
  }
  
  # pois reg coefficient variances (using varbA and varbB)
  ##   Variances
  #    Composition
  Var.E.k <- dEdb%*%varbA%*%t(dEdb)
  seWdx   <- sqrt(diag(Var.E.k))
  dCdbA   <- array(NA, c(length(bA), length(bA)))
  dCdbB   <- array(NA, c(length(bB), length(bB)))
  for (k in 1:length(bA)){
    for (l in 1:length(bB)) { 
      dCdbA[k,l] <- Wdb[k] * svymean(~as.vector(f(bA,xB,offB)*xB[,l]), Bsub)/mxpB + wbA[k,l]*C
      dCdbB[k,l] <- wbB[k,l] * C - Wdb[k]*svymean(~as.vector(f(bB,xB,offB)*xB[,l]),Bsub)/mxpB
    }
  } 
  
  ##  Variances
  #   Coefficients
  Var.C.k <- dCdbA%*%varbA%*%t(dCdbA) + dCdbB%*%varbB%*%t(dCdbB)
  seWdb <- sqrt(diag(Var.C.k))
  
  # get asymptotic variance of E_k and C_k
  # E_k
  dEdb.0 <- NULL
  for (k in 1:length(bA)){
    dEdb.0[k] <- svymean(~as.vector(f(bA,xA,offA)*xA[,k]), Asub)/mxpA - 
                 svymean(~as.vector(f(bA,xB,offB)*xB[,k]), Bsub)/mxpB
  }
  
  #C_k
  dCdb.0A <- NULL
  for (k in 1:length(bA)){
    dCdb.0A[k] <- svymean(~as.vector(f(bA,xB,offB)*xB[,k]),Bsub)/mxpB
  }
  dCdb.0B <- NULL
  for (k in 1:length(bB)){
    dCdb.0B[k] <- svymean(~as.vector(f(bB,xB,offB)*xB[,k]),Bsub)/mxpB
  }
  
  #         1 x k   k x k     k x 1
  var.E0 <- t(dEdb.0)%*%varbA%*%dEdb.0
  s.E0   <- sqrt(var.E0)
  
  var.C0 <- t(dCdb.0A)%*%varbA%*%dCdb.0A + t(dCdb.0B)%*%varbB%*%dCdb.0B 
  s.C0   <- sqrt(var.C0)
  
  R <- E + C
  v.names <- names(bA)
  if (is.null(v.names)) {
    v.names <- paste("Var", 1:length(bA), sep = "") # if names not available
  }
  
    ymA <- svymean(~F(bA, xA, offA), Asub)/mxpA * scale
    ymB <- svymean(~F(bB, xB, offB), Bsub)/mxpB * scale
    
    if (printit==TRUE) {
      cat("","\n")
      cat("Count Regression Decomposition","\n")
      cat("","\n")
      cat("Mean Group A           ", formatC(ymA, flag=" ", dig=4, wid=10, format="f"), "\n")
      cat("Mean Group B           ", formatC(ymB, flag=" ", dig=4, wid=10, format="f"), "\n")
      cat("__________________________________________________________________", "\n")
    cat("Total Difference (E + C)", formatC((E+C)*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Raw Contrib......(E,  C)", formatC(E*scale, flag=" ", dig=4, wid=10, format="f"), 
        formatC(C*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Std. Errors......(E,  C)", formatC(s.E0*scale, flag=" ", dig=4, wid=10, format="f"),
        formatC(s.C0*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Z-values.........(E,  C)", formatC(E/s.E0, flag=" ", dig=4, wid=10, format="f"),
        formatC(C/s.C0, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Pct Contrib......(E,  C)", formatC(100*E/(E+C), flag=" ", dig=4, wid=10, format="f"),
        formatC(100*C/(E+C), flag=" ", dig=4, wid=10, format="f"), "\n" )
    cat("","\n")
    cat("SUMMARY MEASURES","\n")
    cat("        Comparison Group", "\n")
    cat("        ----------------", "\n")
    cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
    cat("__________________________________________________________________", "\n")
    for (i in 1:length(mA)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(mA[i],flag=" ",dig=4,wid=10, format="f"),
          formatC(bA[i],wid=10,dig=4, format="f"),
          formatC(sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), 
          formatC(bA[i]/sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), "\n")
    }
    cat("","\n")
    cat("        Reference Group", "\n")
    cat("        ---------------", "\n")
    cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
    cat("__________________________________________________________________", "\n")
    for (i in 1:length(mB)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(mB[i],flag=" ",   dig=4,wid=10, format="f"),
          formatC(bB[i],wid=10,     dig=4, format="f"),
          formatC(sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), 
          formatC(bB[i]/sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), "\n")
    }
    cat("","\n")
    cat("DUE TO DIFFERENCE IN COMPOSITION (E)", "\n")
    cat("","\n")
    cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "        PCT",  "\n")
    cat("___________________________________________________________________", "\n")
    
    idx <- grepl("(Intercept)", v.names)
    
    for (i in 1:length(bA)){
      if (!idx[i]) {
        cat(formatC(v.names[i],format="s", wid=20),
            formatC(E*Wdx[i]*scale,       flag=" ", dig=4, wid=10, format="f"),
            formatC(seWdx[i]*scale,       flag=" ", dig=4, wid=10, format="f"),
            formatC(E*Wdx[i]/seWdx[i],    flag=" ", dig=4, wid=10, format="f"),
            formatC(100*(E*Wdx[i]/(E+C)), flag="+", dig=1, wid=10, format="f"), "\n")
      }
    }
    
    cat("","\n")
    cat("DUE TO DIFFERENCE IN COEFFICIENTS (C)", "\n")
    cat("","\n")
    cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "         PCT"  ,"\n")
    cat("___________________________________________________________________", "\n")
    for (i in 1:length(bA)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(C*Wdb[i]*scale,       flag=" ",  dig=4, wid=10, format="f"),
          formatC(seWdb[i]*scale,       flag=" ",  dig=4, wid=10, format="f"),
          formatC(C*Wdb[i]/seWdb[i],    flag=" ",  dig=4, wid=10, format="f"), 
          formatC(100*(C*Wdb[i]/(E+C)), flag="+",  dig=1, wid=10, format="f"),"\n")
      
    }
  } # end printit
  
  K        <- length(bA)
  # store these
  bE       <- matrix(E*Wdx,K,1) 
  varbE    <- matrix(Var.E.k,K,K)
  bC       <- matrix(C*Wdb,K,1)
  varbC    <- matrix(Var.C.k,K,K)
  E        <- E[[1]]
  varE     <- var.E0
  C        <- C[[1]]
  varC     <- var.C0
  rownames(bE) <- v.names
  rownames(bC) <- v.names
  colnames(varbE) <- rownames(varbE) <- v.names
  colnames(varbC) <- rownames(varbC) <- v.names
  
  # returns 
  return(list(vnames=v.names, E=E, varE=varE, C=C, varC=varC,
              bE=bE, varbE=varbE, 
              bC=bC, varbC=varbC))
  
}  # end decomp_svycount


###########################################################
# functions specific to complex survey linear models
###########################################################

decomp_svylinear <- function(bA, bB, xA, xB, varbA, varbB, 
                             mA, mB, Asub, Bsub, wA, wB, printit, scale, F, f) {
  

  Wdx.F  <- function(b,x1,x2) {
    # Yun weight function (E composition) 
    # b      = coef
    # x1, x2 = weighted means
    #
    A <- (x1-x2)%*%b
    Wdx <- NULL
    for (i in 1:length(b)){
      Wdx[i] <- (x1[i] - x2[i])*b[i] / A
    }
    return(Wdx)
  }
  
  Wdb.F  <- function(b1,b2,x) {
    # Yun weight function (C coefficients)
    # b1, b2 = coef
    # x      = weighted mean
    #
    A <- x%*%(b1-b2)
    Wdb <- NULL
    for (i in 1:length(x)){
      Wdb[i] <- (x[i]*(b1[i] - b2[i])) / A
    }
    return(Wdb)
  }
  
  dW.F <- function(b,x1,x2) {
    # derivative of Wdx  (K x K result) 
    # b      = coef
    # x1, x2 = weighted mean
    #
    dW <- array(NA,c(length(b), length(b)))
    A <- (x1-x2)%*%b   
    for (k in 1:length(b)) {
      for (l in 1:length(b)) {
        dW[k,l]  <-  as.numeric(k==l) * ((x1[k] - x2[k])/A) - 
          (b[k] * ( (x1[k] - x2[k])*(x1[l] - x2[l]) )/A^2 )
      }
    }
    return(dW)
  }     
  
  dwA.F <- function(b1,b2,x2) {
    # derivative of Wdb part A --  K x K  
    # b1, b2 = coef
    # x2     = weighted mean
    #
    dwA1 <- array(NA, c(length(b1), length(b2)))
    A <- x2%*%(b1-b2)
    for (k in 1:length(b1)){
      for (l in 1:length(b2)) {
        dwA1[k,l] <- as.numeric(k==l) * (x2[k]/A) - (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 
      } 
    }
    return(dwA1)
  }
  
  dwB.F <- function(b1,b2,x2) { 
    # derivative of Wdb part B
    # b1, b2 = coef
    # x2      = weighted mean
    #  
    dwB1 <- array(NA, c(length(b1), length(b2)))
    A <- x2%*%(b1-b2)
    for (k in 1:length(b1)){
      for (l in 1:length(b2)) {
        dwB1[k,l] <- (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 - as.numeric(k==l) * (x2[k]/A) 
      } 
    }
    return(dwB1)
  }
  
  wbA  <-  dwA.F(bA,bB,mB)
  wbB  <-  dwB.F(bA,bB,mB)
  #############
  # Yun weights
  #############
  Wdx  <-  Wdx.F(bA, mA, mB)
  Wdb  <-  Wdb.F(bA, bB, mB)
  ############################################
  # Convention: 1st moment A - 1st moment B
  # group A is treated as reference (standard)
  # group B is treated as comparison
  ############################################
  # decomp total: E & C
  #####################
  # E   <- svymean(~F(bA, xA), Asub)/mean(wA) -  
  #        svymean(~F(bA, xB), Bsub)/mean(wB)
  E   <- sum(F(bA, xA, wA))/sum(wA) -  
         sum(F(bA, xB, wB))/sum(wB)
  # Coefficients + Unexplained: 
  # C   <- svymean(~F(bA, xB), Bsub)/mean(wB) -  
  #        svymean(~F(bB, xB), Bsub)/mean(wB)
    C   <-  sum(F(bA, xB, wB))/sum(wB) -  
            sum(F(bB, xB, wB))/sum(wB)
  # get dEdb for variance estimator
  dWx  <- dW.F(bA, mA, mB)
  dEdb <- array(NA, c(length(bA), length(bA)))
  for (k in 1:length(bA)) {
    for (l in 1:length(bA)) {
      # 
#       df <- svymean(~as.vector(f(wA)*xA[,l]), Asub) / mean(wA) - 
#             svymean(~as.vector(f(wB)*xB[,l]), Bsub) / mean(wA)
       df <- mean(f(wA)*xA[,l])/mean(wA) - 
             mean(f(wB)*xB[,l])/mean(wB)
      dEdb[k,l] <- Wdx[k] * df + dWx[k,l] * E
    }
  }
  
  ##   Variances
  #    Composition (E)
  Var.E.k <- dEdb%*%varbA%*%t(dEdb)
  seWdx   <- sqrt(diag(Var.E.k))
  dCdbA   <- array(NA, c(length(bA), length(bA)))
  dCdbB   <- array(NA, c(length(bB), length(bB)))
  for (k in 1:length(bA)){
    for (l in 1:length(bB)) { 
 #    dCdbA[k,l] <- Wdb[k] * svymean(~as.vector(f(wB)*xB[,l]), Bsub)/mean(wB) + wbA[k,l]*C
 #    dCdbB[k,l] <- wbB[k,l] * C - Wdb[k]*svymean(~as.vector(f(wB)*xB[,l]),Bsub)/mean(wB)
      dCdbA[k,l]  <- Wdb[k] * mean(f(wB)*xB[,l])/mean(wB) + wbA[k,l]*C
      dCdbB[k,l]  <- wbB[k,l] * C - Wdb[k] *mean(f(wB)*xB[,l])/mean(wB)
    }
  } 
  
  ##  Variances
  #   Coefficients (C)
  Var.C.k <- dCdbA%*%varbA%*%t(dCdbA) + dCdbB%*%varbB%*%t(dCdbB)
  seWdb   <- sqrt(diag(Var.C.k))
  
  # get asymptotic variance of E_k and C_k
  # E_k
  dEdb.0 <- NULL
  for (k in 1:length(bA)){
#    dEdb.0[k] <- svymean(~as.vector(f(wA)*xA[,k]), Asub) / mean(wA) - 
#                svymean(~as.vector(f(wB)*xB[,k]), Bsub) / mean(wB)
   dEdb.0[k] <-  mean(f(wA)*xA[,k])/mean(wA) - 
                 mean(f(wB)*xB[,k])/mean(wB)
    
  }
  
  #C_k
  dCdb.0A <- NULL
  for (k in 1:length(bA)){
#      dCdb.0A[k] <- svymean(~as.vector(f(wB)*xB[,k]),Bsub) / mean(wB)
      dCdb.0A[k] <- mean(f(wB)*xB[,k])/mean(wB)
  }
  dCdb.0B <- NULL
  for (k in 1:length(bB)){
#      dCdb.0B[k] <- svymean(~as.vector(f(wB)*xB[,k]),Bsub) / mean(wB)
      dCdb.0B[k] <- mean(f(wB)*xB[,k])/mean(wB)
  }
  
  #         1 x k   k x k     k x 1
  var.E0 <- t(dEdb.0)%*%varbA%*%dEdb.0
  s.E0   <- sqrt(var.E0)
  
  var.C0 <- t(dCdb.0A)%*%varbA%*%dCdb.0A + t(dCdb.0B)%*%varbB%*%dCdb.0B 
  s.C0   <- sqrt(var.C0)
  
  R <- E + C
  v.names <- names(bA)
  if (is.null(v.names)) {
    v.names <- paste("Var", 1:length(bA), sep = "") # if names not available
  }
  # 
  
  ymA <- sum(F(bA, xA, wA))/sum(wA) * scale
  ymB <- sum(F(bB, xB, wB))/sum(wB) * scale
  
  if (printit==TRUE) {
    cat("","\n")
    cat("Linear Regression Decomposition","\n")
    cat("","\n")
    cat("Mean Group A           ", formatC(ymA, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Mean Group B           ", formatC(ymB, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("__________________________________________________________________", "\n")
    cat("Total Difference (E + C)", formatC((E+C)*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Raw Contrib..... (E,  C)", formatC(E*scale, flag=" ", dig=4, wid=10, format="f"), 
        formatC(C*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Std. Errors......(E,  C)", formatC(s.E0*scale, flag=" ", dig=4, wid=10, format="f"),
        formatC(s.C0*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Z-values.........(E,  C)", formatC(E/s.E0, flag=" ", dig=4, wid=10, format="f"),
        formatC(C/s.C0, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Pct Contrib......(E,  C)", formatC(100*E/(E+C), flag=" ", dig=4, wid=10, format="f"),
        formatC(100*C/(E+C), flag=" ", dig=4, wid=10, format="f"), "\n" )
    cat("","\n")
    cat("SUMMARY MEASURES","\n")
    cat("        Comparison Group", "\n")
    cat("        ----------------", "\n")
    cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
    cat("__________________________________________________________________", "\n")
    for (i in 1:length(mA)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(mA[i],flag=" ",dig=4,wid=10, format="f"),
          formatC(bA[i],wid=10,dig=4, format="f"),
          formatC(sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), 
          formatC(bA[i]/sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), "\n")
    }
    cat("","\n")
    cat("        Reference Group", "\n")
    cat("        ---------------", "\n")
    cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
    cat("__________________________________________________________________", "\n")
    for (i in 1:length(mB)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(mB[i],flag=" ",   dig=4,wid=10, format="f"),
          formatC(bB[i],wid=10,     dig=4, format="f"),
          formatC(sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), 
          formatC(bB[i]/sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), "\n")
    }
    cat("","\n")
    cat("DUE TO DIFFERENCE IN COMPOSITION (E)", "\n")
    cat("","\n")
    cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "        PCT",  "\n")
    cat("___________________________________________________________________", "\n")
    
    idx <- grepl("(Intercept)", v.names)
    
    for (i in 1:length(bA)){
      if (!idx[i]) {
        cat(formatC(v.names[i],format="s", wid=20),
            formatC(E*Wdx[i]*scale,       flag=" ", dig=4, wid=10, format="f"),
            formatC(seWdx[i]*scale,       flag=" ", dig=4, wid=10, format="f"),
            formatC(E*Wdx[i]/seWdx[i],    flag=" ", dig=4, wid=10, format="f"),
            formatC(100*(E*Wdx[i]/(E+C)), flag="+", dig=1, wid=10, format="f"), "\n")
      }
    }
    
    cat("","\n")
    cat("DUE TO DIFFERENCE IN COEFFICIENTS (C)", "\n")
    cat("","\n")
    cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "         PCT"  ,"\n")
    cat("___________________________________________________________________", "\n")
    for (i in 1:length(bA)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(C*Wdb[i]*scale,       flag=" ",  dig=4, wid=10, format="f"),
          formatC(seWdb[i]*scale,       flag=" ",  dig=4, wid=10, format="f"),
          formatC(C*Wdb[i]/seWdb[i],    flag=" ",  dig=4, wid=10, format="f"), 
          formatC(100*(C*Wdb[i]/(E+C)), flag="+",  dig=1, wid=10, format="f"),"\n")
      
    }
  } # end printit
  
  K        <- length(bA)
  # store these
  bE       <- matrix(E*Wdx,K,1) 
  varbE    <- matrix(Var.E.k,K,K)
  bC       <- matrix(C*Wdb,K,1)
  varbC    <- matrix(Var.C.k,K,K)
  E        <- E[[1]]
  varE     <- var.E0
  C        <- C[[1]]
  varC     <- var.C0
  rownames(bE) <- v.names
  rownames(bC) <- v.names
  colnames(varbE) <- rownames(varbE) <- v.names
  colnames(varbC) <- rownames(varbC) <- v.names
  
  # returns 
  return(list(vnames=v.names, E=E, varE=varE, 
              C=C, varC=varC,
              bE=bE, varbE=varbE, 
              bC=bC, varbC=varbC))
  
}  # end decomp_linear

#############

# END MAIN

#############