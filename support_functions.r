lincom <- function(model, varnum, multiplier=NULL, ci.level=0.95){
  # If multiplier is NULL, assume all are +1
  if(is.null(multiplier)) multiplier = rep(1, length(varnum))

  # get the lincom variance via a bilinear sum  
  out_var <- 0
  for(i in varnum){
    for(j in varnum){
      out_var <- out_var + vcov(model)[i,j]*multiplier[match(i, varnum)]*multiplier[match(j, varnum)]
    }
  }

  # calculate the lincom coef
  coef <- sum(coef(model)[varnum] * multiplier)

  # calculate Z
  z <- coef / sqrt(out_var)

  # Prepare output
  return(
    c(
      coef = coef, 
      se = sqrt(out_var),
      lci = coef + qnorm((1-ci.level)/2)*sqrt(out_var),
      uci = coef + qnorm((1+ci.level)/2)*sqrt(out_var),
      e_coef = exp( coef ), 
      e_lci = exp( coef + qnorm((1-ci.level)/2)*sqrt(out_var) ),
      e_uci = exp( coef + qnorm((1+ci.level)/2)*sqrt(out_var) ),
      z = z,
      p = pnorm(-abs(z))*2
    )
  )
}



nyu <- function(var){
  var[var=="N"] <- "No"
  var[var=="Y"] <- "Yes"
  var[var=="U"] <- NA
  var[var==""] <- NA
  return(factor(var))
}



lspline <- function( x, knots=NULL, marginal=FALSE, names=NULL ) {
  if(!is.null(names)) {
    .NotYetUsed("names")
  }
  n <- length(x)
  nvars <- length(knots) + 1
  # if( length(knots) > n/2)
  #   stop("too many knots")
  namex <- deparse(substitute(x))
  knots <- sort(knots)
  if(marginal) {
    rval <- cbind(
      x,
      sapply(knots, function(k) ifelse((x - k) > 0, x-k, 0) )
    )
  } else {
    rval <- matrix(0, nrow = n, ncol=nvars)
    rval[,1] <- pmin(x, knots[1])
    rval[,nvars] <- pmax(x, knots[length(knots)]) - knots[length(knots)]
    if(nvars > 2) {
      for(i in seq(2, nvars-1)) {
        rval[,i] <- pmax(  pmin(x, knots[i]), knots[i-1] ) - knots[i-1]
      }
    }
  }
  colnames(rval) <- seq(1, ncol(rval))
  structure(
    rval,
    knots = knots,
    marginal = marginal,
    class = c("lspline", "matrix")
  )
}


## Convert parameters from survreg/Weibull to flexsurvreg/Weibull or WeibullPH
survreg_to_flexsurvreg <- function(survreg_model, output="WeibullPH"){
    ## survreg AFT to flex AFT
    shape_flex_AFT <- 1/survreg_model$scale
    scale_flex_AFT <- exp(survreg_model$coefficients[1])
    coef_flex_AFT <- survreg_model$coefficients[-1]

    ## flex AFT to flex PH
    shape_flex_PH <- shape_flex_AFT
    scale_flex_PH <- scale_flex_AFT ^(-shape_flex_AFT)
    coef_flex_PH <- -coef_flex_AFT * shape_flex_AFT

    if(is.null(output) | output=="WeibullPH") list(shape=shape_flex_PH, scale=as.numeric(scale_flex_PH), coef=coef_flex_PH)
    else list(shape=shape_flex_AFT, scale=as.numeric(scale_flex_AFT), coef=coef_flex_AFT)
}

## Use the params from flexsurvreg/WeibullPH to calculate survival and quantiles
survivalPH <- function(time, scale, shape) exp(-scale * time^shape)
quantilePH <- function(survival, scale, shape) (log(survival)/(-1*scale))^(1/shape)

