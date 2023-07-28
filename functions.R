
## Automatically convert numeric variables to factors based on the number of unique values
preproc <- function(df, n.levels=10, na.values=c(-1,-9), no_factor=NULL) {
		summary(df)

		# Convert na.values to NA
		if(!is.null(na.values)){
			for(i in 1:length(na.values)) {
				df[df==na.values[i]] <- NA
			}
		}

		# define factors
		tgt <- sapply(df, function(x) length(unique(x))<n.levels)
		tgt[names(df) %in% no_factor] <- FALSE
		df[tgt] <- lapply(df[tgt], factor)
		print(names(df[tgt]))
		return(df)
}

		      

## Plot the object derived from tune.rsfrc
# library(akima)
# plot.tune <- function(o, linear = TRUE) {
# 	x <- o$results[,1]
# 	y <- o$results[,2]
# 	z <- o$results[,3]
# 	so <- interp(x=x, y=y, z=z, linear = linear)
# 	idx <- which.min(z)
# 	x0 <- x[idx]
# 	y0 <- y[idx]
# 	filled.contour(x = so$x,
# 	y = so$y,
# 	z = so$z,
# 	xlim = range(so$x, finite = TRUE) + c(-2, 2),
# 	ylim = range(so$y, finite = TRUE) + c(-2, 2),
# 	color.palette =
# 	colorRampPalette(c("yellow", "red")),
# 	xlab = "nodesize",
# 	ylab = "mtry",
# 	main = "OOB error for nodesize and mtry",
# 	key.title = title(main = "OOB error", cex.main = 1),
# 	plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
# 	points(x,y,pch=16,cex=.25)})
# }


## Wrapper for xgb.cv
xgb_cv_wrap <- function(x, data, nfold=5) {
  cv <- xgb.cv(
          params = x,
          data = data,
          nrounds = 2000,
          early_stopping_rounds = 50,
          print_every_n = 50, nfold = nfold)
  out <- data.frame(x, eval=as.numeric(cv$evaluation_log[cv$best_iteration,4]), best_iter=cv$best_iteration)
  return(out)
}

			      

### Just like stata's lincom
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
  z <- coef/sqrt(out_var)

  # Prepare output
  return(
    c(
      coef = coef, 
      se = sqrt(out_var),
      lci = coef + qnorm((1-ci.level)/2)*sqrt(out_var),
      uci = coef + qnorm((1+ci.level)/2)*sqrt(out_var),
      z = z,
      p = pnorm(-abs(z))*2
    )
  )
}
			      

## Van Calster calibration curves -- my implementation
calicurve <- function(pred, obs, graph = TRUE) {
    m0 <- glm(obs ~ 0 + offset(qlogis(pred)), family = binomial)
    m1 <- glm(obs ~ offset(qlogis(pred)), family = binomial)
    m2 <- glm(obs ~ qlogis(pred), family = binomial)

    p <- c(
        1 - pchisq(m0$deviance - m2$deviance, 2),
        1 - pchisq(m0$deviance - m1$deviance, 1),
        1 - pchisq(m1$deviance - m2$deviance, 1)
    )
    names(p) <- c("H0_1", "H0_2", "H0_3")
    print(p)

    if (graph != FALSE) {
        m3 <- predict(loess(obs ~ pred), se = T)
        plot(NULL, xlim = 0:1, ylim = 0:1, asp = 1, xlab = "Predicted Probability", ylab = "Observed Proportion")
        polygon(c(sort(pred), rev(sort(pred))),
            c(
                (m3$fit - qt(0.975, m3$df) * m3$se)[order(eval(pred))],
                (m3$fit + qt(0.975, m3$df) * m3$se)[rev(order(eval(pred)))]
            ),
            border = NA, col = "grey90"
        )
        lines(sort(pred), m3$fit[order(pred)], lwd = 3, xlim = 0:1)
        abline(0, 1, xlim = 0:1, col = "grey50")
    }

    return(list(models = list(model0 = m0, model1 = m1, model2 = m2), p = p))
}

