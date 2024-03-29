library(readr)

## CKID specific: load & preprocess datasets
read_ckid <- function(cdbname, filename=NULL, ver = '/common', keep = NA, drop = NA) {
    keep <- tolower(keep)
    drop <- tolower(drop)

    cdbloc <- paste0('../data/codebook', ver, '/', cdbname, '.ndx')
    if (is.null(filename)) filename <- cdbname
    dataloc <- paste0('../data/data', ver,'/', filename, '.data')
    skip <- grep('----------', readLines(cdbloc))

    line_name <-
    read_table(cdbloc, col_names = FALSE, skip = skip) %>%
    select(X1, X2, X5)

    outfile <-
    as.data.frame(
        read_fwf(dataloc, 
            fwf_widths(c(diff(line_name$X1), line_name$X2[length(line_name$X2)]), tolower(line_name$X5)),
            guess_max=1e6, na = ".") 
        %>% select_if(function(x){!all(is.na(x))}))
	
    outfile <- na_if(outfile, '.')
	
    if (!is.na(keep)) {
 	   outfile <- outfile[, keep]
    } else if (!is.na(drop)) {
 	   outfile <- as.data.frame(outfile %>% select(-drop))
    } else {
  	  outfile <- outfile
    }
    return(outfile)
  }


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


## KIDMAC style conversion
ymd_to_decimal <- function(timestring){
    as.numeric(round((as.Date(timestring) - as.Date("1960-1-1"))/365.25+1960, 3))
}
decimal_to_ymd <- function(input){
    as.Date(round((input-1960)*365.25), origin='1960-1-1')
}
			      
## Alt way. Perhaps a bit more robust?
decimal_date <- function(timestring, format="%m/%d/%Y"){
    y <- (strptime(timestring, format=format)$year+1900)
    days_in_the_year <- strptime(paste0("12/31/",y), format=format)$yday+1
    yday <- (strptime(timestring, format=format)$yday/days_in_the_year)
    round(y + yday, 3)
}


## Plot the object derived from tune.rsfrc
library(akima)
plot.tune <- function(o, linear = TRUE) {
	x <- o$results[,1]
	y <- o$results[,2]
	z <- o$results[,3]
	so <- interp(x=x, y=y, z=z, linear = linear)
	idx <- which.min(z)
	x0 <- x[idx]
	y0 <- y[idx]
	filled.contour(x = so$x,
	y = so$y,
	z = so$z,
	xlim = range(so$x, finite = TRUE) + c(-2, 2),
	ylim = range(so$y, finite = TRUE) + c(-2, 2),
	color.palette =
	colorRampPalette(c("yellow", "red")),
	xlab = "nodesize",
	ylab = "mtry",
	main = "OOB error for nodesize and mtry",
	key.title = title(main = "OOB error", cex.main = 1),
	plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
	points(x,y,pch=16,cex=.25)})
}


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
			      
			      
## Estimate urine Albu-Protein ratio (Ref: Schneider AJKD 2021. https://pubmed.ncbi.nlm.nih.gov/32898620/)
calc_apr <- function(upcr, g, age){   ## upcr takes mg/g. CKiD default is g/g. May need to multiply by 1000.
    ifelse(g==TRUE,
        ## G
        ifelse(upcr>=1000,
            0.724,      # G, upcr>=1000
            1/(1 + 0.382*(upcr/1000)^-0.579)),      # G, upcr<1000
        
        ## NG
        ifelse(upcr>=1000,
            1/(1 + 0.642*0.906^((age-15)*(age<15))),      # NG, upcr>=1000
            ifelse(upcr>=100,
                1/(1 + 0.642*(upcr/1000)^(-0.720)*0.906^((age-15)*(age<15))),     ## NG, 100<=upcr<1000
                1/(1+ 3.369*0.906^((age-15)*(age<15)))         ## NG, upcr<100
            )
        )
    )
}			      
		      
