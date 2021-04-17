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
        select(X2, X5)
      outfile <-
        as.data.frame(read_fwf(dataloc, fwf_widths(line_name$X2, tolower(line_name$X5)),
                               na = ".") %>% select_if(function(x){!all(is.na(x))}))
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

preproc <- function(df, n.levels=10, na.values=c(-1,-9)) {
		summary(df)

		# Convert na.values to NA
		if(!is.null(na.values)){
			for(i in 1:length(na.values)) {
				df[df==na.values[i]] <- NA
			}
		}

		# define factors
		tgt <- sapply(df, function(x) length(unique(x))<n.levels)
		df[tgt] <- lapply(df[tgt], factor)
		print(names(df[tgt]))
		return(df)
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

