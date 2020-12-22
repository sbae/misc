

## Plot the object derived from tune.rsfrc
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

