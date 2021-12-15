#' Plotting biogeographic membership
#' Adam T. Kocsis (Erlangen, 2020-06-17)
#' CC-BY 4.0
#' @param mem Membership vector.
#' @param cols color vector.
#' @param bg Spatial object, background.
#' @param alpha alpha values of region colors
#' @param labels should the labels be plotted
#' @param gri icosa grid used for plotting
biogeoplot <- function(mem, cols=allHex, bg=land, alpha="99", labels=TRUE, gri=gr){
	# empty vector for the colors
	member <- rep(NA, nrow(gri@faces))

	# every entry corresponds to a face, ordered in the grid
	names(member) <- rownames(gri@faces)

	# color every entry
	reorder <- cols[mem] # implies: names(reorder) <- names(mem)

	# assign colors to appropriate face
	member[names(mem)] <- paste0(reorder, alpha)

	# plot empty background
	plot(bg, col="gray70")

	# plot colors
	plot(gri, col=member, add=TRUE)

	if(labels){
		# centroids reordered to match the memberhsip vector
		cent <- centers(gri)[names(mem),]

		# plot the membership
		text(x=cent[,1], y=cent[,2], label=mem, cex=0.6)
	}
}
