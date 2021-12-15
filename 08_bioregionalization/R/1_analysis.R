# Paleobiogeographic analysis of PBDB data (last 10 Ma) 
# Adam T. Kocsis (Erlangen, 2020-06-17)
# CC-BY 4.0

library(icosa)
library(vegan)
library(igraph)
library(chronosphere)

# establish working directory
workdir <- file.path(Sys.getenv("Teaching"), "2020-06-17_biogeo")
setwd(workdir)

####################################################################################
# read in the data
# color scheme
load("data/allHex.RData")

# pbdb dataset prepared earlier
ceno6 <- readRDS("export/pbdb_species_49.rds")
# land polygons
land <- fetch("NaturalEarth")

source("R/methods/plots.R")

####################################################################################

# Geographic griding (icosa)

# create a grid
gr <- hexagrid(7, sp=TRUE)

	# visualize
	plot(gr)
	gridlabs(gr)

# locate the cells
ceno6$cells <- locate(gr, ceno6[, c("paleolng", "paleolat")])

# omit those entries where there is no species name or cells
ceno6use <- ceno6[!is.na(ceno6$cells) & !is.na(ceno6$trinomen),]

# contingency matrix
cont <- table(ceno6use$cells, ceno6use$trinomen)
samp <- cont[1:100,1:100]

# incidence
cont[cont>1] <- 1

# Method 1. Compositional similarity
distmat <- vegdist(cont, method="jaccard")

# clustering
cluster <- hclust(distmat, "ward.D")

# plot this
par(mfrow=c(2,1))
plot(cluster)
abline(h=1.2, col="red")

# cutting the dendrogram-> membership vector
mem <- cutree(cluster, h=1.2)

biogeoplot(mem)

# try setting it to 
# clustering method set to "average"
# cutting height set to 0.99


####################################################################################
dev.off()
# Network- approach
# same contingency matrix 

# transform to a bipartite network
bipartite <- graph_from_incidence_matrix(cont)

# # you can cluster this directly
# infoBi <- cluster_infomap(bipartite)
# 
# # membership 
# info <- membership(infoBi)
# 
# # localities 
# binfoLoc <- info[names(info)%in%rownames(gr@faces)]
# 
# # with network anaysis
# biogeoplot(binfoLoc)

# project the graph (to look at only localities)
graph <- bipartite_projection(bipartite, which="false")
plot(graph)

# you can cluster this directly
infoGraph <- cluster_infomap(graph)

# membership 
info <- membership(infoGraph)

# with network anaysis
biogeoplot(info)

# outliers?

#######################################################
# Simplified code development:
# https://github.com/adamkocsis/obigeo