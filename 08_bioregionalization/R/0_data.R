# Preparing species-level PBDB data (last 10 Ma) for paleobiogeographic analysis
# Adam T. Kocsis (Erlangen, 2020-06-17)
# CC-BY 4.0

library(chronosphere)
library(divDyn)

# establish working directory
workdir <- file.path(Sys.getenv("Teaching"), "2020-06-17_biogeo")
setwd(workdir)


#########################################################################
# 1. get the data
dat <- fetch("pbdb")

# resolve occurrences (as in divDyn)
# omit non-marine

#########################################################################
# 2. taxonomic filtering
# filter records not accepted at least to genus
dat <-dat[dat$accepted_rank %in% c("genus", "species"),]

# omit non-informative genus entries
dat <- dat[dat$genus!="", ]

# taxonomic subset (modern benthic)
dat <- dat[c(
	which(dat$class=="Gastropoda" ),
	which(dat$class=="Bivalvia" ),
	which(dat$phylum=="Bryozoa" ),
	which(dat$phylum=="Echinodermata" ),
	which(dat$phylum=="Brachiopoda" ),
	which(dat$order=="Scleractinia" ),
	which(dat$order=="Decapoda"))
	, ]

# resolve the potential homonymy problem
dat$clgen <- paste(dat$class, dat$genus)


################################################################################
# 2. filter by environment
levels(factor((dat$environment)))

omitEnv <- c(
	"\"floodplain\"",
	"alluvial fan",
	"cave",
	"\"channel\"",
	"channel lag" ,
	"coarse channel fill",
	"crater lake",
	"crevasse splay",
	"dry floodplain",
	"delta plain",
	"dune",
	"eolian indet.",
	"fine channel fill",
	"fissure fill",
	"fluvial indet.",
	"fluvial-lacustrine indet.",
	"fluvial-deltaic indet.",
	"glacial",
	"interdune",
	"karst indet.",
	"lacustrine - large",
	"lacustrine - small",
	"lacustrine delta front",
	"lacustrine delta plain",
	"lacustrine deltaic indet.",
	"lacustrine indet.",
	"lacustrine interdistributary bay",
	"lacustrine prodelta",
	"levee",
	"loess",
	"mire/swamp",
	"pond",
	"sinkhole",
	"spring",
	"tar",
	"terrestrial indet.",
	"wet floodplain"
)

dat<-dat[!dat$environment%in%omitEnv, ]


################################################################################
# 3. stratigraphic resolution
	data(keys)
	# time scales
	data(stages)
	data(tens)
	
	# rename the entries
	colnames(tens)[colnames(tens)=="X10"]<-"name"
	colnames(stages)[colnames(stages)=="stage"]<-"name"

	
#-------------------------------------------------------------------------------
# A.  the bin entries
		
	# a. categorize interval names to bin numbers
		# categorize is the new function of the package
		tenMin<-categorize(dat[,"early_interval"],keys$tenInt)
		tenMax<-categorize(dat[,"late_interval"],keys$tenInt)

		# convert to numeric
		tenMin<-as.numeric(tenMin)
		tenMax<-as.numeric(tenMax)

	# b. resolve max-min interval uncertainty
	# final variable (empty)
		dat$ten <- rep(NA, nrow(dat))

	# use entries, where
		tenCondition <- c(
		# the early and late interval fields indicate the same bin
			which(tenMax==tenMin),
		# or the late_interval field is empty
			which(tenMax==-1))

	# in these entries, use the bin indicated by the early_interval
		dat$ten[tenCondition] <- tenMin[tenCondition]

#####################################-------------------------------------------
# do the same for the stages
# B. the stg entries (lookup)
		stgMin<-categorize(dat[,"early_interval"],keys$stgInt)
		stgMax<-categorize(dat[,"late_interval"],keys$stgInt)

		# convert to numeric
		stgMin<-as.numeric(stgMin)
		stgMax<-as.numeric(stgMax)

	# empty container
		dat$stg <- rep(NA, nrow(dat))

	# select entries, where
		stgCondition <- c(
		# the early and late interval fields indicate the same stg
			which(stgMax==stgMin),
		# or the late_intervarl field is empty
			which(stgMax==-1))

	# in these entries, use the stg indicated by the early_interval
		dat$stg[stgCondition] <- stgMin[stgCondition]

################################################################################
# We only need 10 my bin (ten) No. 49
	datRec <- dat[which(dat$ten==49), ]


################################################################################
# C. Species-level taxonomy
	
	# genus level
	datRec<-datRec[!is.na(datRec$genus),]
	
	# wrong entries replaced by NA
	datRec$clgen[which(datRec$genus=="")] <- NA
	
	# the accepted name should be taken as default - if there is any
	spVec <- datRec$accepted_name
	
	# the identified names vector
	# which occs are identified only at the genus level?
	acceptGenus <- which(datRec$accepted_rank=="genus")
	
	# identified names should be copied over - or only the genus names will be there
	spVec[acceptGenus]<-datRec$identified_name[acceptGenus]

	# placeholder is to be put instead of those entries not identified to a species
	spVec[datRec$identified_rank=="genus"] <- "Dummy sp."
	
	# for sp level
	vec<-levels(factor(spVec))
	#vec<-as.character(vec[1:150])
	
# the new names table
	spNames<-cleansp(vec, debug=TRUE)
	
# render new names to occurrences
	newNames<-spNames$binomen
	names(newNames)<- vec
	
	# reconstruct original order
	pureSP <- newNames[spVec]

# add them to the table
	datRec$trinomen<-paste(datRec$class, pureSP, sep="_")
	datRec$trinomen[is.na(pureSP)]<- NA 

	# number of species level entries
	sum(!is.na(datRec$trinomen))
	
	# save data
	saveRDS(datRec, file="export/pbdb_species_49.rds")

