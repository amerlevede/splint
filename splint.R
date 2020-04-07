#!/usr/bin/Rscript
startTime = Sys.time()

superspline <- function(x, y, M, lambda, island=NULL, newx=x, newisland=island, discont=NULL, na.rm=TRUE, weights=NULL) {
	stopifnot(length(x)>0)
	stopifnot(length(x)==length(y))
	stopifnot(M %in% c(1, 2, 3))
	stopifnot(lambda>0)
	stopifnot(is.null(island) || length(island)==length(x))
	stopifnot(is.null(newisland) || length(newisland)==length(newx))
	stopifnot(all(is.finite(x)))
	stopifnot(is.null(weights) || length(x) == length(weights))
	stopifnot(all(is.finite(weights)))
	stopifnot(all(weights>0))

	# If weights are not given, set to 1
	if (is.null(weights)) {
		weights = rep(1,length(x))
	} else {
		weights = weights / sum(weights) * length(x)
	}

	# If islands are not given, calculate from given discontinuities or fit spline without discontinuities
	if (is.null(island)) {
		if (is.null(discont)) {
			island = rep("All",length(x))
			newisland = rep("All",length(newx))
		} else {
			island = colSums(dist(discont,x,function(x,y){x<y}))
			newisland = colSums(dist(discont,newx,function(x,y){x<y}))
		}
		return (superspline(x=x, y=y, island=island, M=M, lambda=lambda, newx=newx, newisland=newisland, na.rm=na.rm, weights=weights))
	}
	stopifnot(length(island)==length(x))

	# Ignore non-finite predictions
	if (na.rm) {
		return (superspline(x=x[is.finite(y)], y=y[is.finite(y)], island=island[is.finite(y)], M=M, lambda=lambda, newx=newx, newisland=newisland, na.rm=FALSE, weights=weights))
	} else {
		stopifnot(all(is.finite(y)))
	}

	# Normalize x
	n <- length(x)
	minx <- min(x)
	newx = (newx-min(x))/max(x)
	x = (x-min(x))/max(x)
	# Define basis functions
	phi = function(t,j) {t^j/factorial(j)}
	islands = unique(island)
	R = if (M==1) {
		pmin
	} else if (M==2) {
		function(s,t){pmax(s,t)*pmin(s,t)^2/2 - pmin(s,t)^3/6}
	} else if (M==3) {
		function(s,t){pmin(s,t)^5/30 - pmax(s,t)*pmin(s,t)^4/6 + pmax(s,t)^2*pmin(s,t)^3/3}
	}
	# Calculate values
	Sigma = outer(x, x, R)
	basis = function(t,i){cbind(if (M>1) outer(t,1:(M-1),phi) else NULL, outer(i,islands,"=="))}
	T = basis(x,island)
	# Solve solution
	Minv = solve(Sigma + n*lambda*diag(1/weights))
	beta = solve(t(T)%*%Minv%*%T)%*%t(T)%*%Minv%*%y # Coefficients of unpenalized functions
	coef = Minv%*%(y - T%*%beta) # Coefficients of evaluation representations
	# Calculate solution
	result = c(basis(newx,newisland)%*%c(beta) + outer(newx,x,R)%*%c(coef))
	attr(result, "islands") = c(beta)[M:length(beta)]
	names(attr(result,"islands")) = islands
	return(result)
}

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

quotient <- function(x,y) {pmax(x/y,y/x)}

harmonic.mean <- function(x) {1/mean(1/x)}


## Load packages
################################################################################################################################################################
################################################################################################################################################################
packlib = "~/Software/Rpackages/"
R_LIBS_USER = packlib
.libPaths(c(packlib, .libPaths()))
maskfun = function(x){suppressMessages(suppressWarnings(x))}
maskfun(library("ggplot2"))
maskfun(library("memoise"))
maskfun(library("reshape"))
maskfun(library("zoo"))
maskfun(library("grid"))
maskfun(library("splines"))
maskfun(library("quantmod"))
maskfun(library("HMM"))
maskfun(library("Hmisc"))
maskfun(library("proxy"))

## Input variables
################################################################################################################################################################
################################################################################################################################################################
stopifnot(length(commandArgs(T))>=1)
# input : input .tsv file
stopifnot(length(grep('input=(.*)', commandArgs(T)))==1)
input = strsplit(commandArgs(T)[grep('input=(.*)', commandArgs(T))==1],"=")[[1]][2]
# plot : output .pdf file to put cnv plot
stopifnot(length(grep('plot=(.*)', commandArgs(T)))<=1)
plot=F
if (length(grep('plot=(.*)', commandArgs(T)))==1) {
	plot = T
	plotfile = strsplit(commandArgs(T)[grep('plot=(.*)', commandArgs(T))],"=")[[1]][2]
}
# byplot : output .pdf file to put diagnostic plots plot
stopifnot(length(grep('(by)?plot=(.*)', commandArgs(T)))<=1)
byplot=F
if (length(grep('(by)?plot=(.*)', commandArgs(T)))==1) {
	byplot = T
	byplotfile = strsplit(commandArgs(T)[grep('(by)?plot=(.*)', commandArgs(T))],"=")[[1]][2]
}
# forceLine : if True, smiley is not fitted
forceLine = F
if (length(grep('^line$', commandArgs(T)))==1) {
	message("Forcing line ...")
	forceLine = T
}
# tab : output table file of calculated frame properties
stopifnot(length(grep('tab=(.*)', commandArgs(T)))<=1)
writeTab=F
if (length(grep('tab=(.*)', commandArgs(T)))==1) {
	writeTab=T
	tabFile = strsplit(commandArgs(T)[grep('tab=(.*)', commandArgs(T))],"=")[[1]][2]
}
# islands : output table file of calculated island
stopifnot(length(grep('islands=(.*)', commandArgs(T)))<=1)
writeIslands=F
if (length(grep('islands=(.*)', commandArgs(T)))==1) {
	writeIslands=T
	islandsFile = strsplit(commandArgs(T)[grep('islands=(.*)', commandArgs(T))],"=")[[1]][2]
}
#chromosomes : chromosomes to keep or remove from analysis
stopifnot(length(grep('chr=(.*)', commandArgs(T)))<=1)
keepChr=NULL
rmChr="Mito"
if (length(grep('chr=(.*)', commandArgs(T)))==1) {
	chrarg=strsplit(commandArgs(T)[grep('chr=(.*)', commandArgs(T))],"=")[[1]][2]
	chrarg=strsplit(chrarg,",")[[1]]
	keepChr=chrarg[!grepl("^-",chrarg)]
	rmChr=substring(chrarg[grepl("^-",chrarg)],2)
}

## Settings
################################################################################################################################################################
################################################################################################################################################################
maskfun = function(x){suppressMessages(suppressWarnings(x))}
# maskfun = identity

# options(warn=1)

telomereLength 		= 50000		# Frames are telomereLength from the sides of a scaffold will be considered telomeric
								# These regions will not be used for several analysis steps because they are too chaotic
breakSize			= 50000		# k parameter for robust derivative calculation
# breakTreshold		= 0.12		# Peaks in the robust derivative will be identified as breakpoint when the abs exceeds this value
								# Standard deviation of high-coverage sample of reference genome (i.e. no CNVs) without chromosome XII or Mito is 0.01166956
mindist 			= breakSize # Minimum distance allowed between large breakpoints
maxlouche 			= 0.85		# frames with more than maxlouche*100% "louche" mapped reads will be marked as outliers
								# (This calculation is no longer used)
ratioMergeTreshold 	= 0.15 		# adjacent islands with less than ratioMergeTreshold difference in depth ratio (ratio of depth in that region compared to the genomic average depth) will be merged
outlierWeight		= 0.1		# this weight will be given to outliers in the baseline spline fitting procedure


## Import data
################################################################################################################################################################
################################################################################################################################################################
message("Importing data...")

# read input table
# contains reading frames with headers (at least): scaffold, frameStart, frameEnd, depth, louchedepth
fdata 				= read.table(input, header=T)
stopifnot(all(c("scaffold","frameStart","frameEnd","depth")%in%colnames(fdata))) ## use of "louchedepth" depricated
fdata[,"scaffold"]  = factor(fdata[,"scaffold"], levels=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI","Mito"))
fdata 				= fdata[order(as.numeric(fdata[,"scaffold"]),fdata[,"frameStart"]),]
fdata[,"scaffold"] 	= as.character(fdata[,"scaffold"])

# remove specific scaffolds
if (length(keepChr)>=1) fdata = fdata[fdata[,"scaffold"]%in%keepChr,]
fdata 				= fdata[!(fdata[,"scaffold"]%in%rmChr),]
fdata[,"scaffold"] 	= factor(fdata[,"scaffold"], levels=unique(fdata[,"scaffold"]))

# calculate basic information from input data
frameStep 			= min(diff(sort(fdata[fdata[,"scaffold"]==fdata[1,"scaffold"],"frameStart"])))
frameSize 			= max(fdata[,"frameEnd"]-fdata[,"frameStart"])
fdata[,"louchefrac"]= fdata[,"louchedepth"]/fdata[,"depth"]
fdata[!is.finite(fdata[,"louchefrac"]),"louchefrac"] = 1
fdata[,"depth"]		= fdata[,"depth"]/frameSize

basedepth 			= median(fdata[,"depth"])
depthmin 			= quantile(fdata[,"depth"], .1)
depthmax 			= quantile(fdata[,"depth"], .9)

# calculate basic properties of scaffolds
annotateScaffolds <- function(scaffolds=fdata[,"scaffold"], starts=fdata[,"frameStart"], ends=fdata[,"frameEnd"]) {
	n=length(scaffolds)
	stopifnot(n>0)
	stopifnot(n==length(starts))
	stopifnot(n==length(ends))
	stopifnot(!is.unsorted(as.numeric(factor(scaffolds))))

	# calculate properties
	result = data.frame(name=unique(scaffolds))
	result[,"name"]  = as.character(result[,"name"])
	rownames(result) = result[,"name"]
	result[,"start"] = sapply(result[,"name"], function(x){which.max(scaffolds==x)})
	result[,"end"]   = sapply(result[,"name"], function(x){n-which.max(rev(scaffolds)==x)+1})
	result[,"length"]= ends[result[,"end"]]

	result[,"centromere"] = sapply(result[,"name"], function(x){
		switch(x,
			"I"		= (151466+151582)/2/1000000,
			"II"	= (238208+238323)/2/1000000,
			"III"	= (114386+114501)/2/1000000,
			"IV"	= (449712+449821)/2/1000000,
			"V"		= (151988+152104)/2/1000000,
			"VI"	= (148511+148627)/2/1000000,
			"VII"	= (496921+497038)/2/1000000,
			"VIII"	= (105587+105703)/2/1000000,
			"IX"	= (355630+355745)/2/1000000,
			"X" 	= (436308+436425)/2/1000000,
			"XI"	= (440130+440246)/2/1000000,
			"XII"	= (150829+150947)/2/1000000,
			"XIII"	= (268032+268149)/2/1000000,
			"XIV"	= (628759+628875)/2/1000000,
			"XV"	= (326585+326702)/2/1000000,
			"XVI"	= (555958+556073)/2/1000000,
			0)
		})

	return (result)
}
scaffolds = annotateScaffolds()

# add new calculated columns to data
# position is average position of frame
# telomere is boolean and true iff the frame is in a telomeric region, as specified by the "setting" telomereLength
#	telomeres are often left out because they are highly variable
# outlier is a boolean that specifies which frames should be used in some calculations
#	outliers are one or more of the following types:
#		* frames located in highly variable regions
# 		* frames with especially low or high read depth
# (this calculation was removed)		* frames where many of the fraction of "louche" mapped reads exceeds the set maximum
fdata[,"position"]	= ceiling((fdata[,"frameEnd"] + fdata[,"frameStart"])/2)/10^6
fdata[,"telomere"] 	= fdata[,"frameEnd"] < telomereLength | fdata[,"frameStart"] > scaffolds[fdata[,"scaffold"],"length"] - telomereLength
fdata[,"deviation"]	= rollapply(fdata[,"depth"],floor(12500/frameSize),sd,partial=T)
medianLocalSd 		= median(fdata[,"deviation"])
fdata[,"outlier"] 	= fdata[,"deviation"]>2*medianLocalSd |
						fdata[,"depth"] < 1/4*basedepth |
						fdata[,"depth"] > max(c(2*basedepth,1.5*depthmax)) |
						fdata[,"telomere"]
						# fdata[,"louchefrac"] < maxlouche

# clean chromosome IX
# fdata[fdata[,"scaffold"]=="IX" & fdata[,"position"]<0.05, "depth"] = 0

# check data integrity
stopifnot(all(fdata[,"depth"] >= 0))
stopifnot(all(diff(fdata[,"frameStart"])==frameStep | diff(fdata[,"frameStart"]) < 0))
stopifnot(all(fdata[,"frameEnd"] - fdata[,"frameStart"] <= frameSize))
stopifnot(sum(fdata[,"outlier"]) > 0.75 * ncol(fdata))
for (scaf in scaffolds[,"name"]) stopifnot(fdata[which.max(fdata[,"scaffold"]==scaf),"frameStart"] == 0)

# output some stats
message(paste("   coverage depth:",paste(depthmin,basedepth,depthmax,sep=", ")))
message(paste("   frameSize:",frameSize))
message(paste("   frameStep:",frameStep))
message(paste("   median local sd:",median(fdata[!fdata[,"outlier"],"deviation"])))

# Create plot file
if (byplot) {
	graphics.off()
	pdf(byplotfile, width=5*length(scaffolds[,"name"]))
}

# Plot initial read depth landscape
if (byplot) maskfun(print(maskfun(ggplot(fdata, aes(x=position, y=depth, color=outlier))
			+ geom_point()
			+ facet_grid(.~scaffold, scales="free_x", space="free_x")
			+ theme(panel.margin=unit(0,"lines"))
			+ theme_bw()
			+ ggtitle("Initial read depth landscape and outliers"))))

# Plot outlier detection
if (byplot) {
	outliersdtreshold = median(rollapply(fdata[,"depth"],floor(25*500/frameSize),sd,partial=T))
	maskfun(print(maskfun(ggplot(fdata, aes(x=position, y=rollapply(depth,floor(25*500/frameSize),sd,partial=T)))
		+ geom_point()
		+ geom_hline(aes(yintercept=2*outliersdtreshold), color="red")
		+ facet_grid(.~scaffold, scales="free_x", space="free_x")
		+ theme_bw()
		+ ggtitle("Local deviation and outlier treshold"))))
}

## Locate large breakpoints
################################################################################################################################################################
################################################################################################################################################################
# Calculate breakpoints based on robust log derivative of the data
message("Finding breakpoints...")
# minlogdepth = log2(basedepth/10)
minlogdepth = 0
# safely calculate log2 of depth
#	very low and undefined frames are given a lower bound value
log2depth <- function(depth) {
	result = log2(depth)
	result[!is.finite(result)] = minlogdepth
	result = pmax(minlogdepth,result)

	return (result)
}
flatBias = 5
depthquotient <- function(depth1,depth2) {
	high = pmax(depth1/2+depth2/2,pmax(depth1,depth2)-flatBias/2)
	low = pmin(depth1/2+depth2/2,pmin(depth1,depth2)+flatBias/2)
	return (2^(log2depth(high) - log2depth(low)))
}

# Calculate robust derivative
calculateRobdev <- function(depth, position, telomeres, outliers=rep(FALSE,length(depth)), forceLine.=forceLine) {
	stopifnot(length(depth)>0)
	stopifnot(all(depth>=0))
	stopifnot(all(position>0))
	stopifnot(length(depth)==length(position))
	stopifnot(length(depth)==length(telomeres))

	# Calculate robust log average of depth
	logdepth = log2depth(depth)
	logdepth[outliers] = NA
	robdev = rollapply(logdepth,list(0:(breakSize/frameStep)),median,partial=T,na.rm=T) - rollapply(logdepth,list((-breakSize/frameStep):0),median,partial=T,na.rm=T)

	robdev[telomeres] = 0
	robdev[!is.finite(robdev)] = 0

	# correct for linear bias in derivative-like measure
	if (!forceLine.) {
		slope = (lm(robdev~position,weights=1-telomeres)$coefficients)["position"]
		robdev = robdev - position*slope
	}
	robdev = robdev - median(robdev)
	robdev[telomeres] = 0
	return (robdev)
}

for (scaf in scaffolds[,"name"]) {
	fdata[fdata[,"scaffold"]==scaf,"robdev"] = calculateRobdev(
		fdata[fdata[,"scaffold"]==scaf,"depth"],
		fdata[fdata[,"scaffold"]==scaf,"position"],
		fdata[fdata[,"scaffold"]==scaf,"telomere"],
		fdata[fdata[,"scaffold"]==scaf,"outlier"])
}

# Smoothe out wrinkles and small deviations in robdev
# 	this allows to differentiate between single frames that exceed the treshold, and actual sizeable peaks
smootheWrinkles <- function(robdev) {
	robdevav = rollmean(robdev, breakSize/frameStep, fill=0)
	robdevav = robdevav - sign(robdevav)*mapply(min, abs(robdevav), breakTreshold)
	return (robdevav)
}

# Find breakpoints
findBreakpoints <- function(robdev, scaffold, telomeres) {
	stopifnot(length(robdev)==length(scaffold))
	stopifnot(all(is.finite(robdev)))

	# Find treshold value based on median of robust log-derivative
	breakTreshold <<- 2.5*median(abs(robdev))
	# breakTreshold <<- 0.15
	message(paste("   Median abs robust log-derivative:",breakTreshold/2.5))

	robdevav = smootheWrinkles(robdev)

	# Find valleys and peaks outside the treshold value
	breakpoint = rep(0,length(robdevav))
	breakpoint[findPeaks(robdevav)] = 1
	breakpoint[findValleys(robdevav)] = 1
	breakpoint[telomeres] = 0
	breakpoint 	= breakpoint & abs(robdevav) > breakTreshold

	# Remove breakpoints that are too close together
	indices = which(breakpoint)
	tresholdDifference = floor(mindist/frameStep)
	for (x in indices) {
		for (y in indices) {
			if (x==y) next
			if (scaffold[x] != scaffold[y]) next
			if (!breakpoint[x] || !breakpoint[y]) next
			if (abs(x-y) < tresholdDifference) {
				lowestIndex = c(x,y)[which.min(c(robdev[x], robdev[y]))]
				breakpoint[lowestIndex] = FALSE
			}
		}
	}

	return (breakpoint)
}
fdata[,"breakpoint"] = findBreakpoints(fdata[,"robdev"], fdata[,"scaffold"], fdata[,"telomere"])
message(paste("   Breakpoints:",paste(fdata[fdata[,"breakpoint"],"scaffold"],fdata[fdata[,"breakpoint"],"position"],sep=":",collapse=", ")))

# Define CN regions (islands)
for (scaf in scaffolds[,"name"]) fdata[fdata[,"scaffold"]==scaf,"island"] = paste(scaf,1+cumsum(fdata[fdata[,"scaffold"]==scaf,"breakpoint"]),sep=".")


# Output some stats
message(paste("   Islands found:", paste(unique(fdata[,"island"]), collapse=", ")))

# plot robust derivative
if (byplot) maskfun(print(maskfun(ggplot(fdata,aes(x=position,y=robdev, color=island))
		+ geom_line()
		+ geom_line(aes(y=smootheWrinkles(robdev)), linetype=2)
		+ geom_hline(yintercept=-breakTreshold, color="red")
		+ geom_hline(yintercept=breakTreshold, color="red")
		+ geom_vline(aes(xintercept=position*breakpoint),alpha=0.5,color="black")
		+ facet_grid(.~scaffold, scales="free_x", space="free_x")
		+ theme(panel.margin=unit(0,"lines"))
		+ theme_bw()
		+ ggtitle("robust derivative (solid), smoothed robust derivative (dashed) and breakpoint treshold"))))


## Grow small copy number regions
################################################################################################################################################################
################################################################################################################################################################
# Calculate Markov states representing copy number in a region
#	States: * CNV-: represents a lower-than-expected read depth according to given predictions
#			* CNV1: represents the expected read depth
#			* CNV+: represents a higher-than-expected read depth
#			* CNV0: represents a total deletion
#	Symbols:* obs>pred: the observed read depth is higher than the predicted value according to the baseline prediction
#			* obs=pred: the observed read depth is approximately equal to the predicted value (for some allowed ratio delta)
#			* obs<pred: the observed read depth is lower than the predicted value
CNVHMM <- function(observed, predicted, delta=1.2, p=0.05*frameSize/1000) {
	stopifnot(length(observed) == length(predicted))

	# define HMM properties
	#	p represents the `probability' (i.e. in this case penalty) of transitioning between different copy number states
	states 		= c("CNV-","CNV1","CNV+","CNV0")
	symbols 	= c("obs>pred","obs=pred","obs<pred")
	startProbs 	= rep(1/length(states),length(states))
	transProbs 	= diag(rep(1-length(states)*p,length(states)))+p
	emissionProbs = t(matrix(c(
						#+ 		#0 		#-
					 	0.05,	0.45,	0.5,	#CNV-
						0.15,	0.7,	0.15,	#CNV1
						0.5,	0.45,	0.05,	#CNV+
						0.05,	0.05,	0.9),	#CNV0
					ncol=length(states)	))
	HMM 		= initHMM(States=states, Symbols=symbols, startProbs=startProbs, transProbs=transProbs, emissionProbs=emissionProbs)


	# calculate observed symbols from read depths
	#	delta is the maximal quotient between the observed and expected values for them to be seen as approximately equal
	delta <<- delta
	obs 		= gsub("FALSE","obs<pred",gsub("TRUE","obs>pred",as.character(observed>predicted)))
	obs[depthquotient(observed,predicted)<delta] = "obs=pred"
	# obs[observed/predicted > 0 & (observed/predicted < delta & observed/predicted > 1/delta)] = "obs=pred"

	# calculate the most likely Markov path, i.e. copy number sequence
	vit 		= viterbi(HMM, obs)

	return (vit)
}

# find small copy number regions given a read depth sequence
# 	inference is based on HMM state compared to a spline prediction
# 	after finding smaller regions, the process is repeated 
findSmallIslands <- function(depth, position, outlier, islands=rep("0", length(depth)), forceLine.=forceLine) {
	stopifnot(length(depth) == length(position))
	stopifnot(length(depth) == length(outlier))
	stopifnot(all(is.finite(outlier)))
	stopifnot(all(is.finite(1 - (1-outlierWeight)*outlier)))

	# repeat until convergence
	splineRepetitions = 0
	repeat {
		splineRepetitions = splineRepetitions + 1
		stopifnot(splineRepetitions <= 50)

		# calculate predicted depth according to spline
		if (forceLine.) {
			depths = sapply(unique(islands), function(x){median(depth[islands==x])})
			names(depths) = as.character(unique(islands))
			pdepth = depths[as.character(islands)]
		} else {
			pdepth = pmax(1, superspline(
				x = position,
				y = depth,
				island = islands,
				# weights = (1 - (1-outlierWeight)*outlier),
				M=1, lambda=10^-2.5
				))
		}

		# find abnormal copy number regions
		viterbiState = CNVHMM(depth, pdepth)

		# define new subislands
		# newislands = paste(islands, cumsum(c(1,diff(as.numeric(factor(viterbiState)))!=0))*(viterbiState!="CNV1"), sep=".") # Regular regions get kept together
		newislands = paste(islands, cumsum(c(1,diff(as.numeric(factor(viterbiState)))!=0)), sep=".") # All regions are monolythic
		newislands = as.numeric(factor(newislands, levels=unique(newislands)))-1

		# # plot result
		# if (byplot) print(ggplot(data.frame(position=position, depth=depth, pdepth=pdepth, outlier=outlier, viterbiState=viterbiState), aes(x=position))
		# 	+ geom_point(aes(y=depth, color=viterbiState, shape=outlier))
		# 	+ geom_line(aes(y=pdepth), color="red")
		# 	+ geom_line(aes(y=pdepth*delta), color="red", linetype=2, alpha=0.5)
		# 	+ geom_line(aes(y=pdepth/delta), color="red", linetype=2, alpha=0.5)
		# 	# + geom_vline(xintercept=breaks, linetype=2, alpha=0.5)
		# 	)

		# continue if new subislands were found
		if (all(newislands == islands)) break
		islands = newislands
	}

	return (islands)
}


# Define subislands
message("Growing subislands",appendLF=F)
for (i in unique(fdata[,"island"])) {
	message(" ... ",appendLF=F)
	fdata[fdata[,"island"]==i,"subisland"] = findSmallIslands(
		fdata[fdata[,"island"]==i,"depth"],
		fdata[fdata[,"island"]==i,"position"],
		fdata[fdata[,"island"]==i,"outlier"]
		)
	message(i,appendLF=F)
}
message("")
# redefine islands
fdata[,"island"] = as.character(interaction(fdata[,"island"],fdata[,"subisland"]))
fdata[,"subisland"] = NULL

# output some stats
message(paste("   Total grown islands:",length(unique(fdata[,"island"])),"(",round(length(unique(fdata[,"island"]))/nrow(scaffolds)),"per chromosome",")"))

# Plot fully grown subislands
if (byplot) print(ggplot(fdata, aes(x=position, y=depth, shape=outlier, color=factor(cumsum(c(0,diff(as.numeric(factor(island)))!=0))%%7)))
	+ geom_point()
	+ facet_grid(.~scaffold, scales="free_x", space="free_x")
	+ theme_bw()
	+ theme(legend.position="none")
	+ ggtitle("Fully grown subislands"))

## Shift copy number regions
################################################################################################################################################################
################################################################################################################################################################

# calculate a neutral prediction from given predicted depths, islands and intercepts
#	neutrality can be either:
#	 * left-sided: based on intercept value only; that is the pdepth is rescaled so that the left intercept is 1
#  	 * right-sided: based on right-intercept value only; that is the pdepth is rescaled so that the right intercept is 1
# 	 * both: mean of left-sided and right-sided neutrals
neutralPrediction <- function(pdepth, islands, logintercepts, scaffolds., side, scaffoldAnnotation=scaffolds) {
	stopifnot(length(pdepth) == length(islands))
	stopifnot(length(pdepth)==length(scaffolds.))
	stopifnot(all(as.character(unique(islands)) %in% names(logintercepts)))
	stopifnot(all(as.character(unique(scaffolds.)) %in% rownames(scaffoldAnnotation)))
	stopifnot(side %in% c("left","right","both"))
	stopifnot(all(is.finite(pdepth)))
	stopifnot(all(is.finite(logintercepts)))
	stopifnot(all(scaffoldAnnotation[unique(scaffolds.),"end"] <= length(pdepth)))
	islands = as.character(islands)

	if (side == "left") {
		result = pdepth / 2^logintercepts[islands]
	} else if (side == "right") {
		leftneutral = neutralPrediction(pdepth, islands, logintercepts, scaffolds., side="left")
		result = leftneutral / leftneutral[scaffoldAnnotation[scaffolds.,"end"]]
	} else if (side == "both") {
		leftneutral = neutralPrediction(pdepth, islands, logintercepts, scaffolds., side="left")
		rightneutral = neutralPrediction(pdepth, islands, logintercepts, scaffolds., side="right")
		result = (leftneutral + rightneutral) / 2
	}

	return (result)
}

# compute properties of islands given a spline prediction
annotateIslands <- function(pdepth, islands, logintercepts, position=fdata[,"position"], depth=fdata[,"depth"], scaffolds.=fdata[,"scaffold"], scaffoldAnnotation=scaffolds) {
	n = length(position)
	stopifnot(n>0)
	stopifnot(n==length(depth))
	stopifnot(n==length(pdepth))
	stopifnot(n==length(islands))
	stopifnot(n==length(scaffolds.))
	stopifnot(all(as.character(unique(islands)) %in% as.character(names(logintercepts))))
	stopifnot(all(as.character(unique(scaffolds.)) %in% as.character(rownames(scaffoldAnnotation))))
	islands = as.character(islands)

	# create data frame
	result 			= data.frame(name=as.character(unique(islands)))
	result[,"name"] 	= as.character(result[,"name"])
	rownames(result)	= result[,"name"]

	# calculate basic properties
	result[,"start"] 	= sapply(result[,"name"], function(x){which.max(islands==x)})
	result[,"end"] 		= sapply(result[,"name"], function(x){n-which.max(rev(islands)==x)+1})
	result[,"scaffold"] = sapply(result[,"start"], function(x){scaffolds.[x]})
	result[,"size"] 	= sapply(result[,"name"], function(x){sum(islands==x)})

	# check if island has many low-coverage frames (which means it is probably a deletion)
	# 	that is if more than half of frames in the island has low coverage
	# result[,"probableDeletion"] = sapply(result[,"name"], function(x){sum(depth[islands==x]<basedepth/15)}) >= result[,"size"] / 2

	# calculate noise in the island region
	# 	meanNoise is the harmonic mean of the multiplicative deviation
	#	noisePercentage is the percentage of frames that need to be removed to get the meanNoise between a treshold value
	result[,"meanNoise"]= sapply(result[,"name"], function(x){
		islandIndices = islands == x
		islandIndices[c(FALSE,diff(as.numeric(factor(islands)))!=0)] = FALSE
		islandIndices[c(diff(as.numeric(factor(islands)))!=0,FALSE)] = FALSE
		# islandIndices[which.max(islandIndices)] = FALSE
		# islandIndices[length(islandIndices)-which.max(rev(islandIndices))+1] = FALSE
		harmonic.mean(depthquotient(pdepth[islandIndices],depth[islandIndices]))})
	# result[result[,"probableDeletion"],"noise"] = sapply(result[result[,"probableDeletion"],"name"], function(x){sum(depth[islands==x]>basedepth/15)/sum(islands==x)})
	noiseTreshold = 1+2*(median(result[islands,"meanNoise"])-1)
	result[,"noisePercentage"] = sapply(result[,"name"], function(x) {
		islandIndices = islands == x
		islandIndices[c(FALSE,diff(as.numeric(factor(islands)))!=0)] = FALSE
		islandIndices[c(diff(as.numeric(factor(islands)))!=0,FALSE)] = FALSE
		# islandIndices[which.max(islandIndices)] = FALSE
		# islandIndices[length(islandIndices)-which.max(rev(islandIndices))+1] = FALSE

		cumnoise = 1/(cumsum(1/sort(depthquotient(pdepth[islandIndices],depth[islandIndices])))/(1:sum(islandIndices)))
		cumindex = sum(cumnoise<noiseTreshold)
		return (1-cumindex/sum(islandIndices))
		})

	# calculate height properties
	result[,"intercept"] = logintercepts[result[,"name"]]
	# result[,"rightIntercept"] = result[,"intercept"] - log2(pdepth[scaffoldAnnotation[result[,"scaffold"],"start"]]) + log2(pdepth[scaffoldAnnotation[result[,"scaffold"],"end"]])
	leftneutralprediction = neutralPrediction(pdepth, islands, logintercepts, scaffolds., side="left", scaffoldAnnotation=scaffoldAnnotation)
	result[,"rightIntercept"] = result[,"intercept"] + log2(leftneutralprediction[scaffoldAnnotation[result[,"scaffold"],"end"]])
	result[,"doubleIntercept"] = log2(2^result[,"intercept"] + 2^result[,"rightIntercept"])-1
	result[,"ratio"] = 2^(result[,"doubleIntercept"] - neutralIntercept(result, side="both", islands=islands))

	# Also include ratio compared to "neutral" island in scaffold, if it exists
	result[,"ratioToNeutral"] = result[,"ratio"] / sapply(result[,"scaffold"], function(x) {
			neutr = result[,"scaffold"]==x & grepl("Neutral",result[,"name"])
			stopifnot(sum(neutr)<=1)
			return (if (sum(neutr)==0) 1 else result[neutr,"ratio"])
		})

	result[,"ratioToMedian"] = result[,"ratio"] / sapply(result[,"scaffold"], function(x)median(result[as.character(islands[scaffolds.==x]),"ratio"]))

	return (result)
}

# compute neutral value of intercept
#	this value is the median of the genome's "intercept" (left, right ot double-sided) value for the prediction
neutralIntercept <- function(islandAnnotation, side, islands=fdata[,"island"], method="quantile") {
	stopifnot(method %in% c("mean","median","quantile"))
	stopifnot(side %in% c("left","right","both"))

	# read method
	methodfun = if (method=="quantile") (function(x)quantile(x,.70)) else {if (method=="mean") mean else median}


	# calculate intercept
	if (side == "left") {
		result = log2(methodfun(2^(islandAnnotation[islands,"intercept"])))
	} else if (side == "right") {
		result = log2(methodfun(2^(islandAnnotation[islands,"rightIntercept"])))
	} else if (side == "both") {
		result = log2(methodfun(2^(islandAnnotation[islands,"doubleIntercept"])))
	}

	return (result)
}

# predict spline baseline on a single scaffold
scaffoldspline <- function(depth, position, islands, outliers, scaflen, forceLine.=forceLine) {
		stopifnot(all(is.finite(outliers)))
		stopifnot(all(as.numeric(outliers)%in%c(0,1)))

		# Fit spline on logs
		M = if (forceLine.) 1 else {if (scaflen>400000) 3 else 1}
		lambda = if (forceLine.) 10^-3 else (if (scaflen>400000) 10^-6.5 else 10^-2)
		fit = 2^superspline(
			x=position,
			y=log2depth(depth),
			island=islands,
			M=M,
			lambda=lambda,
			weights = pmax(2^minlogdepth,depth) * (1 - (1-outlierWeight)*as.numeric(outliers))
		)

		# Correct fitted spline height
		if (length(unique(islands))>1) {
			depthcorrection = lm(log2(pmax(2^minlogdepth,depth)/pmax(2^minlogdepth,fit)) ~ factor(islands) - 1, weight=(1 - (1-outlierWeight)*as.numeric(outliers)))
			result = fit * 2^predict(depthcorrection)
			stopifnot(length(depthcorrection$coefficients) == length(attr(fit,"islands")))
			attr(result,"islands") = attr(fit,"islands") + (depthcorrection$coefficients)[paste("islands",names(attr(fit,"islands")),sep="")]
			names(attr(result,"islands")) = names(attr(fit,"islands"))
		} else {
			depthcorrection = lm(log2(pmax(2^minlogdepth,depth)/pmax(2^minlogdepth,fit)) ~ 1, weight=(1 - (1-outlierWeight)*as.numeric(outliers)))
			newscale = mean(log2(pmax(2^minlogdepth,depth)/pmax(2^minlogdepth,fit)))
			result = fit * 2^newscale
			attr(result,"islands") = attr(fit,"islands") + newscale
			names(attr(result,"islands")) = as.character(islands[1])
		}
		return (fit)
}

# predict spline baseline curve across multiple chromosomes
genomespline <- function(depth=fdata[,"depth"], position=fdata[,"position"], scaffolds.=fdata[,"scaffold"], islands=fdata[,"island"], outliers=fdata[,"outlier"], scaffoldAnnotation=scaffolds, forceLine.=forceLine) {
	message("Fitting spline on whole genome",appendLF=F)
	stopifnot(length(depth)==length(position))
	stopifnot(length(depth)==length(scaffolds.))
	stopifnot(length(depth)==length(islands))
	stopifnot(length(depth)==length(outliers))
	stopifnot(all(unique(as.character(scaffolds.)) %in% as.character(rownames(scaffoldAnnotation))))

	# define initial results
	logintercepts = c()
	pdepth = rep(NA,length(depth))

	# for each scaffold, calculate prediction and add to result
	for (scaf in unique(scaffolds.)) {
		message(" ... ",appendLF=F)
		fit = scaffoldspline(depth[scaffolds.==scaf], position[scaffolds.==scaf], islands[scaffolds.==scaf], outliers[scaffolds.==scaf], scaffoldAnnotation[scaf,"length"], forceLine.)
		pdepth[scaffolds.==scaf] = fit
		logintercepts[names(attr(fit,"islands"))] = attr(fit,"islands")
		message(scaf,appendLF=F)
		}

	attr(pdepth,"logintercepts") = logintercepts
	message("")
	return (pdepth)
}

# expand a region, using the HMM and returning a new region vector
expandIsland <- function(region, intercept, outlier, depth) {
	stopifnot(length(region)>0)
	stopifnot(length(region)==length(outlier))
	stopifnot(length(region)==length(depth))
	stopifnot(all(is.finite(outlier)))
	stopifnot(all(is.finite(depth)))
	stopifnot(is.finite(intercept))

	result = rep(FALSE,length(region))
	prediction = rep(1,length(region))

	# if region is empty, ignore
	if (sum(region)==1) {

	# if region is not monolythic, repeat procedure in each sub-region
	} else if (sum(diff(region)==1) > 1) {
	subparts = cumsum(c(region[1],diff(region)!=0))
	for (i in unique(subparts)) {
		if (!(i%%2)) next
		result = result | expandIsland(subparts==i, intercept, outlier, depth)
		prediction = pmax(1, superspline(
			x = (1:length(region))[region],
			y = depth[region],
			island = NULL,
			weights = (1 - (1-outlierWeight)*outlier[region]),
			newx = 1:length(region),
			newisland = NULL,
			M=1, lambda=10^-2
			))
	}

	# else expand using HMM
	} else {
		prediction = pmax(1, superspline(
			x = (1:length(region))[region],
			y = depth[region],
			island = NULL,
			weights = (1 - (1-outlierWeight)*outlier[region]),
			newx = 1:length(region),
			newisland = NULL,
			M=1, lambda=10^-3
			))

		# prediction = neutralprediction * 2^intercept
		viterbiState = CNVHMM(depth, prediction, delta=1.2, p=0.2)

		margin = floor(5000/frameSize)
		islandstart = NA
		islandend = NA
		newislandend = NA
		newislandstart = NA

		if (sum(region & viterbiState=="CNV1")>0) {
			islandstart = which.max(region)
			islandend 	= length(region) - which.max(rev(region)) + 1

			result[region] = TRUE			
			result[1:length(region)>=islandend & cumsum(1:length(region)>=islandend & viterbiState!="CNV1")==0] = TRUE
			result[1:length(region)<=islandstart & rev(cumsum(rev(1:length(region)<=islandstart & viterbiState!="CNV1")))==0] = TRUE

			# newislandstart = which.max((1:length(region)) < islandend & viterbiState=="CNV1")
			# newislandend = length(region) - which.max(rev((1:length(region)) > islandstart & viterbiState=="CNV1")) + 1
			# endsearchstart = max(islandend-margin, islandstart)
			# startsearchstart = min(islandend+margin, islandend)
			# newislandend = (endsearchstart) + which.min((viterbiState=="CNV1")[-(1:(endsearchstart))]) - 1
			# newislandstart = (startsearchstart) - which.min(rev((viterbiState=="CNV1")[1:(startsearchstart)])) + 2

			# if (newislandend > newislandstart) {
			# 	result[endsearchstart:newislandend] = TRUE
			# 	result[newislandstart:startsearchstart] = TRUE
			# }
		}

		# # plot inferred islands
		# 	print(ggplot(data.frame(region=region, intercept=intercept, depth=depth, x=1:length(depth), viterbiState=viterbiState, prediction=prediction, result=result), aes(x=x, y=depth, color=viterbiState, alpha=result))
		# 	+ geom_point(aes(shape=region))
		# 	+ geom_line(aes(y=prediction), color="red")
		# 	+ geom_line(aes(y=prediction*delta), color="red", linetype=2, alpha=0.5)
		# 	+ geom_line(aes(y=prediction/delta), color="red", linetype=2, alpha=0.5)
		# 	)
	}

	attr(result,"prediction") = prediction
	return (result)
}

updateBaseline <- function() {
	pdepth = genomespline()
	fdata[,"pdepth"] <<- pdepth
	fdata[,"neutralpred"] <<- neutralPrediction(pdepth, fdata[,"island"], attr(pdepth,"logintercepts"), fdata[,"scaffold"], side="both")

	return (annotateIslands(pdepth, fdata[,"island"], attr(pdepth, "logintercepts")))
}

			# ## Merge adjacent similar islands
			# ################################################################################################################################################################
			# ################################################################################################################################################################
			# adjacent <- function(island1, island2, islandAnnotation=islands) {
			# 	return (as.numeric(islandAnnotation[island1,"scaffold"] == islandAnnotation[island2,"scaffold"] &
			# 		((islandAnnotation[island1,"start"]==islandAnnotation[island2,"end"]+1) | (islandAnnotation[island2,"start"]==islandAnnotation[island1,"end"]+1))))
			# }

			# islands = updateBaseline()

			# # find most similar adjacent islands and merge
			# #	repeat until all islands are too dissimilar to merge
			# message("Merging islands ... ",appendLF=F)
			# distances = outer(as.character(islands[,"name"]),as.character(islands[,"name"]), function(x,y){abs(islands[x,"ratio"]-islands[y,"ratio"]) + (100+ratioMergeTreshold)*(!adjacent(x,y))})
			# stopifnot(all(is.finite(distances)))
			# while (min(distances) <= ratioMergeTreshold) {
			# 	# find islands to merge
			# 	mergeIndices = which(distances==min(distances),arr.ind=T)[1,]
			# 	island1 = islands[mergeIndices,"name"][which.min(abs(islands[mergeIndices,"ratio"]-1))]
			# 	island2 = rev(islands[mergeIndices,"name"])[which.max(abs(rev(islands[mergeIndices,"ratio"])-1))]

			# 	# print message
			# 	message(paste(island2," -> ",island1," ... ",sep=""),appendLF=F)

			# 	# recalculate spline on this scaffold
			# 	fdata[fdata[,"island"]==island2,"island"] = island1
			# 	newpdepth = scaffoldspline(
			# 		fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"depth"],
			# 		fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"position"],
			# 		fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"island"],
			# 		fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"outlier"],
			# 		scaffolds[islands[island1,"scaffold"],"length"]
			# 		)

			# 	# redefine islands
			# 	fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"pdepth"] = newpdepth
			# 	logintercepts = islands[,"intercept"]
			# 	names(logintercepts) = islands[,"name"]
			# 	logintercepts[names(attr(newpdepth,"islands"))] = attr(newpdepth,"islands")
			# 	islands = annotateIslands(fdata[,"pdepth"], fdata[,"island"], logintercepts)

			# 	# redefine distance matrix for next iteration
			# 	distances = outer(as.character(islands[,"name"]),as.character(islands[,"name"]), function(x,y){abs(islands[x,"ratio"]-islands[y,"ratio"]) + (100+ratioMergeTreshold)*(!adjacent(x,y))})
			# 	stopifnot(all(is.finite(distances)))
			# }
			# message("Done.")

# expand islands
#	create matrix of expanded regions
#	islands where each frame can be allocated to another island are "absorbed" into those other islands
# 	other frames allocated to multiple islands are kept in favor of their original region, or they are allocated randomly
islands = updateBaseline()
message(paste("   Neutral intercept:",neutralIntercept(islands, side="both"), "(", 2^neutralIntercept(islands, side="both"), ")"))
message("Shifting islands ...")

# define second island set
fdata[,"island2"] = fdata[,"island"]

shifts = 0
repeat {
	shifts = shifts+1

	# try to extend within each scaffold
	for (scaf in scaffolds[,"name"]) {

		# if there is only one island, there is no shifting to do
		stopifnot(sum(islands[,"scaffold"]==scaf)>=1)
		if (sum(islands[,"scaffold"]==scaf)==1) next

		# for each scaffold, compute the extensions of all the islands
		newAllocationMatrix = c()
		predictionMatrix = c()
		for (i in islands[islands[,"scaffold"]==scaf,"name"]) {
			expansion = expandIsland(
				fdata[fdata[,"scaffold"]==scaf,"island2"]==i,
				islands[i,"intercept"],
				fdata[fdata[,"scaffold"]==scaf,"outlier"],
				fdata[fdata[,"scaffold"]==scaf,"depth"]
				)
			newAllocationMatrix = cbind(newAllocationMatrix, expansion)
			predictionMatrix = cbind(predictionMatrix, attr(expansion,"prediction"))
		}
		colnames(newAllocationMatrix) = as.character(islands[islands[,"scaffold"]==scaf,"name"])
		colnames(predictionMatrix) = as.character(islands[islands[,"scaffold"]==scaf,"name"])

		# # shift border frames that belong better to the next island
		# unanimousIndices = apply(newAllocationMatrix,1,sum)==1

		# if (sum(unanimousIndices)>1) {
		# 	unanimousIslands = apply(newAllocationMatrix[unanimousIndices,],1,which.max)
		# 	unanimousIslands = colnames(newAllocationMatrix)[unanimousIslands]

		# 	fdata[fdata[,"scaffold"]==scaf,][unanimousIndices,"island2"] = unanimousIslands
		# }

		# merge islands that can be fully absorbed
		for (i in as.character(islands[islands[,"scaffold"]==scaf,"name"])) {
			otherColumns = !(colnames(newAllocationMatrix)==i)
			if (sum(fdata[,"island2"]==i)==0) next
			if (sum(otherColumns)<=1) next
			if (all(rowSums(data.matrix(newAllocationMatrix[fdata[fdata[,"scaffold"]==scaf,"island2"]==i, otherColumns]))>=1)) {
				newIslands = apply(newAllocationMatrix[fdata[fdata[,"scaffold"]==scaf,"island2"]==i, otherColumns],1,which.max)
				fdata[fdata[,"island2"]==i,"island2"] = colnames(newAllocationMatrix[,otherColumns])[newIslands]
			}
		}

		newAllocationMatrix = newAllocationMatrix[,unique(fdata[fdata[,"scaffold"]==scaf,"island2"])]
		predicationMatrix = predictionMatrix[,unique(fdata[fdata[,"scaffold"]==scaf,"island2"])]

		# shift border frames in contested region to the optimal location
		contestedIndices = apply(newAllocationMatrix,1,sum)==2
		leftContestedIslands = apply(newAllocationMatrix,1,which.max)
		rightContestedIslands = apply(newAllocationMatrix,1,function(x)length(x)-which.max(rev(x))+1)
		parts = factor(as.numeric(interaction(leftContestedIslands,rightContestedIslands))*contestedIndices)

			# print(ggplot(fdata[fdata[,"scaffold"]==scaf,], aes(x=position, y=depth, color=parts))
			# 	+ geom_point()
			# )

		# print(newAllocationMatrix)

		for (part in unique(parts)) {
			if (!all(contestedIndices[parts==part])) next
			stopifnot(length(unique(leftContestedIslands[parts==part]))==1)
			stopifnot(length(unique(rightContestedIslands[parts==part]))==1)
			leftIsland = leftContestedIslands[parts==part][1]
			rightIsland = rightContestedIslands[parts==part][1]

			observation = fdata[fdata[,"scaffold"]==scaf,"depth"][parts==part]
			leftPrediction = predictionMatrix[parts==part,leftIsland]
			rightPrediction = predictionMatrix[parts==part,rightIsland]

			# penalty for fit error
			leftError = depthquotient(leftPrediction, observation)
			rightError = depthquotient(rightPrediction, observation)

			# penalty for fit unsmoothness
			# if (length(observation)<10) {
			# 	leftUnsmoothness = rep(0,length(observation))
			# 	rightUnsmoothness = rep(0,length(observation))
			# } else {
				# leftUnsmoothness = c(0,diff(filter(leftPrediction,rep(1/5,5),sides=2), lag=3),0)
				# rightUnsmoothness = c(0,diff(filter(rightPrediction,rep(1/5,5),sides=2), lag=3),0)
				# leftUnsmoothness[!is.finite(leftUnsmoothness)] = 0
				# rightUnsmoothness[!is.finite(rightUnsmoothness)] = 0
				# leftUnsmoothness = (c(leftUnsmoothness[1],leftUnsmoothness)+c(leftUnsmoothness,rev(leftUnsmoothness)[1]))/2
				# rightUnsmoothness = (c(rightUnsmoothness[1],rightUnsmoothness)+c(rightUnsmoothness,rev(rightUnsmoothness)[1]))/2
			# }
			leftUnsmoothness = abs(diff(leftPrediction))
			rightUnsmoothness = abs(diff(rightPrediction))
			leftUnsmoothness = (c(leftUnsmoothness[1],leftUnsmoothness)+c(leftUnsmoothness,rev(leftUnsmoothness)[1]))/2
			rightUnsmoothness = (c(rightUnsmoothness[1],rightUnsmoothness)+c(rightUnsmoothness,rev(rightUnsmoothness)[1]))/2

			# penalty for island assignment inertia
			leftInertia = (as.character(fdata[fdata[,"scaffold"]==scaf,"island2"][parts==part])!=colnames(newAllocationMatrix)[leftIsland])
			rightInertia = (as.character(fdata[fdata[,"scaffold"]==scaf,"island2"][parts==part])!=colnames(newAllocationMatrix)[rightIsland])

			# total penalty
			lambda = mean(c(abs(leftError-observation),abs(rightError-observation)))/mean(c(leftUnsmoothness,rightUnsmoothness))*5
			# inertia = max(0, 0.01*shifts - 0.03) + mean(c(leftUnsmoothness[leftUnsmoothness!=0],rightUnsmoothness[rightUnsmoothness!=0]))
			inertia = max(lambda*c(leftUnsmoothness, rightUnsmoothness))
			if (!is.finite(lambda)) { lambda=0; inertia = max(0,0.01*shifts - 0.03) }
			leftPenalty = leftError #+ lambda*leftUnsmoothness + inertia*leftInertia
			rightPenalty = rightError #+ lambda*rightUnsmoothness + inertia*rightInertia

			triangle = outer(0:length(observation), 1:length(observation), "<")
			totalPenalties = triangle%*%rightPenalty + (1-triangle)%*%leftPenalty

			index = which.min(totalPenalties)-1
			fdata[fdata[,"scaffold"]==scaf,][parts==part,"island2"] = c(rep(colnames(newAllocationMatrix)[leftIsland],index), rep(colnames(newAllocationMatrix)[rightIsland], length(observation)-index))
		}
	}

	if (all(fdata[,"island"]==fdata[,"island2"])) break

	message(paste(shifts,".   frames shifted:",sum(fdata[,"island"]!=fdata[,"island2"])))
	message(paste("   islands absorbed:",
		if (length(unique(fdata[,"island2"]))==length(unique(fdata[,"island"]))) "none" else paste(unique(fdata[,"island"])[!(unique(fdata[,"island"]) %in% unique(fdata[,"island2"]))], collapse=", ")
		))

	# plot island shifts
	if (byplot) maskfun(print(maskfun(ggplot(fdata, aes(x=position, y=depth, alpha=island2!=island, shape=island2!=island, color=factor(cumsum(c(0,diff(as.numeric(factor(island2)))!=0))%%7)))
				+ geom_point()
				+ facet_grid(.~scaffold, scales="free_x", space="free_x")
				+ theme(panel.margin=unit(0,"lines"))
				+ theme_bw()
				+ scale_alpha_discrete(range=c(0.4,1))
				+ theme(legend.position="none")
				+ ggtitle("island shifts"))))

	# redefine islands
	fdata[,"island"] = fdata[,"island2"]
	islands = islands[unique(fdata[,"island2"]),]
	if (shifts > 25) break;

}
islands = updateBaseline()
fdata[,"island2"] = NULL

## Update outliers
################################################################################################################################################################
################################################################################################################################################################
# Update outlier annotation
#	new outliers are frames with:
#		- high local standard deviation, but not only counted within an island (removing natural high standard deviation in breakpoints)
#		- high amount of louche frames
#		- breakpoint frames (as they often have intermediary values, belonging half to the left and half to the right copy number region)
# 	very low depth frames in predominantly low-depth frames are now explicitly not called outliers, as they are probably deletions
message("Updating outliers ...")
# for (i in unique(fdata[,"island"])) {
# 	normaobs = fdata[fdata[,"island"]==i,"depth"] / fdata[fdata[,"island"]==i,"pdepth"]
# 	fdata[fdata[,"island"]==i,"deviation2"] = rollapply(normaobs, floor(12500/frameSize), sd, partial=T)
# 	fdata[!is.finite(fdata[,"deviation2"]),"deviation2"] = Inf
# }
fdata[,"deviation2"] = rollapply(fdata[,"depth"]/(fdata[,"neutralpred"]*islands[fdata[,"island"],"ratio"]), floor(12500/frameSize),sd,partial=T)
medianLocalSd2 = median(fdata[,"deviation2"])
message(paste("   Median deviation2:",medianLocalSd2))
stopifnot(is.finite(medianLocalSd2))
# fdata[,"outlier2"] = fdata[,"deviation2"]>1.5*medianLocalSd2 |
# 						# fdata[,"depth"] > max(c(2*basedepth,1.5*depthmax)) |
# 						# fdata[,"louchefrac"] < maxlouche |
# 						c(FALSE,diff(as.numeric(factor(fdata[,"island"])))!=0) | c(diff(as.numeric(factor(fdata[,"island"])))!=0,FALSE)
# 						# fdata[,"telomere"]
fdata[,"outlier2"] = fdata[,"outlier"]

# fdata[fdata[,"depth"]<basedepth/4 && islands[fdata[,"island"],"probableDeletion"],"outlier2"] = FALSE

# plot new outliers
if (byplot) maskfun(print(maskfun(ggplot(fdata, aes(x=position, y=deviation2))
		+ geom_point()
		+ geom_hline(aes(yintercept=2*medianLocalSd2), color="red")
		+ facet_grid(.~scaffold, scales="free_x", space="free_x")
		+ theme_bw()
		+ ggtitle("Local deviation and outlier treshold, second"))))

# if (byplot) maskfun(print(maskfun(ggplot(fdata, aes(x=position, y=depth))
# 			+ geom_point(aes(color=outlier2, alpha=outlier))
# 			+ geom_line(aes(y=pdepth), color="red")
# 			+ facet_grid(.~scaffold, scales="free_x", space="free_x")
# 			+ theme(panel.margin=unit(0,"lines"))
# 			+ theme_bw()
# 			+ scale_alpha_discrete(range=c(0.4,1))
# 			+ ggtitle("outlier update"))))

# redefine outliers
stopifnot(all(is.finite(fdata[,"outlier2"])))
fdata[,"outlier"] = fdata[,"outlier2"]
fdata[,"outlier2"] = NULL

## Merge adjacent similar islands
################################################################################################################################################################
################################################################################################################################################################
adjacencyTreshold = 20000
adjacent <- function(island1, island2, islandAnnotation=islands) {
	# i1 = c(island1,island2)[which.min(c(islandAnnotation[island1,"start"], islandAnnotation[island2,"start"]))]
	# i2 = c(island1,island2)[which.max(c(islandAnnotation[island1,"start"], islandAnnotation[island2,"start"]))]
	s1 = islandAnnotation[island1,"start"]
	s2 = islandAnnotation[island2,"start"]
	e1 = islandAnnotation[island1,"end"]
	e2 = islandAnnotation[island2,"end"]
	return (as.numeric(
		island1 != island2 &
		islandAnnotation[island1,"scaffold"] == islandAnnotation[island2,"scaffold"] &
		((s1 > s2 & s1-e2 < adjacencyTreshold/frameStep ) |
		 (s2 > s1 & s2-e1 < adjacencyTreshold/frameStep))
		# (s2-e1 < adjacencyTreshold/frameStep | s1-e2 < adjacencyTreshold/frameStep | (!all(sort(s1,e1,s2)==c(s1,e))))
		# (islandAnnotation[island2,"start"] - islandAnnotation[island1,c("end")] < adjacencyTreshold/frameStep |
		# 	islandAnnotation[island1,"start"] - islandAnnotation[island2,c("end")] < adjacencyTreshold/frameStep |
		# 	(islandAnnotation[island1,"end"] < islandAnnotation[island2,"end"] & islandAnnotation[island1,"end"] > islandAnnotation[island1,"start"]) |
		# 	(islandAnnotation[island2,"end"] < islandAnnotation[island1,"end"] & islandAnnotation[island2,"end"] > islandAnnotation[island2,"start"])
		# 	)
		# min(c(islandAnnotation[island1,"start"]-islandAnnotation[island2,"end"],islandAnnotation[island2,"start"]-islandAnnotation[island1,"end"])) < adjacencyTreshold/frameStep
		))
		# ((islandAnnotation[island1,"start"]==islandAnnotation[island2,"end"]+1) | (islandAnnotation[island2,"start"]==islandAnnotation[island1,"end"]+1))))
}

# find most similar adjacent islands and merge
#	repeat until all islands are too dissimilar to merge
message("Merging islands ... ",appendLF=F)
distances = outer(as.character(islands[,"name"]),as.character(islands[,"name"]), function(x,y){abs(islands[x,"ratio"]-islands[y,"ratio"]) + (100+ratioMergeTreshold)*(!adjacent(x,y))})
stopifnot(all(is.finite(distances)))
while (min(distances) <= ratioMergeTreshold) {
	# find islands to merge
	mergeIndices = which(distances==min(distances),arr.ind=T)[1,]
	island1 = islands[mergeIndices,"name"][which.min(abs(islands[mergeIndices,"ratio"]-1))]
	island2 = rev(islands[mergeIndices,"name"])[which.max(abs(rev(islands[mergeIndices,"ratio"])-1))]

	# print message
	message(paste(island2," -> ",island1," ... ",sep=""),appendLF=F)

	# move islands
	fdata[fdata[,"island"]==island2,"island"] = island1

	# recalculate spline
	newpdepth = scaffoldspline(
		fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"depth"],
		fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"position"],
		fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"island"],
		fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"outlier"],
		scaffolds[islands[island1,"scaffold"],"length"]
		)

	# redefine islands
	fdata[fdata[,"scaffold"]==islands[island1,"scaffold"],"pdepth"] = newpdepth
	logintercepts = islands[,"intercept"]
	names(logintercepts) = islands[,"name"]
	logintercepts[names(attr(newpdepth,"islands"))] = attr(newpdepth,"islands")
	islands = annotateIslands(fdata[,"pdepth"], fdata[,"island"], logintercepts)

	# redefine distance matrix for next iteration
	distances = outer(as.character(islands[,"name"]),as.character(islands[,"name"]), function(x,y){abs(islands[x,"ratio"]-islands[y,"ratio"]) + (100+ratioMergeTreshold)*(!adjacent(x,y))})
	stopifnot(all(is.finite(distances)))
}
message("Done.")

## Find and eliminate small islands
################################################################################################################################################################
################################################################################################################################################################
maxCNRsize = 5000
message("Cleaning islands ... ",appendLF=F)
logintercepts = islands[,"intercept"]
names(logintercepts) = islands[,"name"]
fdata[,"neutralpred"] = fdata[,"pdepth"] / islands[fdata[,"island"],"ratio"]
fdata[,"island2"] = fdata[,"island"]
for (i in islands[,"name"]) {
	# If the island is small, redistribute its frames among neighbouring islands so that the total distance between predicted and observed values is minimal
	if (islands[i,"size"]*frameSize < maxCNRsize) {
		message(paste(i,"..."),appendLF=F)
		# stop if small island is the only island in scaffold (should not happen)
		stopifnot(!((islands[i,"start"]==1 || fdata[islands[i,"start"]-1,"scaffold"] != islands[i,"scaffold"]) &&
			(islands[i,"end"]==nrow(fdata) || fdata[islands[i,"end"]+1,"scaffold"] != islands[i,"scaffold"])))
		# if the island is at the beginning of a scaffold, join to right island
		if (islands[i,"start"]==1 || (fdata[islands[i,"start"]-1,"scaffold"]) != islands[i,"scaffold"]) {
			rightIsland = fdata[islands[i,"end"]+1,"island2"]
			fdata[fdata[,"island2"]==i,"island2"] = rightIsland
		# if the island is at the end of a scaffold, join to the left island
		} else if (islands[i,"end"]==nrow(fdata) || (fdata[islands[i,"end"]+1,"scaffold"]) != islands[i,"scaffold"]) {
			leftIsland = fdata[islands[i,"start"]-1,"island2"]
			fdata[fdata[,"island2"]==i,"island2"] = leftIsland
		# otherwise, find optimal distribution
		} else {
			leftIsland = fdata[islands[i,"start"]-1,"island2"]
			rightIsland = fdata[islands[i,"end"]+1,"island2"]
			leftPrediction = pmax(10,fdata[fdata[,"island2"]==i,"neutralpred"] * islands[leftIsland,"ratio"])
			rightPrediction = pmax(10,fdata[fdata[,"island2"]==i,"neutralpred"] * islands[rightIsland,"ratio"])
			observation = pmax(10,fdata[fdata[,"island2"]==i,"depth"])

			# calculate penalty comparing observation to left and right islands
			leftPenalty = depthquotient(leftPrediction, observation)
			rightPenalty = depthquotient(rightPrediction, observation)

			# print(leftPenalty)
			# print(rightPenalty)

			# calculate actual penalty for each possible location of the splitting point
			triangle = outer(0:length(observation), 1:length(observation), "<")
			totalPenalties = triangle%*%rightPenalty + (1-triangle)%*%leftPenalty

			# find and implement the optimal splitting point
			index = which.min(totalPenalties)-1
			fdata[fdata[,"island2"]==i,"island2"] = c(rep(leftIsland,index), rep(rightIsland, length(observation)-index))
		}
	}

	# island is noisy
}
message("")

# redefine islands
fdata[,"island"] = fdata[,"island2"]
fdata[,"island2"] = NULL

# update baseline fit
islands = updateBaseline()


# ## Neutralize islands
# ################################################################################################################################################################
# ################################################################################################################################################################
# Find all islands that are close to the genome-wide baseline and annotate them as "neutral"
#	this will merge 
message("Neutralizing islands ...")
neutralizationTreshold = 1.15
neutralIslands = islands[abs(islands[,"ratio"]-1)<neutralizationTreshold-1,"name"]
if (length(neutralIslands) == 0) {
	message("   islands neutralized: none")
} else {
	fdata[fdata[,"island"]%in%neutralIslands,"island"] = paste(fdata[fdata[,"island"]%in%neutralIslands,"scaffold"],"Neutral",sep=".")
	pdepth = genomespline()
	fdata[,"pdepth"] = pdepth
	islands2 = annotateIslands(pdepth, fdata[,"island"], attr(pdepth,"logintercepts"))
	message(paste("   islands neutralized:",
		paste(islands[,"name"][!(islands[,"name"] %in% islands2[,"name"])], collapse=", ")
		))
	islands=islands2
}

# annotate final noisy islands
noiseTreshold = 1-(1-median(islands[fdata[,"island"],"meanNoise"]))*2
message(paste("   Noise treshold:",noiseTreshold))

logintercepts = islands[,"intercept"]
names(logintercepts) = islands[,"name"]
fdata[,"neutralpred"] = fdata[,"pdepth"] / islands[fdata[,"island"],"ratio"]
showSign <- function(x){paste(gsub("FALSE","",gsub("TRUE","+",as.character(x>0))),x,sep="")}
# plot new islands
scaffolds[,"scaffold"] = scaffolds[,"name"]
if (byplot) maskfun(print(maskfun(ggplot(fdata, aes(x=position, y=depth))
			+ geom_vline(data=scaffolds, aes(xintercept=centromere))
			+ geom_point(aes(shape=outlier, color=factor(cumsum(c(0,diff(as.numeric(factor(island)))!=0))%%7)))
			+ geom_line(aes(y=pdepth), color="red")
			+ geom_line(aes(y=pdepth+flatBias), color="red", linetype=2, alpha=0.5)
			+ geom_line(aes(y=pdepth-flatBias), color="red", linetype=2, alpha=0.5)
			+ geom_line(aes(y=neutralpred),color="green")
			+ geom_hline(yintercept=2^neutralIntercept(islands, "both"), color="black", linetype=2, alpha=0.5)
			+ geom_text(data=islands, aes(
				x=fdata[start,"position"],
				y=fdata[start,"pdepth"],
				# label=paste(name,paste(showSign(round(100*ratio)-100),"%",sep=""),paste(round(noise2*100),"%",sep=""),sep="\n")), hjust=0, vjust=0, size=1.5, color="black")
				label=paste(name,paste(showSign(round(100*ratio)-100),"%",sep=""),islands[,"meanNoise"], paste(round(islands[,"noisePercentage"]*100),"%",sep=""),sep="\n")),
				hjust=0, vjust=0, size=1.5, color="black")
			+ facet_grid(.~scaffold, scales="free_x", space="free_x")
			+ theme(panel.margin=unit(0,"lines"))
			+ theme_bw()
			+ scale_alpha_discrete(range=c(0.4,1))
			+ ggtitle("islands after merging")
			+ theme(legend.position="none"))))

# plot island height distribution
binsize=0.1
if (byplot) print(ggplot(islands, aes(x=ratio))
	+ geom_bar(aes(weight=size, y=..density.., fill=scaffold), color="black", binwidth=binsize, stat="bin")
	+ geom_density(data=fdata, aes(x=islands[island,"ratio"]), fill="black", alpha=0.5)
	+ scale_x_continuous(breaks=c(seq(0,2,0.25),seq(0,2,1/3)))
	+ ggtitle("island height distribution")
	+ theme_bw()
	)

binsize=0.1
if (byplot) print(ggplot(islands, aes(x=ratioToNeutral))
	+ geom_bar(aes(weight=size, y=..density.., fill=scaffold), color="black", binwidth=binsize, stat="bin")
	+ geom_density(data=fdata, aes(x=islands[island,"ratio"]), fill="black", alpha=0.5)
	+ scale_x_continuous(breaks=c(seq(0,2,0.25),seq(0,2,1/3)))
	+ ggtitle("island height distribution")
	+ theme_bw()
	)

# Stop plot
if (byplot) graphics.off()

# If required, output fdata table file
if (writeTab) write.table(x=fdata, file=tabFile, quote=F, sep="\t", row.names=F, col.names=T)

# If required, output islands table file
if (writeIslands) write.table(x=islands, file=islandsFile, quote=F, sep="\t", row.names=T, col.names=T)


endTime = Sys.time()

message(paste("Time needed:",as.numeric(difftime(endTime,startTime,units="mins")),"minutes"))
