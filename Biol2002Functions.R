#R Functions from Biology 2002/3004. Written by Elise Morten with comments and edits by J. Teshera-Levye

#Libraries to import
library(plyr)
library(car)
library(lme4)
library(ape)
library(vegan)
library(reshape2)
#*************
#New functions by J. Teshera-Levye
#

#Reading an OTU into a data frame from text file; pretty superfluous but exists as convenience function
#Wrapped in a tryCatch to produce friendlier errors.
read.OTUTable <- function(filename){
	out <- tryCatch({read.table(filename, comment='', header = T,
						sep="\t", as.is = TRUE, skip = 1, check.names = F)},
					error=function(cond){
						message("You got an error message. Here's the original error: ")
						message(cond)
						message("Check your input. Is everything spelled correctly? Are you in the right directory?")
						return(NULL)
					},
					warning=function(cond){
						message("You got a warning message.")
						message(cond)
						message("Check your OTU table file and try again.")
						return(NULL)
						})
	return(out)
	}


#Reading in a meta data table; again, superfluous but clearner
read.metadata <- function(filename){
	out <- tryCatch({read.table(filename, header=T, sep="\t", as.is=TRUE, check.names=F)},
				error=function(cond){
						message("You got an error message. Here's the original error: ")
						message(cond)
						message("Check your input. Is everything spelled correctly? Are you in the right directory?")
						return(NULL)
					},
					warning=function(cond){
						message("You got a warning message.")
						message(cond)
						message("Check your OTU table file and try again.")
						return(NULL)
						})
	return(out)
}

#Alternate alpha diversity data frame organizer - takes in alpha diversity file and metadata, returns single data frame with everything, which can then be used to plot, aggregate, summarize, however students want
#Takes in either data frames or strings as file paths
#First columns need to be sample ID, even if they aren't called that
#Prints messages if messages is T
alpha_organize <- function(alpha_data, meta_data, messages=T) {
	#Reads in alpha diversity data from qiime output
	if(is.data.frame(alpha_data)){
		alpha <- alpha_data
	} else if (is.character(alpha_data)) {
		alpha <- read.table(alpha_data, sep = "\t", header = TRUE, as.is = TRUE)
	} else {
		message("Alpha diversity of type not supported. Please use a file path or data frame and try again.")
		return(NULL)
	}

	metrics <- names(alpha)[-1]
	names(alpha)[1] <- SampleID

	#Reads in or assigns metadata
	if(is.data.frame(meta_data)){
		metadata <- meta_data
	} else if (is.character(meta_data)) {
		metadata <- read.table(meta_data, sep = "\t", header = TRUE, as.is = TRUE)
	} else {
		message("Meta data table of type not supported. Please use a file path or data frame and try again.")
		return(NULL)
	}

	if(names(metadata)[1] != "SampleID") {names(metadata)[1] <- "SampleID"}

	covariates <- names(metadata)[-1]

	out <- merge(alpha, metadata, by = "SampleID")

	if (messages==T){
		print("Your covariates are:")
		print(covariates)
		print("Your diversity metrics are:")
		print(metrics)
	}
	return(out)
}

#*******************

#Elise's functions, modified

#Taxa summary plots
#Changes made: simplified taxa level assignment code
#Now takes data frames instead of file names
#Eliminated redundant copies, more logical variable names
taxa_lev_sum <- function (otu_table, metadata_table, minpct=0, minind=1, cov, mincovpct=0, output=F, taxa=F, TL){
	#Sets up assigment of taxa names to numeric value
	taxaLevels <- c(1:7)
	names(taxaLevels) <- c("kingdom", "phylum", "class", "order", "family", "genus","species")
	if(is.element(TL, names(taxaLevels))){
		level <- taxaLevels[[TL]]
	} else{
		message("Invalid taxa level. Options are kingdom, phylum, class, order, family, genus, species.")
		return(NULL)}

	#Setting up otu table to work with
	if(is.data.frame(otu_table)){
		otu <- otu_table[,-1]
	} else if (is.character(otu_table)) {
		otu <- read.table(otu_table, sep = "\t", header = TRUE, as.is = TRUE, comment = "", skip = 1, check.names=F)
		otu <- otu[,-1]
	} else {
		message("OTU table is of type not supported. Please use a file path or data frame and try again.")
		return(NULL)
	}

	n.col <- ncol(otu)
	n.samples <- ncol(otu) - 1
	otu <- otu[,c(n.col,1:n.samples)]  #Puts the taxa column first
	otu[,2:n.col] <- lapply(otu[,2:n.col], as.numeric) #I think this step but good in case somehow nonnumeric.

	#Set up taxa names at correct level
	names_split <- array(dim=c(length(otu$taxonomy), level))
	otu_names <- as.character(otu$taxonomy)
	otu_names <- as.matrix(otu_names)
	for (i in 1:nrow(otu_names)){
		names_split[i,] <- head(strsplit(otu_names[i], "; ", fixed=T)[[1]], n=level)
	}
	otu[,1] <- ldply(apply(names_split, 1, function(x) data.frame(x=paste(x[1:level], sep="", collapse="; "))))

	#Sets up some output dataframes
	#OTU.m is counts by taxonomic group (so reduces the number of rows based on that taxonomic level)
	otu.m <- aggregate(otu[,2:ncol(otu)], by=list(otu$taxonomy), FUN=sum)
	names(otu.m)[1] <- "taxonomy"

	otu.p <- as.matrix(otu.m[,-1])
	rownames(otu.p) <- otu.m$taxonomy

    for(i in 1:ncol(otu.p)){
        otu.p[,i] <- 100*otu.p[,i]/sum(otu.p[,i])
        }

	cutoff <- rowSums(otu.p >= minpct) >= minind
    otu.cut <- as.matrix(otu.p[cutoff, ])
    otuNames <- otu.m$taxonomy[cutoff]

	#Set up the metadata to make sure it only includes samples in the OTU table
	if(is.data.frame(meta_data)){
		metadata <- meta_data
	} else if (is.character(meta_data)) {
		metadata <- read.table(meta_data, sep = "\t", header = TRUE, as.is = TRUE)
	} else {
		message("Meta data table of type not supported. Please use a file path or data frame and try again.")
		return(NULL)
	}

	metadata$SampleID <- as.character(metadata$SampleID)
    metadata <- metadata[which(is.element(metadata$SampleID, colnames(otu.cut))),]
	mm <- match(colnames(otu.cut)[2:ncol(otu.cut)], as.character(metadata$SampleID))
    metadata <- metadata[mm,]

	#Error condition for bad input
	if(is.element(cov, names(metadata))==FALSE){
		message("Invalid covariate. Check spelling etc")
		return(NULL)}
	#Get just the covariate of interest
	metadata.cov <- as.data.frame(na.omit(metadata[,cov]), row.names = metadata$SampleID)
    cov_lev <- as.vector(levels(factor(metadata.cov[,1])))

	#Get things set up in desired order, and make transpose to work with
	otu.cov <- otu.cut[,colnames(otu.cut) %in% rownames(metadata.cov)]
    mm <- match(rownames(metadata.cov), colnames(otu.cov))
    otu.cov <- otu.cov[,mm]
    colnames(otu.cov) <- metadata.cov[,1]
	otu.cov.t <- t(otu.cov)
    otus.len <- length(otu.cov.t[1,])

	#Aggregate relative abundances by covariate level; set up new table for formatting output and stuff
	cov_otu <- aggregate(otu.cov.t[,1:otus.len], by=list(rownames(otu.cov.t)), FUN=mean, rownames=T)
	cov_otu.t <- t(cov_otu[-1])
    colnames(cov_otu.t) <- c(cov_lev)
	cutoff <- rowSums(cov_otu.t >= mincovpct) >= 1
    cov_otu.cut <- cov_otu.t[cutoff, ]

	if (taxa != F){
        cov_otu.cut <- cov_otu.cut[grep(taxa, rownames(cov_otu.cut)),]
        cov_otu.cut <- as.matrix(cov_otu.cut)
        print(rownames(cov_otu.cut))
        }

	otuNames <- as.matrix(rownames(cov_otu.cut))
    cov_otu.cut <- as.matrix(cov_otu.cut)

	#This section does a bunch of string parsing to get all the names correct.
    names=vector()
    for (i in 1:nrow(otuNames)) {
        names[i] <- strsplit(otuNames[i,], "; ", fixed=T)
		}
    names_tail <- array(dim=c(length(names),level))
    for (i in 1:nrow(otuNames)) {
        names_tail[i,] <- tail(strsplit(otuNames[i,], "; ", fixed=T)[[1]], n=level)
		}
    names_tail.df <- as.data.frame(names_tail)
    names_tail.df[,length(names_tail.df[1,])+1] <- ldply(apply(names_tail, 1, function(x) data.frame(x=paste(x[1:level], sep="", collapse="; "))))
    names.n <- as.vector(names_tail.df[,length(names_tail.df[1,])])
    names.N <- vector()
    names.end <- vector()
    for (i in 1:length(names.n)){
        names.N[i] <- gsub('; .__',' ', as.character(names.n[i]))
        names.N[i] <- gsub('k__','', as.character(names.N[i]))
        names.end[i] <- paste(tail(strsplit(names.N[i], " ", fixed=T)[[1]], n=2), collapse=" ")
        }
    rownames(cov_otu.cut) <- names.N
    names <- as.vector(length(names_tail.df[1,]))

	RBcolors <- sample(rainbow(length(names.N)))

	if (output != F){
		cov_otu.return <- cbind(Taxon = rownames(cov_otu.cut), cov_otu.cut)
		write.table(cov_otu.return, output, sep='\t', row.names=F)}

	print("RBcolors = color vector with the same length as the number of taxa; taxa_cov = matrix of taxa abundances grouped by the covariate specified in the function; taxa_names = the names of the taxa in the filtered data set as the last two levels; cov_names = the names of factors in the specified covariate")

	return(list(RBcolors = RBcolors,
                taxa_cov = cov_otu.cut,
                taxa_names = names.end,
                cov_names = cov_lev))
}

#Functions used in ANOVA tutorial (no changes made)
boxCoxNormalize <- function(x){
  if(length(x)<=1){
    x;
  } else{
        if(sum(x<=0)>0){
            x=x-min(x)+0.0001;
        }
        p=powerTransform(x)[6]$start;
        xnorm=(x^p -1)/p
        xnorm
    }
}

#OTU composition summary by taxa level; revisions are similar to taxa_lev_sum
otu_tax_lev <- function(otu_table, TL){
	#Sets up assigment of taxa names to numeric value
	taxaLevels <- c(1:7)
	names(taxaLevels) <- c("kingdom", "phylum", "class", "order", "family", "genus","species")
	if(is.element(TL, names(taxaLevels))){
		level <- taxaLevels[[TL]]
	} else{
		message("Invalid taxa level. Options are kingdom, phylum, class, order, family, genus, species.")
		return(NULL)}

		if(is.data.frame(otu_table)){
			otu <- otu_table[,-1]
		} else if (is.character(alpha_data)) {
			otu <- read.table(otu_table, sep = "\t", header = TRUE, as.is = TRUE, comment = "", skip = 1, check.names=F)
			otu <- otu[,-1]
		} else {
			message("OTU table is of type not supported. Please use a file path or data frame and try again.")
			return(NULL)
		}

    TAX <- ncol(otu)
    LC <- ncol(otu) - 1
    otu <- otu[,c(TAX,1:LC)] #Sticks the taxonomy vector at the front
    otu[,2:TAX] <- lapply(otu[,2:TAX], as.numeric)

	#Sets up names to be at correct taxa level
    names_split <- array(dim=c(length(otu$taxonomy), level))
    otu_names <- as.character(otu$taxonomy)
    otu_names <- as.matrix(otu_names)
    for (i in 1:nrow(otu_names)){
        names_split[i,] <- head(strsplit(otu_names[i], "; ", fixed=T)[[1]], n=level)
    }
    otu[,1] <- ldply(apply(names_split, 1, function(x) data.frame(x=paste(x[1:level], sep="", collapse="; "))))

	#Does the aggregation and sets up matrix to be returned.
	otu.merged <- aggregate(otu[,2:ncol(otu)], by=list(otu$taxonomy), FUN=sum)
    rownames(otu.merged) <- otu.merged[,1]
    otuNames <- otu.merged[,1]
    otu.merged <- as.matrix(otu.merged[,-1])
    return(list(otu_table = otu.merged, otuNames = as.character(otuNames)))
}

#Organizes alpha diversity values based on covariates. This function doesn't calculate anything new, it just makes the format from R a bit more digetstible. I changed the input so now it can take either data frames or files; it checks if the covariates input actually exist
ALPHA_cov = function (alpha_data, meta_data, Metric="PD_whole_tree", cov1, cov2=F){

	#Reads in alpha diversity data from qiime output
	if(is.data.frame(alpha_data)){
		alpha <- alpha_data
	} else if (is.character(alpha_data)) {
		alpha <- read.table(alpha_data, sep = "\t", header = TRUE, as.is = TRUE)
	} else {
		message("Alpha diversity of type not supported. Please use a file path or data frame and try again.")
		return(NULL)
	}

	if (!(is.element(Metric, names(alpha)))){
		message("Metric not found - please check spelling and your metric list.")
		return(NULL)
	}

	#I'm not sure why colors are hard-coded; I might change this to use color palettes instead later
	colorL <- c('purple2','gold','limegreen','red','gray','royalblue2','darkorange1','greenyellow','turquoise2', 'magenta4', 'yellow', 'blue', 'slateblue1', 'springgreen', 'orchid3', 'tan4', 'palevioletred1', 'maroon')

	#Changing row names to be the first column (ie, sample ids)
	IDs <- alpha[,1]
    alpha <- alpha[,-1]
    rownames(alpha) <- IDs

	#Read in metadata
	if(is.data.frame(meta_data)){
		metadata <- meta_data
	} else if (is.character(meta_data)) {
		metadata <- read.table(meta_data, sep = "\t", header = TRUE, as.is = TRUE)
	} else {
		message("Meta data table of type not supported. Please use a file path or data frame and try again.")
		return(NULL)
	}

    #Subsets data to match correct samples and puts in order
	metadata <- metadata[as.character(metadata$SampleID) %in% rownames(alpha),]
    mm <- match(rownames(alpha), as.character(metadata$SampleID))
    metadata <- metadata[mm,]

	if (is.element(cov1, names(metadata))){
		Cov <- metadata[,cov1]
	} else{
		message("cov1 is not found. Please enter a valid covariate.")
		return(NULL)
	}

	#Sets up covariate 1 metadata
	covs <- as.vector(levels(factor(Cov)))
    metadata.cov <- as.data.frame(Cov)
    rownames(metadata.cov) <- metadata$SampleID
    metadata.cov <- as.data.frame(na.omit(metadata.cov))

	#Sets up alpha div matrix to match up with covariate 1
    alpha <- alpha[rownames(alpha) %in% rownames(metadata.cov),]
    mm <- match(rownames(metadata.cov), rownames(alpha))
    alpha <- alpha[mm,]
    alpha <- as.matrix(alpha)
    rownames(alpha) = Cov

	#Sets up colors
	colorV <- c()
    for(i in seq_len(NROW(metadata.cov[,1]))){
        for (j in 1:length(covs)){
            if ((metadata.cov[i,]) == covs[j]) {colorV[i]=paste(colorL[j])}
    }}

	#Does everything above for cov2 if included.
	if (cov2 != F){
		if (is.element(cov2, names(metadata))){
			Cov2 <- metadata[, cov2]
		} else {
			message("cov2 is not found. Please enter a valid covariate.")
			return(NULL)
		}

		covs2 <- as.vector(levels(factor(Cov2)))
		metadata.cov2 <- as.data.frame(Cov2)
		rownames(metadata.cov2) <- metadata$SampleID
		metadata.cov2 <- as.data.frame(na.omit(metadata.cov2))

		colorV2 <- list()
		for(i in seq_len(nrow(metadata.cov2[,1]))){
			for (j in 1:length(covs2)){
				if ((metadata.cov2[i,]) == covs2[j]) {colorV2[i]=paste(colorL[j])}
			}
		}
	} else {
		Cov2 = "NA"
		colorV2 <- "not included"
	}


	#Sets up output
    alphaT <- alpha #All the metrics, rows named for cov 1 values.
	alpha_DI <- alphaT[,grep(Metric, colnames(alphaT))] #Limit to metric of choice

	#I'm really unclear what this section does
	alpha_list <- list()
    for(i in 1:length(unique(Cov))){
        alpha_list[[i]] = length(unique(Cov))
    }
    for(i in seq(covs)){
        alpha_list[[i]]=alpha[grep(covs[i], rownames(alpha)),]
        }
    alpha_list_single <- list()
    for(i in 1:length(as.vector((levels(factor(Cov)))))){
        alpha_list_single[[i]]=alpha_list[[i]][,grep(Metric, colnames(alpha_list[[1]]))]
        }
    max.len <- max(sapply(alpha_list_single, length))
    alpha.c <- lapply(alpha_list_single, function(x) {c(x, rep(NA, max.len - length(x)))})

	#This output is all the diversity values
	alphaM <- matrix(unlist(alpha.c, use.names=T), ncol=length(alpha_list), byrow=F, dimnames=list(NULL,as.vector(levels(factor(Cov)))))
    rows <- max.len

	print("alphaM: alpha diversity matrix for chosen covariate 1 and metric")
    print("alphaT: alpha diversity matrix including all diversity metrics in input table as columns, where rownames are the sample specific factor for variable 1")
    print("alphaD: vector of alpha diversity values for the chosen diversity metric for each sample where names are factor IDs for variable 1")
    print("Cov1: vector of sample-associated factors for variable 1")
	print("Cov2: vector of sample-associated factors for variable 1")
	print("var.names: vector of factors for variable 1 in the same order as columns of alphaM")
	print("rows: the number of rows for jitter")
	print("colorV1: color vector where each color codes for the sample-associated factor in variable 1")
	print("colorV2: color vector where each color codes for the sample-associated factor in variable 2")

	return(list(alphaM = alphaM,
                alphaT = alphaT,
                alphaD = alpha_DI,
                Cov1 = Cov,
				Cov2 = Cov2,
				var.names = covs,
				rows = rows,
				colorV1 = colorV,
				colorV2 = colorV2))
}

#PCOA function
#This one's fairly straightforward - made the same adjustments I've made on others. Can take in either file or data frame, and checks if covariate is valid. This one I changed the color assignment a little, so the actual values aren't hardcoded. Could be made even more flexible.
#Also pulling out library calls assuming; this will be run as a part of a package
#Requires vegan and ape.
PCOA <- function (beta_data, meta_data, cov){

	if (is.data.frame(beta_data)){
		beta <- beta_data
	} else if (is.character(beta_data)) {
		beta <- read.table(beta_data, sep='\t', , header=T, as.is=TRUE, check.names=F)
	} else {
		message("Beta diversity data is wrong format - please use a dataframe or a file name/path.")
		return(NULL)}

	if (is.data.frame(meta_data)){
		metadata <- meta_data
	} else if (is.character(beta_data)) {
		metadata <- read.table(meta_data, sep='\t', , header=T, as.is=TRUE, check.names=F, comment="")
	} else {
		message("Metadata is wrong format - please use a dataframe or a file name/path.")
		return(NULL)}

	#The rest of the function assumes metadata$SampleID exists. Return error if it doesn't.
	if (!(is.element("SampleID", names(metadata)))){
		message("SampleID needs to be a metadata column. Check your metadata table/header names and try again.")
		return(NULL)
	}

	IDs <- beta[,1]
    beta.m <- as.matrix(beta[-1], dimnames=list(IDs, IDs))

	#Actually does the PCOA calculation.
    PCOA <- pcoa(beta.m)$vectors
    rownames(PCOA) <- IDs

	#Puts rows of metadata file in same order as PCOA matrix
    metadata = metadata[as.character(metadata$SampleID) %in% rownames(PCOA),]
    mm = match(rownames(PCOA), as.character(metadata$SampleID))
    metadata = metadata[mm,]

	#Pulls out just covariate column of interest
	Cov <- metadata[,cov]
    covs <- as.vector(levels(factor(Cov)))
    metadata.cov <- as.data.frame(Cov)
    rownames(metadata.cov) <- metadata$SampleID

	#Sets up colors
    colorL <- rainbow(length(covs))
    colorV <- c()
    for(i in seq_len(NROW(metadata.cov[,1]))){
        for (j in 1:length(covs)){
            if ((metadata.cov[i,]) == covs[j]) {colorV[i]=paste(colorL[j])}
    }}
    print("PC = matrix of axes of variation; Cov = vector for sample identities for chosen variable; var.names = names of factors for chosen variable; and colorV = color vector for chosen variable")
    return(list(PC = PCOA,
                Cov = Cov,
                var.names = unique(Cov),
                colorV = colorV))
}

#BETA_cov function for testing statistical significance
#Pulling package reshape to the beginning of file.
BETA_cov <- function (beta_data, meta_data, cov, output=F){

	#Read in data files (either data frame or filename)
	if (is.data.frame(beta_data)){
		beta <- beta_data
	} else if (is.character(beta_data)) {
		beta <- read.table(beta_data, sep='\t', , header=T, as.is=TRUE, check.names=F)
	} else {
		message("Beta diversity data is wrong format - please use a dataframe or a file name/path.")
		return(NULL)}

	if (is.data.frame(meta_data)){
		metadata <- meta_data
	} else if (is.character(beta_data)) {
		metadata <- read.table(meta_data, sep='\t', , header=T, as.is=TRUE, check.names=F, comment="")
	} else {
		message("Metadata is wrong format - please use a dataframe or a file name/path.")
		return(NULL)}

    IDs <- beta[,1]
    beta.m <- as.matrix(beta[,-1], dimnames = list(IDs, IDs))

	#Data cleanup to match up with metadata file
    beta.m[lower.tri(beta.m, diag=F)] <- NA
    metadata <- metadata[as.character(metadata$SampleID) %in% rownames(beta.m),]
    mm <- match(rownames(beta.m), as.character(metadata$SampleID))
    metadata <- metadata[mm,]

	#Pull out metadata of interest
	Cov <- metadata[,cov]
    covs <- as.vector(unique(Cov))
    metadata.cov <- as.data.frame(Cov)
    rownames(metadata.cov) <- metadata$SampleID
    metadata.cov <- as.data.frame(na.omit(metadata.cov))

	#Set up beta.m to match up with covariate data
    beta.m <- beta.m[rownames(beta.m) %in% rownames(metadata.cov),]
    mm <- match(rownames(metadata.cov), rownames(beta.m))
    beta.m <- beta.m[mm,]
    rownames(beta.m) <- Cov
    colnames(beta.m) <- Cov

	#Bunch of data manipulation that I think something like ddply would do better. In general, somewhat confusing, weird data formatting practices.
	beta.T <- beta.m
    lists_cov=list()
    for(i in 1:length(as.vector(levels(factor(Cov))))){
        lists_cov[[i]] = length(as.vector(levels(factor(Cov))))
    }
    for(i in seq(covs)){
        lists_cov[[i]]=beta.T[grep(covs[i], rownames(beta.T)),]
    }
    lists_within_cov = list()
    for(i in 1:length(as.vector(levels(factor(Cov))))){
        lists_within_cov[[i]] = lists_cov[[i]][,grep(covs[i], colnames(lists_cov[[i]]))]
    }
    lists_intra_means = list()
    for(i in 1:length(as.vector(levels(factor(Cov))))){
        lists_intra_means[[i]] = mean(lists_within_cov[[i]], na.rm=TRUE)
    }
    means = as.array(lists_intra_means)
    rownames(means) <- covs

	#Color vector
	colorL <- rainbow(len(covs))
    colorV <- c()
    for(i in seq_len(NROW(metadata.cov[,1]))){
        for (j in 1:length(covs)){
            if ((metadata.cov[i,]) == covs[j]) {colorV[i]=paste(colorL[j])}
    }}

	#Continuing very confusing data organization
    lists_intra_V1 = list()
    lists_intra_V2 = list()
    lengths = vector()
    for(i in 1:length(as.vector(levels(factor(Cov))))){
        lists_intra_V1[[i]] = as.vector(lists_within_cov[[i]])
        rownames(lists_intra_V1[[i]]) = NULL
        lists_intra_V2[[i]] = as.vector(lists_within_cov[[i]])[!is.na(lists_within_cov[[i]])]
        lengths[i]=length(lists_intra_V1[[i]])
    }
    Max = max(lengths)
    covsL = length(covs)
    Cov_BD = array(dim=c(Max,covsL), dimnames=list(c(), c(covs)))
    colnames(Cov_BD)
    for(i in 1:length(as.vector(levels(factor(Cov))))){
        Cov_BD[1:lengths[i],i] = lists_intra_V1[[i]]
    }
    Cov_BD <- as.data.frame(Cov_BD)
    beta_WT <- matrix(ncol=covsL, nrow=covsL)
    colnames(beta_WT) <- covs
    rownames(beta_WT) <- covs

	#For all the messy formatting etc, this is ultimately the core of what this function does.
    for (i in 1:ncol(Cov_BD)) {
        for (j in 1:ncol(Cov_BD)){
            if (j < i){
				beta_WT[[j,i]] = wilcox.test(Cov_BD[,i], Cov_BD[,j], na.rm=T)$p.value
				} else {beta_WT[j,i] <- NA}}}

	#Returning output.
	if (output != F){
		write.table(beta_WT, output, sep='\t', row.names=FALSE)}
    print("beta.C = beta diversity data by covariate; beta_sig = matrix of p-values; Cov = variable being analyzed; var.names = factors for variable 1; colorV1 = color vector for variable 1")
    return(list(beta.C = lists_within_cov,
				beta.sig = beta_WT,
				Cov = Cov,
				var.names = covs,
				colorV1 = colorV))
}
