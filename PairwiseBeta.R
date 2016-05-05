#Code for plotting and analyzing between-group UniFrac distance. This code uses the tutorial HMP data - replace column and file names with your own covariates of interest. 
library(reshape2)
library(ggplot2)

#Read in your beta diversity matrix - 
beta <- read.table("unweighted_unifrac_dm.txt", header=F, skip = 1, sep="\t",comment="")
SampleIDs <- as.character(beta[,1])
beta <- as.matrix(beta[,-1])
beta[lower.tri(beta, diag=F)] <- NA

metadata <- read.table("HMP_5BS_mapping.txt", header=T, sep="\t",comment="")
if (names(metadata)[1] != "SampleID") {names(metadata)[1] <- "SampleID"}

#Play around with your metadata file however you want (delete unnecessary columns, etc) - make sure you know the names of the columns you want! 

#Now we take our matrix and turn it into a data frame formatted for the rest of our analysis
beta.M <- data.frame(beta)
names(beta.M) <- sampleIDs
beta.M$sub1 <- sampleIDs

Beta <- melt(beta.M, id="sub1", variable.name = "sub2")
Beta <- na.omit(Beta)

Beta$Pairs <- rep(NA, nrow(Beta)) #Makes a blank column for the categories you're interested in. 

#This loop goes through and grabs the covariate values for each of the samples being compared. 
for (i in 1:nrow(Beta)) {
  traits <- c(as.character(metadata$Sex[which(metadata$SampleID == Beta$sub1[i])]), as.character(metadata$Sex[which(metadata$SampleID == Beta$sub2[i])]))
  traits <- sort(traits) #This way, you won't have the problem of "a-b" and "b-a" being named different things
  Beta$Pairs[i] <- paste(traits[1], traits[2], sep = "-")
  
  #If you don't actually care about the different values, but only if they are the same or different, here's a way you could do that. Un-comment to run. 
  #Beta$StatComp[i] <- ifelse(traits[1] == traits[2], "Same","Different")
}

#What if there's some other category you might want to filter by - what if it's only relavent to compare individuals at the same body site? We can make a flag for if there's a match or not. 

Beta$Match <- rep(NA, nrow(Beta))

for (i in 1:nrow(Beta)) {
  traits <- c(as.character(metadata$BodySite[which(metadata$SampleID == Beta$sub1[i])]), as.character(metadata$BodySite[which(metadata$SampleID == Beta$sub2[i])]))
  Beta$Match[i] <- ifelse(traits[1] == traits[2], TRUE,FALSE)
}

#Now you can do stuff with your data: 
summary(aov(value~Pairs, data = Beta))
#same thing, with a filter
summary(aov(value~Pairs, data = subset(Beta, Beta$Match == TRUE)))

#If you want to have a table with the means per category, you can use aggregate: 
means <- aggregate(Beta$value, by = list(Beta$Pairs, Beta$Match), FUN=mean)
names(means) <- c("Pairs", "SiteMatch?","Value")
means

ggplot(Beta, aes(x = Pairs, y = value)) + geom_boxplot()
ggplot(subset(Beta, Beta$Match == TRUE), aes(x = Pairs, y = value)) + geom_boxplot()

#We can even do a pairwise test for significance. This uses the package "lsmeans"
library(lsmeans)
m <- lm(value ~ Pairs, data = Beta)
meancomp <- lsmeans(m, "Pairs")
contrast(meancomp, method="pairwise")

#Play around with subsets, plotting, and tests to see what will work the best for you! 