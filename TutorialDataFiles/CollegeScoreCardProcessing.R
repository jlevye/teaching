college_full <- read.csv("Most-Recent-Cohorts-All-Data-Elements.csv")

metadata <- read.csv("CollegeScorecardDataDictionary-09-08-2015.csv")

#There are 1728 variables in this data set. This is the list of those I think might be fun to play with. I've pulled out a few dozen that are fun to play with. 

keep <- c(4,6,12,15,16,25,28,32:35,289,300,331,333,334, 336:340,382, 453,566,571,577,684, 1474)
college <- college_full[,keep]

#Only keep current schools, main branches (then get rid of those columns)
college <- subset(college, college$main == 1 & college$CURROPER == 1)
college <- college[,!(names(college) %in% c("main","CURROPER"))]

#Let's make some nicer names for some of these:

newNames <- c("Name","State","HighestDegree","Ownership","HBCU","Tribal","Mens","Womens","Religious","AdmitRate","Size","AvgCostYear","InstateTuition","OutStateTuition","RevenuePerStudent","SpendingPerStudent","FacultySalary","FullTimeFaculty","PellGrant","FedLoan","Debt","AnyLoan","PctFemale","MedFamIncome","MedEarnings10Years","CompletionRate8Years")
names(college) <- newNames

#Some of the factor levels here aren't very descriptive. Let's fix them. Sometimes they're in weird orders, so there's some manipulation to do. 

college$HighestDegree <- as.factor(college$HighestDegree)
levels(college$HighestDegree) = c("None","Certificate","Associate","Bachelors","Graduate")
college$Ownership <- as.factor(college$Ownership)
levels(college$Ownership) = c("Public","Private","ForProfit")

levels(college$Religious) <- c("Evangelical","Mainline","Evangelical","Mainline","Evangelical","Mainline","Pentecostal","Evangelical","Catholic","Evangelical","Evangelical","Evangelical","Evangelical","Evangelical","Evangelical","Evangelical","Evangelical","Multi","Anabaptist","Mainline","Evangelical","Pentecostal","Mainline","Mainline","Mainline","Mainline","Mainline","Evangelical","Mainline","Pentecostal","Mainline","Evangelical","Mainline","Mainline","Evangelical","Mainline","Mainline","Mainline","Evangelical","Anabaptist","Mainline","Mainline","Mainline","Evangelical","Mainline","Mainline","Mainline","Mainline","Jewish","Mainline","Evangelical","Evangelical","None","Evangelical","Orthodox","Orthodox","Unitarian","Mormon","Pentecostal","Mainline","Other","NULL")

levels(college$HBCU) <- c("No","Yes","NULL")
levels(college$Tribal) <- c("No", "Yes","NULL")

college$Gender <- ifelse(college$Mens == "1", "Men",ifelse(college$Womens == "1", "Women", "Coed"))
college$Gender <- as.factor(college$Gender)
college <- college[,!(names(college) %in% c("Mens","Womens"))]

#The rest of the data should be numeric, but are currently factors. 
for (i in 9:25) {
  college[,i] <- as.numeric(college[,i])
}

write.csv("college","CollegeScoreCardReduced.csv")
