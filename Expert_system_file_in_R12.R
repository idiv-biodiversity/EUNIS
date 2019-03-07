###          R code of theCzech Expert system - Bruelheide H. & Tichý L.  2019    ###
###   The result was compared and optimised with the output from JUICE program    ###
###                             Running version 11                                 ###

# The core of the expert system are machine-readable assignment rules that
# decide whether a particular plot is assigned to a particular class (in this case
# a EUNIS habitat). Each assignment rule is a logical membership formula (i.e. the whole sentence, 
# such as the formula for "B16a Atlantic and Baltic coastal dune scrub":

# (<#TC Shrubs GR 25> NOT (<#TC Trees GR 25> OR <#TC Native-light-canopy-trees GR 15>)) AND 
# (<$$C DUNES_BOHN EQ Y> AND (<$$C COAST_EEA EQ ATL> OR <$$C COAST_EEA EQ BAL>))

# Each membership formula consists of several logical membership expressions (e.g. <#TC Shrubs GR 25>), 
# which start and end with angle brackets ("<", ">") and are combined by formal logic.
# Each logical membership expression can be evaluated to "T" (TRUE) or "F" (FALSE). In consequence, their 
# combination by logical operators "AND", "OR" and "NOT" also result in "T" or "F".
# Each logical expression has a left-hand and (in most cases also a) right-hand condition, called membership
# conditions, which can be evaluated with the result of a numeric value. The left-hand and right-hand membership conditions
# are compared using the logical operators "GR" (greater), "GE" (greater or equal) and 
# "EQ" (equal). 

# The following code extracts all membership formulas (Step ###), all membership expressions from membership formulas (Step ###),
# and all membership conditions from membership expressions (Step ###).

# Then these logical comparisons are evaluated in reverse order, i.e. from membership conditions to
# expressions to formulas. In other words, the assignment rules are resolved from inside to
# outside.

# Membership conditions from the left- and right-hand side of a membership expression are joined and 
# evaluated together, resulting in numerical values (stored in plot.group, a numeric array of the dimensions of
# 1,564,846 plots x 745 left- and right-hand side membership conditions).


# ... making use of logical matrices
# ... Different to sequential approach. Allows vector-oriented evaluation, much quicker in R
# ... can be used for any expert system, highly flexible.



#remowing all memory
rm(list=ls()) 

library(data.table)
library(fastmatch)
library(dplyr) 
library(readr)

match <- fmatch

setwd("C:\\Daten\\Cocktail2.0\\CzechSystem") #HB
setwd("G:\\0_lubos\\Prace\\2018\\2018 EXPERT SYSTEM") #HB
setwd("~/EUNIS") #server

# Speed up the calculation by using parallel computing
library(parallel) 
detectCores() #20 
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)

#Checking the time duration of the code aq - partial, aq2 - total
aq <- Sys.time()
aq2 <- Sys.time()

### Functions ###
# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)
# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)
# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

###############################################################################
### Step 1: read in vegetation data, header data and the expert system file ###
###############################################################################
# Table input (JUICE menu File > Export > Table > To Database Files)
#DT1  <- read.fwf("table2.txt", widths=(c(9,11,4,11,8)), 
#                 encoding="UTF-8", colClasses = c(rep("numeric",4),"character"))
#DT1  <- read_fwf("table.txt", fwf_widths(widths = c(9,11,4,13,6)), 
#                 col_types = "ddddc")
# EVA_2018_09_24_merged.exp
DT1  <- read_fwf("table.txt", fwf_widths(widths = c(9,11,4,13,6)), 
                 col_types = "ddddc")
#str(DT1) #32957910 obs. of  5 variables:
#DT1 <- DT1[1:108170,-5]
DT1 <- DT1[,-5]
names(DT1) <- c("PlotID","SpeciesID","Layer", "Cover")
DT1 <- data.table(DT1)

#species1  <- read.fwf("species2.txt", widths=c(8,50,1,99),skip=1, 
#                      encoding="UTF-8", colClasses = c("numeric","character","numeric","character"))
species1  <- read_fwf("species.txt", fwf_widths(widths=c(8,55,1,60,8)), col_types="dcddc",skip=1)
str(species1)
species1 <- species1[,-c(4,5)]
names(species1) <- c("SpeciesID", "Taxon_name", "Layer")
species1$Taxon_name <- trim.trailing(species1$Taxon_name)
index1 <- fmatch(DT1$SpeciesID, species1$SpeciesID)
any(is.na(index1)) # F
DT1$Taxon_name <- species1$Taxon_name[index1]
setkey(DT1, PlotID)
plot.names <- unique(DT1$PlotID)
plot.number <- length(plot.names)
#1564846

# Environmental variables
header.definitions <- read.csv("EVA_2018_09_24_merged.str", 
              encoding="UTF-8", stringsAsFactors = F, sep="\t", header=F)
# header.definitions[1,1]
header.definitions2 <- matrix(NA, nrow=dim(header.definitions)[[1]], ncol=2)
for (i in 1:dim(header.definitions)[[1]]){
  header.definitions2[i,1] <- as.numeric(tstrsplit(trim.leading(header.definitions[i,1]), " ")[[1]])
}
header.definitions2
for (i in 2:dim(header.definitions2)[[1]]){
  header.definitions2[i-1,2] <- header.definitions2[i,1] -header.definitions2[i-1,1]
}
header.definitions2

header.definitions2[26,2] <- 728 - sum(header.definitions2[,2], na.rm=T)
sum(header.definitions2[,2])
# 728
head1  <- read_fwf("head.txt", fwf_widths(widths = header.definitions2[,2]), 
                 col_types = "dddcccddddddddccccccccdddc")
'head1  <- read_fwf("head.txt", fwf_widths(widths = 
        c(7,8,10,22,10,11,6,8,550,4,12,4,3,10,10,8,39)), 
                 col_types = "dddcccddcccccdddc")
Shorthead: 1,7 
PlotID: 8,8
TV2_rel_no: 16,10
Country: 26,22
Nr_Table: 48, 10
Date: 58, 11
Area: 69, 6
Altitude (m): 75, 8
x1: 82, 550
BIOREG_EEA: 633, 4
ECOREG_WWF: 637, 12
COAST_EEA: 649, 4
DUNES_BOHN: 653,3
DEG_LON: 656, 10
DEG_LAT: 666, 10
Location_uncertainty: 676, 8
Dataset: 684, 39

7+8+10+22+10+11+6+8+550+4+11+5+3+11+11+12+39
# 728
names(head1) <- c("Shorthead","PlotID","TV2_rel_no",
                  "Country","Nr_Table","Date","Area","Altitude (m)","x1",
"BIOREG_EEA", "ECOREG_WWF", "COAST_EEA",
"DUNES_BOHN", "DEG_LON", "DEG_LAT", 
"Location_uncertainty", "Dataset")
'
str(head1)
head1 <- as.data.frame(head1)
table(head1[,26])
names(head1)
names(head1) <- c("ReleveNumber","PlotID","TV2_rel_no",
                  "Country","Nr_Table","Date","Area","Altitude (m)",
                  "Aspect","Slope","Cover_tree_layer","Cover_shrub_layer",
                  "Cover_herb_layer", "Cover_moss_layer","Mosses_identified",
                  "Lichens_identified", "Locality", "Alliance",
                  "BIOREG_EEA", "ECOREG_WWF", "COAST_EEA",
                  "DUNES_BOHN", "DEG_LON", "DEG_LAT", 
                  "Location_uncertainty", "Dataset")
str(head1)
summary(head1$'Altitude (m)')
table(head1$COAST_EEA)
summary(head1$DEG_LON)


expert1  <- read.csv("Expert-system-EUNIS-2019-01-21_recoded_HB.txt", 
                     encoding="UTF-8", stringsAsFactors = F, sep="\t")
str(expert1) #31517 obs. of  1 variable:
expert1 <- expert1$SECTION.1..Species.aggregation
start <- which(substr(expert1,1,25)=="SECTION 2: Species groups")
end <- which(substr(expert1,1,14)=="SECTION 2: End")
i <- start+1
expert2 <- expert1[(start+1):(end-1)]
index.group.names <- which(substr(expert2,1,1)!=" ")
number.groups <- length(index.group.names) # 316
index.group.names <- c(index.group.names,length(expert2))

for (i in 1: number.groups){
  g1 <- trim.leading(expert2[(index.group.names[i]+1):(index.group.names[i+1]-1)])
  g2 <- list(x=g1)
  names(g2) <- expert2[index.group.names[i]]
  if (i==1){
    groups <- g2
  } else {
    groups <- append(groups,g2)
  }
}

#Data import
Sys.time()-aq
aq<-Sys.time()

###########################################
### Step 2: read in membership formulas ###
###########################################
start <- which(substr(expert1,1,28)=="SECTION 3: Group definitions")
end <- which(substr(expert1,1,14)=="SECTION 3: End")
expert3 <- expert1[(start+1):(end-1)]
# holds formula names and formulas.
# I refer to the whole line as a formula (=sentence).
# Everything in a bracket (<...>) is an expression, which has to be T or F.
# Every expression can have one to many conditions, which are also T or F. 
# (Actually the conditions are expressions themselves, but inner expressions.)
# Expressions consist of a left- and right-hand condition, which are compared
# by a logical operator (GR, GE, EQ).

length(expert3) #632
# The membership expressions contain CR (carriage return characters), thus
# the condition may be distributed across more than one line
membership.formula.names <- NULL
# Name of the formula
membership.expressions <- NULL
# inner expressions, between < and >
membership.formulas <- NULL
# whole expressions
i <- 0
while (i < length(expert3)){
  i <- i + 1
  if (substr(expert3[i],1,3)!="---"){
    # these conditions are out-commented and not used
    if (!grepl("[^0-9]", substr(expert3[i],1,1))){
      # checks whether line is a formula name
      membership.formula.names <- c(membership.formula.names,expert3[i])
    } else {
      # then the line is a formula or an out-commented formula
      c <- expert3[i]
      while (grepl("[^0-9]", substr(expert3[i+1],1,1)) & substr(expert3[i+1],1,1)!="-" &
             i < length(expert3)){
        # checks whether the next line is still the current formula
        # if not, two lines are pasted together
        # !grepl("[^0-9]", x)  checks whether the first character of the next line is numeric
        # substr(expert3[i+1],1,1)!="-" checks whether the first character of the next line is an out-commented line
        i <- i +1
        c <- paste(c, expert3[i],sep=" ")
      }
      a <- gregexpr("<",c, fixed=T)[[1]]
      b <- gregexpr(">",c, fixed=T)[[1]]
      if (a[1]>0) {
        membership.formulas <- c(membership.formulas,c)
        membership.expressions2 <- array("",length(a))
        for (j in 1: length(a)){
          membership.expressions2[j] <- substr(c,a[j]+1,b[j]-1)
        }
        membership.expressions <- c(membership.expressions,membership.expressions2)
      }
    }
  }
}
write.csv(membership.formula.names, "membership.formula.names.csv")
write.csv(membership.formulas, "membership.formulas.csv")
write.csv(membership.expressions, "membership.expressions.csv")


# Expert system decoding 1
Sys.time()-aq
aq<-Sys.time()

########################################################################################
### Step 3: Insert right-hand side variables where are no right-hand side conditions ###
########################################################################################
# Before the expressions can be split at the logical operators GR, GE or EQ
# into a left-hand side and right-hand side condition, some right hand-sides
# have to be complemented, because they are not interpretable 
# without the left-hand side. 
# For example, #T$ occurs without group name on the right hand side and means 
# total cover except the species on the left-hand side. 
# In this case, the group name after #T$ on the right-hand side is inserted.

# the only logical operator occurring in #T$ groups in the EUNIS system is GR
# Note, that it has to be checked whether this is also the case in other expert systems

# check whether "#T$" is actually on the right hand side
index3 <- which(regexpr("#T$",membership.expressions, fixed=T)>0 & 
                trim(tstrsplit(membership.expressions,"GR", fixed=T)[[2]])=="#T$")

# not all of them are unique
b <- unique(membership.expressions[index3])
a <- tstrsplit(b,"GR", fixed=T)
a[[1]] <- trim(a[[1]])


for (i in 1:length(b)){
  index4 <- which(regexpr(b[i],
                          membership.formulas, fixed=T)>0)
  # to handle the left-hand side expression <#TC Dry-and-wet-heath-shrubs|#TC Dry-heath-shrubs GR #T$>
  # "#T$" has also to be inserted inside the string
  a[[1]][i] <- gsub("#TC","#T$",a[[1]][i], fixed=T)
  membership.formulas[index4] <- gsub(b[i],
                                      paste(b[i],
                                            substr(a[[1]][i],4,nchar(a[[1]][i])), sep=""),
                                      membership.formulas[index4],
                                      fixed=T)
}
# now also change this in the membership expressions
a <- tstrsplit(membership.expressions[index3],"GR", fixed=T)
a[[1]] <- trim(a[[1]])
# to handle the left-hand side expression <#TC Dry-and-wet-heath-shrubs|#TC Dry-heath-shrubs GR #T$>
# "#T$" has also to be inserted inside the string
a[[1]] <- gsub("#TC","#T$",a[[1]], fixed=T)

membership.expressions[index3] <- paste(membership.expressions[index3],
                                        substr(a[[1]],4,nchar(a[[1]])),sep="")

# we identify all conditions without a right-hand side
# <##D Diagnostic species group>	
# The number of species of the diagnostic species group is greater than 
# the number of species of any other diagnostic species group defined in section 2. 
# <##C Diagnostic species   group>	
# The total cover of the diagnostic species group is greater than the
# total cover of any other diagnostic species group defined in section 2.
# <##Q Diagnostic species group>	
# The sum of the square root cover of the diagnostic species group is greater 
# than the square root cover of any other diagnostic species group defined in section 2.

# these conditions are replaced with "GR NON", meaning that either number, cover or 
# sum of squared cover is greater than on the left-hand side
is.not.right.hand.side <- regexpr("GR",membership.expressions, fixed=T)==-1 &
  regexpr("GE",membership.expressions, fixed=T)==-1 &
  regexpr("EQ",membership.expressions, fixed=T)==-1

# We take only those that do not have numeric condition (i.e. #01, #02 etc.)
index9 <- as.numeric(substr(membership.expressions[is.not.right.hand.side],3,3))
# the warning "NAs introduced by coercion" can be ignored!

a <- membership.expressions[is.not.right.hand.side][is.na(index9)]
any(duplicated(a)) #T
# duplicated membership expressions result in inserting "GR NON" in the loop below
# more than once, thus duplicates are removed
a <- unique(a)

for (i in 1:length(a)){
  index4 <- which(regexpr(a[i], membership.formulas, fixed=T)>0)
  membership.formulas[index4] <- gsub(a[i],
            paste(a[i], "GR NON",a[i], sep=" "),
            membership.formulas[index4], fixed=T)
}

membership.expressions[is.not.right.hand.side][is.na(index9)] <- paste(membership.expressions[is.not.right.hand.side][is.na(index9)],
                                                                       "GR NON", membership.expressions[is.not.right.hand.side][is.na(index9)],sep=" ")
index9 <- regexpr("GR",membership.expressions, fixed=T)==-1 &
  regexpr("GE",membership.expressions, fixed=T)==-1 &
  regexpr("EQ",membership.expressions, fixed=T)==-1
membership.expressions[index9]
# there are 142 expressions without a right-hand side condition
# but these are all numeric expressions


# Expert system decoding 2
Sys.time()-aq
aq<-Sys.time()

#################################################################################################
### Step 4: evaluation of all membership conditions, irrespective of left- or right-hand side ###
#################################################################################################
membership.conditions2 <- unlist(strsplit(membership.expressions,"GR"))
membership.conditions2 <- unlist(strsplit(membership.conditions2,"GE"))
membership.conditions2 <- unlist(strsplit(membership.conditions2,"EQ"))
membership.conditions2 <- trim(membership.conditions2)
membership.conditions2 <- sort(unique(membership.conditions2))
#remove numerical elements
membership.conditions2 <- membership.conditions2[-which(!is.na(as.numeric(membership.conditions2)))]
'Warning message:
In which(!is.na(as.numeric(membership.conditions2))) :
NAs introduced by coercion
-> Warning can be ignored'
#membership.conditions2 #688
length(membership.conditions2)
# 745
write.csv(membership.conditions2, "membership.conditions2.csv")

plot.group <- array(0, c(plot.number, length(membership.conditions2)),
                    dimnames=list(plot.names, membership.conditions2))
# array to collect all membership conditions of the left- and right-hand side

# Expert system decoding 3
Sys.time()-aq
aq<-Sys.time()

### go through all different types of conditions ###
# "###": number of species of a group
index3 <- which(substr(membership.conditions2,1,3)=="###")
# Note that "###" is identical with "##D"
membership.conditions3 <- membership.conditions2
membership.conditions3 <- gsub("###","##D",membership.conditions3 )

for (i in index3){
  index2 <- DT1$Taxon_name %in% groups[[c(which(names(groups)==membership.conditions2[i]),
                                        which(names(groups)==membership.conditions3[i]))]]
  # concatenating the results from membership.conditions2 and membership.conditions3 accounts
  # for looking for "###" as well as "##D" at the beginning of the group name
  result <- DT1[index2, list(x=length(Taxon_name)),by=PlotID]
  index4 <- fmatch(dimnames(plot.group)[[1]],result$PlotID)
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[i])
  plot.group[,index5] <- result$x[index4]
}

'##Q square root of the sum of cover'
index3 <- which(substr(membership.conditions2,1,3)=="##Q")
' ##Q is identical with ##D'
membership.conditions3 <- membership.conditions2
membership.conditions3 <- gsub("##Q","##D",membership.conditions3 )
if (length(index3)>0) {
  ip1<-as.vector(dimnames(plot.group)[[1]])
  ip2<-as.vector(dimnames(plot.group)[[2]])
  nmg<-names(groups)
  pepa<-foreach (i = 1:length(index3),.packages=c("data.table", "fastmatch"), .combine=cbind)%dopar%{
    index2 <- DT1$Taxon_name %chin% groups[[which(nmg==membership.conditions3[index3[i]])]]
    result <- DT1[index2, list(x=sum(Cover^0.5)),by=PlotID]
    index4 <- fmatch(ip1,result$PlotID)
    result$x[index4]
  }
  for(i in 1:length(index3)){
    index5 <- which(ip2==membership.conditions2[index3[i]])
    plot.group[,index5]<-pepa[,i] 
  }
}


'##C total cover of the group'
index3 <- which(substr(membership.conditions2,1,3)=="##C")
i <- index3[1]
' ##C is identical with ##D'
membership.conditions3 <- membership.conditions2
membership.conditions3 <- gsub("##C","##D",membership.conditions3 )
if (length(index3)>0) {
  ip1<-as.vector(dimnames(plot.group)[[1]])
  ip2<-as.vector(dimnames(plot.group)[[2]])
  nmg<-names(groups)
  pepa<-foreach (i = 1:length(index3),.packages=c("data.table","fastmatch"), .combine=cbind)%dopar%{
    index2 <- DT1$Taxon_name %chin% groups[[which(nmg==membership.conditions3[index3[i]])]]
    result <- DT1[index2, list(x=sum(Cover)),by=PlotID]
    index4 <- fmatch(ip1,result$PlotID)
    result$x[index4]
  }
  for(i in 1:length(index3)){
    index5 <- which(ip2==membership.conditions2[index3[i]])
    plot.group[,index5]<-pepa[,i] 
  }
  pepa<-NULL
}  

'#01 to #99 minimum number of species which have to be present in a group'
# this is the only column that is evaluated directly without later 
# comparison
index3 <- which(!is.na(as.numeric(substr(membership.conditions2,2,3))) & 
                  substr(membership.conditions2,1,1)=="#")
'Warning message:
In which(!is.na(as.numeric(substr(membership.conditions2, 2, 3)))) :
NAs introduced by coercion
-> Warning can be ignored'
' #01 to #99 is identical with ###'
membership.conditions3 <- membership.conditions2
substr(membership.conditions3[index3],1,3) <- "###"
for (i in index3){
  if (length(which(names(groups)==membership.conditions3[i]))>0){
    index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==membership.conditions3[i])]]
  } else {
    substr(membership.conditions3[i],3,3) <- "D"
    index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==membership.conditions3[i])]]
  }
  result <- DT1[index2, list(x=length(Taxon_name)),by=PlotID]
  result$y <- ifelse(result$x>=as.numeric(substr(membership.conditions2[i],2,3)),1,0)
  index4 <- fmatch(dimnames(plot.group)[[1]],result$PlotID)
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[i])
  plot.group[,index5] <- result$y[index4]
}

# Set of conditions 1
Sys.time()-aq
aq<-Sys.time()

' #T$ total cover of all other species except those on the left-hand side'
total.cover <- function(cover.values){
  a <- length(cover.values)
  b <- cover.values[1]
  if (a > 1){
    for (i in 2: length(cover.values)){
      b <- b+(100-b)*cover.values[i]/100
    }
  }
  return(b)
}

# check whether there is "#T$" as single condition also on the left-hand side
# this occurs only with "GE" as operator:  "#T$ GE 30"
# and means that total cover is greater or equal than 30
index3a <- which(regexpr("#T$",membership.expressions, fixed=T)>0 & 
                  trim(tstrsplit(membership.expressions,"GE", fixed=T)[[1]])=="#T$")
membership.expressions[index3a]
# this occurs 24 times.

# identify all conditions where "#T$" is mentioned without group name.
# In these cases "#T$" stands on the right hand side.
index3a <- which(membership.conditions2=="#T$")
index3a # 480, i.e. this only occurs once
membership.conditions2[index3a]
result <- DT1[, list(x=total.cover(Cover)),by=PlotID]
index4 <- fmatch(dimnames(plot.group)[[1]],result$PlotID)
index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[index3a])
plot.group[,index5] <- result$x[index4]


# we exclude "|" and "EXCEPT" because we handle these cases
# separately below
index3 <- which(substr(membership.conditions2,1,3)=="#T$" &
                  regexpr("|",membership.conditions2, fixed=T)==-1 &
                  regexpr("EXCEPT",membership.conditions2, fixed=T)==-1)
# This does only work because we have included the left-hand side 
# to the right hand side after #T$ 
' #T$ is identical with #TC and #TC is calculated from ###'
# remove the case where "#T$" occurs alone on the left-hand side
index3 <- index3[-which(index3==index3a)]
membership.conditions3 <- membership.conditions2
substr(membership.conditions3[index3],1,3) <- "###"
for (i in index3){
  if (length(which(names(groups)==membership.conditions3[i]))>0){
    index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==membership.conditions3[i])]]
  } else {
    substr(membership.conditions3[i],3,3) <- "D"
    index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==membership.conditions3[i])]]
  }
  #now exclude all species of this group from total cover calculation
  result <- DT1[!index2, list(x=total.cover(Cover)),by=PlotID]
  index4 <- fmatch(dimnames(plot.group)[[1]],result$PlotID)
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[i])
  plot.group[,index5] <- result$x[index4]
}

# Set of conditions 2
Sys.time()-aq
aq<-Sys.time()

' NON conditions'
# we included them to have a right-hand side in all logical expressions
index3 <- which(substr(membership.conditions2,1,3)=="NON")
# This does only work because we have included the the left-hand side 
# to the right hand side after #T$ 
#membership.conditions2[index3]
# 145 conditions
# 1, 2 and 3 is species number (N), also coded as "D"
# C is cover (C), this does not occur any more (previously in "### Taiga-herbs")
# I have removed this statement as it was out-commented
# Q is sum of squared cover (Q)
index6 <- which(substr(names(groups),3,3)=="D")
# these are all differential groups used for a comparison

# find the corresponding variables 
membership.conditions3 <- membership.conditions2
membership.conditions3[index3] <- substr(membership.conditions2[index3],5,
                                         nchar(membership.conditions2[index3]))

' most #01 etc. groups are is calculated from ##D'
substr(membership.conditions3[index3],1,3) <- "##D"

a <- match(membership.conditions3[index3],names(groups)[index6])
any(is.na(a)) # F, all membership.conditions3 are found in the group names
b <- match(names(groups)[index6],membership.conditions3[index3])
any(is.na(b)) # T, there is one NA
which(is.na(b))
# 84
membership.conditions3[index3][84]
# "##D +03 I1-Arable-land-and-market-gardens"

# we need to calculate species number (N), cover (C) and sum of squared
# cover (Q) of all these groups 
plot.group.non <- array(0, c(plot.number, length(index6), 3),
                        dimnames=list(plot.names, names(groups)[index6], c("N","C","Q")))
for (i in 1: length(index6)){
  index2 <- DT1$Taxon_name %fin% unlist(groups[index6][i], use.names=FALSE)
  result <- DT1[index2, list(x=length(Taxon_name)),by=PlotID]
  index4 <- fmatch(dimnames(plot.group.non)[[1]],result$PlotID)
  index5 <- which(dimnames(plot.group.non)[[2]]==names(groups[index6][i]))
  plot.group.non[,index5,1] <- result$x[index4]
  result <- DT1[index2, list(x=sum(Cover)),by=PlotID]
  plot.group.non[,index5,2] <- result$x[index4]
  result <- DT1[index2, list(x=sum(Cover^0.5)),by=PlotID]
  plot.group.non[,index5,3] <- result$x[index4]
}

# Set of conditions 3
Sys.time()-aq
aq<-Sys.time()

index7 <- fmatch(membership.conditions3[index3],names(groups))
# although there are no NAs, there might be some ###
any(is.na(index7)) # F, all membership.conditions3 are found in the group names
#substr(membership.conditions3[index3][is.na(index7)],1,3) <- "###"
#index7 <- fmatch(membership.conditions3[index3],names(groups))

# now calculate the NON condition
#membership.conditions2[i]
#"NON ##Q X01-Papaveretea-rhoeadis"
#membership.conditions3[i]
"##D X01-Papaveretea-rhoeadis"
c <- 0
index5 <- vector(length=length(index3)) #145 
pgna <- dimnames(plot.group.non)[[2]]
pga <- dimnames(plot.group)[[2]]
result<-foreach (j=1:length(index3), .combine=cbind)%dopar%{
  a <- substr(membership.conditions2[index3[j]],7,7)
  b <- 0
  if (a=="N" | a=="D"){
    # in membership.conditions2 all groups are coded "D", rather than "N"
    # thus we ask for both options
    b <- 1  
  } else {
    if (a=="C"){
      b <- 2
    } else {
      # then it is "Q"
      b <- 3
    }
  }

  # find the name of the group in dimnames(plot.group.non)[[2]]
  index9 <- which(pgna==membership.conditions3[index3[j]])
  # find maximum species number across any other group
  # b is: 1==N, 2==C, 3==Q

  # now calculate the maximum number of species in any other group
  # except the target group
  lada<-plot.group.non[,-index9,b]
  apply(lada,1,FUN=max, na.rm=T)
}
# There are warnings() for maximum calculations for which there are only NA
# These can be ignored
for (i in 1:length(index5)){
  index5 <- which(pga==membership.conditions2[index3[i]])
  plot.group[,index5] <- result[,i]
}
any(is.na(colSums(plot.group[,index7]))) #F
any(colSums(plot.group[,index7])==0) #F
# all data columns are filled

which(dimnames(plot.group)[[2]]=="##Q +03 R1Q-Inland-sanddrift-and-dune-with-siliceous-grassland")
# 438
sum(plot.group[,438]) # 365877.9
which(dimnames(plot.group)[[2]]=="NON ##Q +03 R1Q-Inland-sanddrift-and-dune-with-siliceous-grassland")
# 270
sum(plot.group[,270]) #37325102
membership.expressions[61]
membership.expressions[436]
# "##Q +03 R1Q-Inland-sanddrift-and-dune-with-siliceous-grassland GR NON ##Q +03 R1Q-Inland-sanddrift-and-dune-with-siliceous-grassland"
# occurs two times

# Set of conditions 4
Sys.time()-aq
aq<-Sys.time()

' #TC total cover of all species in this group'
# we exclude "|" and "EXCEPT" because we handle these cases
# separately below
index3 <- which(substr(membership.conditions2,1,3)=="#TC" &
                  regexpr("|",membership.conditions2, fixed=T)==-1 &
                  regexpr("EXCEPT",membership.conditions2, fixed=T)==-1)
' #TC is identical with ###'
membership.conditions3 <- membership.conditions2
substr(membership.conditions3[index3],1,3) <- "###"

# also includes | operator, which is a combinaton of the two groups
# the species groups have to be combined
# "#TC Atlantic-heath-shrubs|#TC Lowland-to-alpine-heath-shrubs"

for (i in index3){
   if (length(which(names(groups)==membership.conditions3[i]))>0){
    index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==membership.conditions3[i])]]
  } else {
    substr(membership.conditions3[i],3,3) <- "D"
    index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==membership.conditions3[i])]]
  }
  result <- DT1[index2, list(x=total.cover(Cover)),by=PlotID]
  index4 <- fmatch(dimnames(plot.group)[[1]],result$PlotID)
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[i])
  plot.group[,index5] <- result$x[index4]
  #remove all situations, where zero is GE or EQ to zero
  plot.group[which(is.na(plot.group[,index5])),index5]<- -1
}

# Set of conditions 5
Sys.time()-aq
aq<-Sys.time()

### the rest are species conditions
# find all species names by taking out "#", "$" and "NON" commands,
# as well as header fields
unique(head1$COAST_EEA)
# "MED" "N"   NA    "BLA" "ATL" "ARC" "BAL"
unique(head1$Country)
sort(unique(head1$Dataset))
#"Swedish_National_Forest_Inventory"
na.omit(unique(head1$DUNES_BOHN))
# [1] "N" "Y"
unique(head1$ECOREG_WWF)

index3 <- which(substr(membership.conditions2,1,1)!="#" & 
                  substr(membership.conditions2,1,1)!="$" &
                  substr(membership.conditions2,1,3)!="NON" &
                  !membership.conditions2 %in% unique(head1$COAST_EEA)&
                  !membership.conditions2 %in% unique(head1$Country)&
                  !membership.conditions2 %in% unique(head1$Dataset)&
                  !membership.conditions2 %in% na.omit(unique(head1$DUNES_BOHN))&
                  !membership.conditions2 %in% unique(head1$ECOREG_WWF))
membership.conditions2[index3]
# only species left

for (i in index3){
  membership.conditions2[i]
  index2 <- DT1$Taxon_name %fin% membership.conditions2[i]
  result <- DT1[index2, list(x=mean(Cover)),by=PlotID]
  index4 <- fmatch(dimnames(plot.group)[[1]],result$PlotID)
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[i])
  plot.group[,index5] <- result$x[index4]
}

# Set of conditions 6
Sys.time()-aq
aq<-Sys.time()

### now handle remaining cases of | and EXCEPT
index3 <- which(regexpr("|",membership.conditions2, fixed=T)>0
                | regexpr("EXCEPT",membership.conditions2, fixed=T)>0)
# most (66 out of 67) conditions have #TC commands
# however there is one with #T$

# first: deal with "#TC"
index3 <- which((regexpr("|",membership.conditions2, fixed=T)>0
                | regexpr("EXCEPT",membership.conditions2, fixed=T)>0) &
                  substr(membership.conditions2,1,3)!="#T$" )
' #TC is identical with ###'
membership.conditions3 <- membership.conditions2
membership.conditions3 <- gsub("#TC", "###",membership.conditions3)

#plot.group.internal <- array(NA, c(plot.number, length(membership.conditions3)),
#                             dimnames=list(plot.names, membership.conditions3))
for (i in index3){
  if (regexpr("EXCEPT",membership.conditions3[i], fixed=T)<0){
    # then the operator is "|", i.e. all all groups are to be combined
    a <- unlist(strsplit(membership.conditions3[i],"|", fixed=T))
    b <- NULL
    for (j in 1: length(a)){
      if (length(which(names(groups)==a[j]))>0){
        index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==a[j])]]
      } else {
        substr(a[j],3,3) <- "D"
        index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==a[j])]]
      }
      c <- groups[[which(names(groups)==a[j])]]
      # combine both groups
      b <- unique(c(b,c))
    }
  } else {
    # 1. divide condition at EXCEPT in left- and right hand side
    # strsplit at EXCEPT
    d <- unlist(strsplit(membership.conditions3[i],"EXCEPT", fixed=T), use.names=FALSE)
    # 2. combine all groups of left- and right-hand side
    # left-hand side
    a <- unlist(strsplit(d[1],"|", fixed=T),use.names=FALSE)
    a <- trim(a)
    b1 <- NULL
    for (j in 1: length(a)){
      if (length(which(names(groups)==a[j]))>0){
        index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==a[j])]]
      } else {
        substr(a[j],3,3) <- "D"
        index2 <- DT1$Taxon_name %fin% groups[[which(names(groups)==a[j])]]
      }
      c <- groups[[which(names(groups)==a[j])]]
      # combine both groups
      b1 <- unique(c(b1,c))
    }
    # right-hand side
    a <- unlist(strsplit(d[2],"|", fixed=T), use.names=FALSE)
    a <- trim(a)
    b2 <- NULL
    for (j in 1: length(a)){
      if (length(which(names(groups)==a[j]))>0){
        index2 <- DT1$Taxon_name %fin% groups[which(names(groups)==a[j])]
      } else {
        substr(a[j],3,3) <- "D"
        index2 <- DT1$Taxon_name %fin% groups[which(names(groups)==a[j])]
      }
      c <- groups[[which(names(groups)==a[j])]]
      # combine both groups
      b2 <- unique(c(b2,c))
    }
    # 3. substract right- hand group from left hand group.
    # now calculate b1 EXCEPT b2
    index7 <- fmatch(b1,b2)
    # only retain species in b1 that are not in b2
    b <- b1[is.na(index7)]
  }
  index2 <- DT1$Taxon_name %fin% b
  result <- DT1[index2, list(x=total.cover(Cover)),by=PlotID]
  index4 <- fmatch(dimnames(plot.group)[[1]],result$PlotID)
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[i])
  plot.group[,index5] <- result$x[index4]
}

# second: deal with "#T$"
index3 <- which((regexpr("|",membership.conditions2, fixed=T)>0
                 | regexpr("EXCEPT",membership.conditions2, fixed=T)>0) &
                  substr(membership.conditions2,1,3)=="#T$" )
membership.conditions2[index3]
# there is only one
# "#T$ Dry-and-wet-heath-shrubs|#T$ Dry-heath-shrubs"
# which should be evaluated as total cover of these two groups
# is higher than of any other group

membership.conditions3 <- membership.conditions2
membership.conditions3[index3] <- gsub("#T$","###", membership.conditions3[index3], fixed=T)
a <- strsplit(membership.conditions3[index3],"|", fixed=T)
# this is a list of groups that have to be combined
length(a)
# 1, currently only one element, but programmed for more
for (i in length(a)){
  #unique(unlist(groups[which(names(groups) %in% a[[i]])]))
  index2 <- DT1$Taxon_name %fin% unique(unlist(groups[which(names(groups) %in% a[[i]])]))
  #now exclude all species of this group from total cover calculation
  result <- DT1[!index2, list(x=total.cover(Cover)),by=PlotID]
  index4 <- fmatch(dimnames(plot.group)[[1]],result$PlotID)
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[index3[i]])
  plot.group[,index5] <- result$x[index4]
}

### now evaluate header data
# There are some formulas that are based on header data.
# There are two group types: "$$C" (character) and "$$N" (numeric). 
# "$$C" is combined with the operator "EQ", 
# "$$N" is combined with all three operators "GR", "GE", "EQ". 

# The name of the group (e.g. "$$C COAST_EEA") is the name of column in the header data table. 
# there are five "##C" fields used in the code: 
# $$C COAST_EEA
# $$C Country
# $$C Dataset
# $$C DUNES_BOHN
# $$C ECOREG_WWF
# there are three "##N" fields used in the code: 
# $$N Altitude (m)
# $$N DEG_LAT
# $$N DEG_LON
# Numeric is used mainly as <$$N DEG_LAT> GR 40 or <$$N DEG_LON GR -10>.

unique(head1$COAST_EEA)
# "MED" "N"   NA    "BLA" "ATL" "ARC" "BAL"
unique(head1$Country)
sort(unique(head1$Dataset))
#"Swedish_National_Forest_Inventory"
na.omit(unique(head1$DUNES_BOHN))
# [1] "N" "Y"
unique(head1$ECOREG_WWF)

index3 <- which(substr(membership.conditions2,1,3)=="$$C")
membership.conditions2[index3]
# here plot.group is filled with the level number
# and the right hand side of the formula is replaced by the level number as well.
for (i in index3){
  # filling in the lef-hand side
  a <- grep(substr(membership.conditions2[i],5,nchar(membership.conditions2[i])),names(head1))
  # number of head column
  b <- as.factor(head1[,a])
  levels(b)
  #c <- as.numeric(b)
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[i])
  plot.group[,index5] <- as.numeric(b)
  
  # filling in the right-hand side
  index5 <- which(dimnames(plot.group)[[2]] %in% levels(b))
  for (j in index5){
    index7 <- match(membership.conditions2[j],levels(b))
    plot.group[,index5] <- index7
  }
}

index3 <- which(substr(membership.conditions2,1,3)=="$$N")
membership.conditions2[index3]
# here plot.group is filled with original content of the field
# there is no right-hand side of the formula, because this is numeric
for (i in index3){
  # filling in the lef-hand side
  a <- grep(substr(membership.conditions2[i],5,nchar(membership.conditions2[i])),names(head1), fixed=T)
  # number of head column
  index5 <- which(dimnames(plot.group)[[2]]==membership.conditions2[i])
  plot.group[,index5] <- head1[,a]

}


# Set of conditions 6 - this part can be faster...
Sys.time()-aq
aq<-Sys.time()

########################################################################
### Step 5: Evaluation of all conditions in a membership expression  ###
########################################################################
a <- colSums(plot.group, na.rm=T)
a[is.na(a)] # named numeric(0)
# all columns have values!
any(is.na(plot.group)) #T
plot.group[is.na(plot.group)] <- 0
any(plot.group==-Inf) #T
plot.group[which(plot.group==-Inf)]<- 0

# In the next step, the membership expressions are evaluated. This is done by
# referring to col1, col2, .... instead of membership condition names.
# For this reason, dimnames are renamed
dimnames(plot.group)[[2]] <- paste("col", seq(1:dim(plot.group)[[2]]),sep="") 

membership.expressions3 <- membership.expressions
# go through the loop of single components of the expressions,
# both left-hand and right-hand side, which is 
# in membership.conditions2 and replace them with col1, col2 etc.
# This has to be done in descending length of membership.condition2
index9 <- order(nchar(membership.conditions2), decreasing =T)

# then replace all others
for (i in 1:length(membership.conditions2)){
  membership.expressions3 <- gsub(membership.conditions2[index9][i],
                                  dimnames(plot.group)[[2]][index9][i],membership.expressions3,
                                  fixed=T)
  # dimnames in plot.group2 are col1, col2 etc.
}

# replace comparisons with R operators
membership.expressions3 <- gsub("GR",">",membership.expressions3)
membership.expressions3 <- gsub("GE",">=",membership.expressions3)
membership.expressions3 <- gsub("EQ","==",membership.expressions3)
# -> then we can use parse and eval

# the logical.matrix1 holds the results from the evaluation of the membership
# conditions in plot.group for every expression.
# The key of this step is the eval(parse(text=x)) command that allows 
# to interpret text as logical expressions in R. 
# This is still the most time-consuming part, which, however, 
# could be speeded up (e.g. by using sweep on the matrix). 
# However, I have not yet tried this out.

# -> then we can use parse and eval
lina<-data.frame(plot.group)
logical.matrix1<-foreach(i=1:length(membership.expressions3),
  .combine='cbind') %dopar% with(lina,
   eval(parse(text=membership.expressions3[i])))
dimnames(logical.matrix1)[[1]] <- plot.names
dimnames(logical.matrix1)[[2]] <- membership.expressions

# the process is very time-intensive. 
# It could be speeded up with "sweep" and "sapply"

#conditions in an expression
Sys.time()-aq
aq<-Sys.time()

#####################################################
### Step 6: Evaluation of all membership formulas ###
#####################################################
# exchange them against col1, col2, etc.
col.vector2 <- paste("col", seq(1:dim(logical.matrix1)[[2]]),sep="") 
dimnames(logical.matrix1)[[2]] <- col.vector2

# we replace the inner expressions with T and F
membership.formulas3 <- membership.formulas
index9 <- order(nchar(membership.expressions), decreasing =T)
membership.expressions[index9]
for (i in 1:length(membership.expressions)){
  membership.formulas3 <- gsub(membership.expressions[index9][i],
                               col.vector2[index9][i],membership.formulas3,
                               fixed=T)
  # dimnames in plot.group2 are col1, col2 etc.
}


# remove the angle brackets around expressions
membership.formulas3 <- gsub("<","",membership.formulas3)
membership.formulas3 <- gsub(">","",membership.formulas3)
# replace logical operators with R operators
membership.formulas3 <- gsub("AND","&",membership.formulas3)
membership.formulas3 <- gsub("OR","|",membership.formulas3)
membership.formulas3 <- gsub("NOT","&!",membership.formulas3)

lm1_df<-data.frame(logical.matrix1)
logical.matrix2 <- foreach(i=1:length(membership.formulas3),
            .combine='cbind') %dopar% with(lm1_df,
            eval(parse(text=membership.formulas3[i])))

'# from here test code
which(regexpr("col1375 &! (col737 | col1374))",membership.formulas3, fixed=T)>0)
membership.formulas3[293]
membership.formulas3[293] <- "col1375 &! (col737 | col1374)"
membership.formulas[293]
membership.formulas[293] <- "<##Q +10 D-Mires GR NON ##Q +10 D-Mires> NOT (<#TC Trees GR 15> OR <#TC Shrubs GR 15>)"

which(regexpr("GR",membership.formulas3, fixed=T)>0)
#  12 103 210 211 275 276
membership.formulas3[which(regexpr("GR",membership.formulas3, fixed=T)>0)]
membership.formulas[which(regexpr("GR",membership.formulas3, fixed=T)>0)]
membership.expressions3[which(regexpr("GR",membership.expressions3, fixed=T)>0)] # none
membership.expressions[which(regexpr("GR NON",membership.expressions, fixed=T)>0)] # 148
'
dimnames(logical.matrix2)[[2]] <- membership.formulas

# Full assignment
Sys.time()-aq
aq<-Sys.time()

'which(regexpr("##Q +03 R1Q-Inland-sanddrift-and-dune-with-siliceous-grassland",membership.formulas, fixed=T)>0)
# 12 103
membership.formulas[12]
membership.formulas[103]

length(membership.formulas)
# 298
length(membership.formulas[which(regexpr("GR",membership.formulas, fixed=T)>0)])
# 297
length(membership.formulas[which(regexpr("GR",membership.formulas3, fixed=T)>0)])
# 6, one of them:
membership.formulas[which(regexpr("GR",membership.formulas3, fixed=T)>0)]
# ##Q +03 R1Q-Inland-sanddrift-and-dune-with-siliceous-grassland
membership.formulas3[which(regexpr("GR",membership.formulas3, fixed=T)>0)]
membership.expressions[which(regexpr("GR",membership.expressions, fixed=T)>0)]
membership.expressions3[which(regexpr("GR",membership.expressions3, fixed=T)>0)] # none
membership.expressions3[which(regexpr(">",membership.expressions3, fixed=T)>0)]
'

#########################################################
### Step 7: Calculating the vector of plot assignment ###
#########################################################
a <- substr(membership.formula.names,12,length(membership.formula.names))
hier <- substr(membership.formula.names,1,1)
a <- tstrsplit(a, " ")
membership.formula.names.short <- a[[1]]
result.classification <- array("",  dim(logical.matrix2)[[1]],
                               dimnames=list(dimnames(logical.matrix2)[[1]]))
result.classification.10 <- array("",  c(dim(logical.matrix2)[[1]],10),
                                  dimnames=list(dimnames(logical.matrix2)[[1]],seq(1:10)))
for (i in 1: dim(logical.matrix2)[[1]]){
  index8 <- which(logical.matrix2[i,]==T)
  if (length(index8)==0){
    result.classification[i] <- "?"
  } else {
    freq.hier<-table(hier[index8])
    maxi.hier.leve<-max(names(freq.hier))
    freq.hier.maxi.leve<-freq.hier[maxi.hier.leve]
    
    if (freq.hier.maxi.leve==1){
      result.classification[i] <- membership.formula.names.short[index8[which(hier[index8]==maxi.hier.leve)]]  
    } else {
      result.classification[i] <- "+"
    }
    if (length(index8)>10){
      result.classification.10[i,1:10] <- membership.formula.names.short[index8[1:10]]  
    } else {
      result.classification.10[i,1:length(index8)] <- membership.formula.names.short[index8]  
    }
  }
}
result.classification[1:10]
result.classification.10[1:10,]
write.csv(result.classification, "result_classification_v3.csv", row.names = T)
write.csv(result.classification.10, "result_classification_10_v3.csv", row.names = T)

# Final time
Sys.time()-aq
Sys.time()-aq2
stopCluster(cl)

# to do
# check for $25 on the right hand side!
