source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("affy")
library(affy)
library(stringr) 
library(caret)

# untar the raw data file downloaded from NCBI omnibus
untar("GSE40967_RAW.tar", exdir = "GSE40967_RAW")
#  move the files into a subfolder named "GSE40967_RAW" to reduce clutter

# then read the data, the data in the folder is compressed, so set compress option to TRUE
affy.data = ReadAffy(celfile.path = "GSE40967_RAW", compress = T)

# convert aff.data to instance of ExpressionSet with MAS 5.0 method
eset.mas5 = mas5(affy.data)

# retrieve matrix of expression values
exprSet.nologs = exprs(eset.mas5)

# modify sample names so the only thing icluded is GSM#
samps = colnames(exprSet.nologs)
new_colnms = str_sub(samps, 1, 10)
new_colnms = sub("_","",new_colnms)
colnames(exprSet.nologs) = new_colnms
colon_data = data.frame(t(exprSet.nologs))
colnames(colon_data) = str_sub(colnames(colon_data),2,str_length(colnames(colon_data)))

##############################################
# Assemble class file
##############################################
path = "GSE39582_header.txt"
con = file(path, "r")
data = read.delim(con, header = F, sep = "!")

# these are the entries where sample name and the corresponing class are found for the colon dataset
# editing will be needed here for different datasets
names = toString(data$V2[50])
classes = toString(data$V2[88])

# split these entries into lists then put into a table so that we have each sample name with the corresponding class
line = unlist(strsplit(names,"\t"))
line2 = unlist(strsplit(classes,"\t"))
line = line[-1]
line2 = line2[-1]
line2 = str_sub(line2, str_length(line2), str_length(line2))
rms = which(line2=="A")
class_table = data.frame(sample = line[-rms], class = factor(line2[-rms]))

# write this table out so it can be examined seperately or read in seperately should class_table be delete/corrupted
write.table(class_table, quote = F, file = "colon_classes.txt", sep = "\t", row.names = F, col.names = F)

# only hold onto data that there is information on the class of
keep = which(rownames(colon_data) %in% class_table$sample)
colon_data = colon_data[keep,]



##################################################################################
# the following section of code gets rid of genes that are detected in fewer than 20% of all samples
# there were no probes that met this criteria for this data set
# this step will more often be necessary for RNA-seq data
##################################################################################
lowDetects = c()
for (i in 1:length(colon_data)){ 
  colon_data[,i] = as.numeric(colon_data[,i])
  if (sum(colon_data[,i]!=0) < (0.2 * length(colon_data[,i]))){
    lowDetects = c(lowDetects,i)
  }
}
colon_data = colon_data[,-lowDetects]

#################################################################
# Transform the data and write to file
#################################################################

# use log base 2 so that an increase of 1 in the data represents a 2-fold change in expression level
# 2 fold change is typically considered biologically relevant
# add 1 to the raw expression values so that the minimum transformed expression value is 0
colon_data = log(colon_data + 1, 2)
to_export = t(colon_data)

write.table(to_export,"CC_mas5_matrix_colon.txt", quote=F, sep = "\t")
rm(to_export)

#################################################################
# now gene-probe mapping
#################################################################
biocLite("hgu133plus2.db") # this is platform specific and will need adjustment depending on platform of dataset
library(hgu133plus2.db)
# create the mapper in a dataframe
mapper = as.matrix(unlist(as.list(hgu133plus2SYMBOL)))
mapper[,2] = rownames(mapper)
rownames(mapper) = NULL
colnames(mapper) = c("entrez","probe")

# remove any probes from the dataset that are not in the mapper
keep = which(colnames(colon_data) %in% mapper$probe)
keep = factor(mapper$probe[keep])
keep = which(colnames(colon_data) %in% keep)
colon_data = colon_data[,keep]

# remove the probes that do not have an associated entrez gene from both the data and the mapper table
rms = which(is.na(mapper$entrez))
rms = factor(mapper$probe[rms])
leaveOut = which(colnames(colon_data) %in% rms)
leaveOutMapper = which(mapper$probe %in% colnames(colon_data)[leaveOut])
colon_data = colon_data[,-leaveOut]
mapper = mapper[-leaveOutMapper,]
mapper = mapper[-which(is.na(mapper$entrez)),]
keepm = which(mapper$probe %in% colnames(colon_data))
mapper = mapper[keepm,]


################################################################################
# find the max value of each gene for each sample then condense in python
################################################################################
write.csv(file = "colon_data.csv",x =colon_data)
write.csv(file = "mapper.csv",x = mapper)

# there are multiple probes corresponding to each gene, all with different expression values
# find the max values for each gene with python script 'col.py'
# make sure findMax.py is in the working directory
system("python findMax.py") 

new_data = read.csv("mod_colon.csv",header = T, sep = ",", row.names = rownames(colon_data))
colnames(new_data) = sub("X", "",colnames(new_data))

rownames(new_data) = rownames(colon_data)

new_data = new_data[,-1]


old_colon = colon_data
colon_data = new_data
rm(new_data)


##########################################################################
# create 10-folds
##########################################################################
# partitions the data into 10 folds and creates the training sets and test sets
target1 <- rownames(colon_data)
class_table = class_table[match(target1, class_table$sample),]
Class = class_table
rownames(Class) = Class$sample
Class$sample = NULL
flds <- createFolds(Class$class, k = 10, list = F)
samp_range = 1:length(rownames(colon_data))
it = c("01","02","03","04","05","06","07","08","09","10")
j = 1
for (i in it){
  nameTrain <- paste("train",i,sep = "")
  nameTest <- paste("test",i,sep="")
  assign(nameTest, colon_data[which(flds == j),])
  assign(nameTrain, colon_data[which(flds != j),])
  j = j+1
}

##########################################################################
# Transform trains to z-scores
##########################################################################
for (i in 1:length(it)){
  i = 1
  nameTrain <- paste("train",it[i],sep = "")
  nameTest <- paste("test",it[i],sep="")
  dat = get(nameTrain)
  datT = get(nameTest)
  for (j in 1:length(colnames(dat))){
    std = sd(dat[,j])
    mu = mean(dat[,j])
    dat[,j] = (dat[,j] - mu)/std
    datT[,j] = (datT[,j] - mu)/std
  }
  assign(nameTrain, dat)
  assign(nameTest, datT)

}

#################################################################
# Feature selection
#################################################################
#library(pamr)
#pam.data = pamr.from.excel("CC_mas5_matrix_colon.txt",825, sample.labels = FALSE)
library(FSelector)

# The filter methods below are done split into partitions and then the significance statistics are calculated for each gene in each partition because doing all of the genes at once uses a lot of memory
# since the features are independent in univariate filtering, the results are the same as they would be if the whole dataset was done at once
########################
# Chi-squared FS
########################
x2weights = data.frame()
i = 1
inc = 1000
while (i<length(colnames(colon_data))){
  subset = colon_data[,i:(i+min(inc - 1,length(colnames(colon_data))-i))]
  subset$Class = class_table$class
  subset = subset[,c(length(subset),1:(length(subset)-1))]
  x2weights = rbind(x2weights, (FSelector::chi.squared(Class ~., subset)))
  i = i + inc
}


target <- c(order(x2weights$attr_importance, decreasing = T))
x2_ranking = data.frame(whole_dataset = rownames(x2weights)[target])


k = 1
for (j in it){
  nameTrain <- paste("train",j,sep = "")

  dat = get(nameTrain)
  thisClasses = Class[which(flds != k),]
  x2weights = data.frame()
  i = 1
  inc = 2000
  while (i<length(colnames(dat))){
    subset = dat[,i:(i+min(inc - 1,length(colnames(dat))-i))]
    subset$Class = thisClasses
    subset = subset[,c(length(subset),1:(length(subset)-1))]
    x2weights = rbind(x2weights, (FSelector::chi.squared(Class ~., subset)))
    i = i + inc
  }
  cat("Completetd X2 for" , nameTrain, "\n")
  target <- c(order(x2weights$attr_importance, decreasing = T))

  x2_ranking[,nameTrain] = rownames(x2weights)[target]
  k = k + 1
}


write.csv(x2_ranking, file = "chi_squared_genename.csv")
##############################
# Symmetrical Uncertainty FS
##############################
su.weights = data.frame()
i = 1
inc = 1000
while (i<length(colnames(colon_data))){
  subset = colon_data[,i:(i+min(inc - 1,length(colnames(colon_data))-i))]
  subset$Class = class_table$class
  subset = subset[,c(length(subset),1:(length(subset)-1))]
  su.weights = rbind(su.weights, (FSelector::symmetrical.uncertainty(Class ~., subset)))
  i = i + inc
}


target <- c(order(su.weights$attr_importance, decreasing = T))
su_ranking = data.frame(whole_dataset = rownames(su.weights)[target])


k = 1
for (j in it){
  nameTrain <- paste("train",j,sep = "")

  dat = get(nameTrain)
  thisClasses = Class[which(flds != k),]
  su.weights = data.frame()
  i = 1
  inc = 1000
  while (i<length(colnames(dat))){
    subset = dat[,i:(i+min(inc - 1,length(colnames(dat))-i))]
    subset$Class = thisClasses
    subset = subset[,c(length(subset),1:(length(subset)-1))]
    su.weights = rbind(su.weights, (FSelector::symmetrical.uncertainty(Class ~., subset)))
    i = i + inc
  }
  cat("Completetd SU for" , nameTrain, "\n")
  target <- c(order(su.weights$attr_importance, decreasing = T))

  su_ranking[,nameTrain] = rownames(su.weights)[target]
  k = k + 1
}

write.csv(su_ranking, file = "SU_genename.csv")

##############################
# Information Gain FS
##############################
ig.weights = data.frame()
i = 1
inc = 1000
while (i<length(colnames(colon_data))){
  subset = colon_data[,i:(i+min(inc - 1,length(colnames(colon_data))-i))]
  subset$Class = class_table$class
  subset = subset[,c(length(subset),1:(length(subset)-1))]
  ig.weights = rbind(ig.weights, (FSelector::information.gain(Class ~., subset)))
  i = i + inc
}


target <- c(order(ig.weights$attr_importance, decreasing = T))
ig_ranking = data.frame(whole_dataset = rownames(ig.weights)[target])


k = 1
for (j in it){
  nameTrain <- paste("train",j,sep = "")

  dat = get(nameTrain)
  thisClasses = Class[which(flds != k),]
  ig.weights = data.frame()
  i = 1
  inc = 1000
  while (i<length(colnames(dat))){
    subset = dat[,i:(i+min(inc - 1,length(colnames(dat))-i))]
    subset$Class = thisClasses
    subset = subset[,c(length(subset),1:(length(subset)-1))]
    ig.weights = rbind(ig.weights, (FSelector::information.gain(Class ~., subset)))
    i = i + inc
  }
  cat("Completetd IG for" , nameTrain, "\n")
  target <- c(order(ig.weights$attr_importance, decreasing = T))

  ig_ranking[,nameTrain] = rownames(ig.weights)[target]
  k = k + 1
}

write.csv(ig_ranking, file = "IG_genename.csv")



#################
# PAM FS
#################
library(pamr)
pam.data = list(x=t(data.matrix(colon_data, rownames.force = F)),y=class_table$class,geneids=colnames(colon_data),genenames=colnames(colon_data))
pam.train = pamr.train(pam.data)
genelist = pamr.listgenes(pam.train, pam.data,1,genenames=T)
wd = genelist[1:300,2]

pam_ranking = data.frame(whole_dataset = wd)

k = 1
for (j in it){
  nameTrain <- paste("train",j,sep = "")

  dat = get(nameTrain)
  thisClasses = Class[which(flds != k),]

  pam.data = list(x=t(data.matrix(dat, rownames.force = F)),y=thisClasses,geneids=mapper$entrez,genenames=colnames(colon_data))
  pam.train = pamr.train(pam.data, n.threshold=100)
  genelist = pamr.listgenes(pam.train, pam.data,1,genenames=T)

  cat("Completetd PAM for" , nameTrain, "\n")


  pam_ranking[,nameTrain] = genelist[1:300,2]
  k = k + 1
}


write.csv(pam_ranking, file = "PAM_genename.csv")


############################################################
# MEMORY MANAGEMENT
############################################################
# At this point there are a lot of extra variable hanging out in the workspace, remove the unneccesary ones to free up RAM
# the code below is used to cut the the training and test sets down to only the genes deemed neccesary by feature selection

distinct_genes = as.vector(genes_for_train)[which(!duplicated(unlist(as.vector(genes_for_train))))]

for(t in 1:10){
  trainer = get(trains[t])
  tester = get(tests[t])
  trainer = trainer[,match(distinct_genes,colnames(trainer))]
  tester = tester[,match(distinct_genes,colnames(tester))]
  assign(trains[t],trainer)
  assign(tests[t],tester)
}
rm(trainer)
rm(tester)

# after executing the above code, save the workspace, restart R studio, and proceed to PNBC.R














