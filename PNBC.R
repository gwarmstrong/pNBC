# set working directory and load workspace
# it can free up a lot of RAM by rm'ing unneccesary varaibles after preprocessing, then saving the workspace,
# then closing and reopening Rstudio and loading the workspace where you saved the neccessary data

# setwd("filepath/of/workdir")
# load("filepath/of/workspace") 

# these flags are needed to compile the cpp package properly
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

#install.packages(c("Rcpp","stringr","caret","rlist","minet","RWeka","parallel","pryr"))
library(Rcpp)
library(stringr)
library(caret)
library(rlist)
library(minet)
library(RWeka)
library(parallel)
library(pryr)

# compile the cpp source code
sourceCpp("rNBC.cpp",verbose = F)


# initialize for later use
trains = c("train01","train02","train03","train04","train05","train06","train07","train08","train09","train10")
tests = c("test01","test02","test03","test04","test05","test06","test07","test08","test09","test10")
fs_gene = list(x2_ranking[,2:11],su_ranking[,2:11], ig_ranking[,2:11],pam_ranking[,2:11])

# loads the maximum number of features you intend on using to genes_for_train
max_feat_num = 300
num_fs = 4
genes_for_train = array(NA, c(num_fs,length(trains),max_feat_num))
for(fs in 1:num_fs){
  for(t in 1:length(trains)){
    genes_for_train[fs,t,] = fs_gene[[fs]][1:max_feat_num,t]
  }
}

rm(fs_gene)




############################################################################
# Functions for classification
############################################################################

# use the regression coefficients supplied by reg to predict all of the values for (samples x genes) matrix data
wiClass.nbc <- function(data, reg){
  prelim <- cbind(as.matrix(data),rep(1,nrow(data))) %*% reg
  return(cbind(as.matrix(data),rep(1,nrow(data))) %*% reg)
}

# compute and return the squared residuals for of a given set of observed and predicted values of the same size
resid <- function(samples, pred.profile){
  res <- samples - pred.profile
  return(abs(res)^2)
}

# returns the index of the lowest value, useful for determining the class of a sample in predict.nbc
classOfMin <- function(vec){
  which(vec == min(vec))
}

# returns the set of predictions for each sample in data
# samples are predicted to be in the class with the lowest average residuals
predict.nbc <- function(nbc, data, num_classes, nrn){
  ns = length(rownames(data))
  predsByRow = data.frame(matrix(ncol = ns, nrow = 0))
  resMinFinder = vector(mode="list", length = num_classes)
  for (c in 1:num_classes){
    resForClass = data.frame(matrix(ncol = ns, nrow = nrn))
    for(n in 1:nrn){
      pred_vals <- wiClass.nbc(as.matrix(data), nbc[[c]][[n]])
      res_vals <- resid(as.matrix(data), pred_vals)
      resForClass[n,] = rowSums(res_vals)
    }
    resMinFinder[[c]] = resForClass
  }
  
  min_for_class = data.frame(matrix(ncol = ns, nrow = num_classes))
  for (c in 1:num_classes){
    min_for_class[c,] = apply(resMinFinder[[c]], 2, mean)
  }
  
  
  predictions = as.vector(apply(min_for_class, 2, classOfMin))
  return(predictions)
}

# returns the set of predictions for each sample in data when the residuals of gene[i] should be weighted by weight[i]
# chooses the class with lowest average weighted residuals
predictWeights.nbc <- function(nbc, data, num_classes, nrn, weights){
  if(length(weights)!= length(colnames(data))){
    warning("incorrect number of weights")
  }
  ns = length(rownames(data))
  predsByRow = data.frame(matrix(ncol = ns, nrow = 0))
  resMinFinder = vector(mode="list", length = num_classes)
  for (c in 1:num_classes){
    resForClass = data.frame(matrix(ncol = ns, nrow = nrn))
    for(n in 1:nrn){
      pred_vals <- wiClass.nbc(as.matrix(data), nbc[[c]][[n]])
      res_vals <- resid(as.matrix(data), pred_vals)
      resForClass[n,] = res_vals %*% as.matrix(weights,ncol=1)
    }
    resMinFinder[[c]] = resForClass
  }
  
  min_for_class = data.frame(matrix(ncol = ns, nrow = num_classes))
  for (c in 1:num_classes){
    min_for_class[c,] = apply(resMinFinder[[c]], 2, mean)
  }
  
  
  predictions = as.vector(apply(min_for_class, 2, classOfMin))
  return(predictions)
}

# returns the set of predictions for each sample in data when only one network is saved for each class
predictSingle.nbc <- function(nbc, data, num_classes){
  ns = length(rownames(data))
  predsByRow = data.frame(matrix(ncol = ns, nrow = 0))
  resMinFinder = vector(mode="list", length = num_classes)
  for (c in 1:num_classes){
    resForClass = data.frame(matrix(ncol = ns, nrow = 1))
    pred_vals <- wiClass.nbc(as.matrix(data), nbc[[c]])
    res_vals <- resid(as.matrix(data), pred_vals)
    resForClass = rowSums(res_vals)
    resMinFinder[[c]] = resForClass
  }
  
  
  min_for_class = data.frame(matrix(ncol = ns, nrow = num_classes))
  for (c in 1:num_classes){
    min_for_class[c,] = resMinFinder[[c]]
  }
  predictions = as.vector(apply(min_for_class, 2, classOfMin))
  return(predictions)
}

# trains nrn probabilistic networks for each class, mim_thresh and dpi_thresh are the thresholds for the aracne construction. part and classes are the training data and class labels
nbc.train <- function(part, classes, nrn, mim_thresh, dpi_thresh){
  trimmed = part
  classTrim = classes
  edge_list = vector(mode = "list", length = num_classes)
  results = vector(mode = "list", length = num_classes)
  
  # construct gene regulatory network with aracne for the entire training data
  mimtrim <- build.mim(trimmed)
  mimtrim[which(mimtrim<mim_thresh)] = 0
  nettrim = aracne(mimtrim, dpi_thresh)
  rm(mimtrim)
  low_net <- nettrim
  rm(nettrim)
  low_net[upper.tri(low_net, diag = T)] = 0
  connections = which(low_net>0, arr.ind = T)
  edges = data.frame(matrix(ncol = 2, nrow = length(rownames(connections))))
  edges[,1:2] = connections[,2:1]
  colnames(edges) = c("gene_a","gene_b")
  
  
  # calculates the probability (pearson correlation) that each edge exists in each class
  for (c in 1:num_classes){
    new_col = c()
    for (i in 1:length(edges$gene_a)){
      new_col = c(new_col,cor(trimmed[which((classTrim==c)),edges$gene_a[i]],trimmed[which((classTrim==c)),edges$gene_b[i]]))
    }
    edges = cbind(edges,new_col)
    colnames(edges)[c+2] = paste("c",c)
  }
  
  results = vector(mode = "list", length = num_classes)
  
  # uses the graphs from above to train regressions for the nrn probabilistic networks
  for (c in 1:num_classes){
    results[[c]]<- nbcTrainN(as.matrix(edges[,c(1,2,c+2)]),as.matrix(trimmed[which((classTrim==c)),]),nrn)
  }
  
  
  return(results)
}


# generates nrn probabilistic networks, based on the aracne constructed network, for each class and returns then uses the 
# cpp package to train the networks, and returns the regression coefficients on each network
# uses a distinct aracne network for each class
nbc.trainCAracne <- function(part, classes, nrn, mim_thresh, dpi_thresh) {
  trimmed = part
  classTrim = classes
  edge_list = vector(mode = "list", length = num_classes)
  results = vector(mode = "list", length = num_classes)
  
  # constructs aracne networks for each class then calculates the edge probabilities
  for (c in 1:length(levels(factor(classTrim)))){
    trim_for_class = trimmed[which((classTrim==c)),]
    
    # do aracne network construction
    mimtrim <- build.mim(trim_for_class)
    mimtrim[which(mimtrim<mim_thresh)] = 0
    nettrim = aracne(mimtrim, dpi_thresh)
    low_net <- nettrim
    low_net[upper.tri(low_net, diag = T)] = 0
    connections = which(low_net>0, arr.ind = T)
    
    edges = data.frame(matrix(ncol = 2, nrow = length(rownames(connections))))
    if(length(connections>0)){
      edges[,1:2] = connections[,2:1]
    }
    colnames(edges) = c("gene_a","gene_b")
    
    # calculates edge probabilities based on pearson correlation
    new_col = vector(mode = "numeric", length = 0)
    if(length(rownames(edges))>0){
      for (i in 1:length(edges$gene_a)){
        new_col = c(new_col,abs(cor(trim_for_class[,edges$gene_a[i]],trim_for_class[,edges$gene_b[i]])))
      }
    }
    edges = cbind(edges,new_col)
    colnames(edges)[3] = paste("c",c)
    edge_list[[c]] = edges
  }
  
  rm(mimtrim)
  rm(nettrim)
  
  # trains the networks using the cpp package
  for (c in 1:num_classes){
    results[[c]]<- nbcTrainN(as.matrix(edge_list[[c]]),as.matrix(trimmed[which((classTrim==c)),]),nrn)
  }
  
  return(results)
}

# trains pnbc classifier, generating nrn random networks, but only keeps the best performing network from each class
nbc.trainSingle <- function(part, classes, nrn, mim_thresh, dpi_thresh) {
  data = part
  #get the trained networks from nbc.train
  nbc = nbc.train(part, classes, nrn, mim_thresh, dpi_thresh)
  
  results = vector(mode = "list", length = num_classes)
  
  ns = length(rownames(data))
  predsByRow = data.frame(matrix(ncol = ns, nrow = 0))
  resMinFinder = vector(mode="list", length = num_classes)
  net_to_keep = vector(mode = "numeric", length = num_classes)
  # find the network with the lowest average residual for each class and keep this network
  for (c in 1:num_classes){
    resForClass = data.frame(matrix(ncol = ns, nrow = nrn))
    for(n in 1:nrn){
      pred_vals <- wiClass.nbc(as.matrix(data), nbc[[c]][[n]])
      res_vals <- resid(as.matrix(data), pred_vals)
      resForClass[n,] = rowSums(res_vals)
      
    }
    net_to_keep[c] = which(rowSums(resForClass) == min(rowSums(resForClass)))[1]
    results[[c]] = nbc[[c]][[net_to_keep[c]]]
  }
  return(results)
}



