#!/usr/bin/env Rscript

library(protoclust)

bench_protoclust <- function (file, h) {
    startTime <- Sys.time()
    myData <- read.csv(file, header=FALSE)  # read data
    myMatrix <- matrix(unlist(myData), ncol = dim(myData)[1], byrow=TRUE)
    myDist <- as.dist(myMatrix)
    hc <- protoclust(myDist)  # perform clustering
    cut <- protocut(hc, h=h)  # cut the dendrogram at height h, or at k clusters (use k=k keyargument)
    endTime <- Sys.time()

    duration <- as.numeric(endTime - startTime, units='secs')
    centers <- cut$protos  # get the centers of the clusters
    return(list(duration, centers))  # return the time and the number of clusters
}

cmd <- paste(commandArgs(), collapse=" ")
args <- commandArgs(trailingOnly = TRUE)
file <- unlist(args[1])
h <- as.numeric(unlist(args[2]))
results <- bench_protoclust(file, h=h)  # run the benchmark with given dataset, h and k.
duration <- results[[1]]
centers <- results[[2]]
# print the results
cat(duration, centers, sep=",") 
