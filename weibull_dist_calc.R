#!/usr/local/bin/Rscript
## args: eta_val, length
args = (commandArgs(trailingOnly=TRUE))
# read in info (as from command line)
eta_val = as.numeric(args[1])
length = as.numeric(args[2])






deltas = log10(lengths)
ns_minus_1 = pmax(round(lengths/eta_val/gamma(1/deltas + 1)) - 1, 0)
xis = lapply(ns_minus_1, function(n) {diff(sort(c(runif(n), 0, 1)))})
xis_transformed = mapply(function(x, d) {x^(1/d)}, xis, deltas, SIMPLIFY = F)
delta_is = mapply(function(len, x_t) {round(len*x_t/sum(x_t))}, lengths, xis_transformed, SIMPLIFY = F)












































statistic<-read.table(paste('output_data/',fname,"_top",args[1],"_bot",args[2],"_test_stats.csv",sep = ""),header=TRUE, sep=",", dec=".")
#Naming each column to analyze separately
A= statistic [,1]
U= statistic [,2]
C= statistic [,3]
G= statistic [,4]
#G= statistic [,3]
#C= statistic [,4]
#Getting pvlaues for each base with indicated degrees of freedom (passed in from command line)
pvalueA=2*pt(-abs(A),args[3])
pvalueU=2*pt(-abs(U),args[3])
pvalueC=2*pt(-abs(C),args[3])
pvalueG=2*pt(-abs(G),args[3])
#Combining all pvalues into same data frame in order to export together 
finalpvalue<-data.frame(pvalueA,pvalueU,pvalueC,pvalueG)
#Export data frame of pvalues into csv file
#write.table(finalpvalue, file="output_data/pvalue.csv",sep=",",row.names=F)
write.table(finalpvalue, file=paste("output_data/",fname,"_top",args[1],"_bot",args[2],"_pvalues.csv",sep = ""),sep=",",row.names=F)
paste(args[3], "degrees of freedom")
paste("Data written to output_data/",fname,"_top",args[1],"_bot",args[2],"_pvalues.csv",sep = "")

