#!/usr/local/bin/Rscript
library(stats)

# read in info (as from command line)
## args: d_prime, r_prime, transcription_rate (as a list in python)
args = (commandArgs(trailingOnly=TRUE))
print("running sum square equation solver")

# ####################
# # FOR TESTING: ##
# args <- "[0.7333333333333333, 0.0002881702719330355, 1500.0]"
# ####################


# Passing in an array of lengths and splice types from python script
q <- gsub("\\[|\\]" ,"", unlist(strsplit(args, split="]"))) # splits input into separate lists

q <-(as.numeric(unlist(strsplit(q[1], split=", ")))) ## TODO: need to include these in analysis
# Root solver function
D_prime = q[1]
R_prime = q[2]
txnrate = q[3]
hold.row <- c(NA, NA)
f <- function(h){ ((h*(1 - 2^(-D_prime[1]/(h*txnrate)))) - R_prime[1])^2}# + ((h*(1 - 2^(-D_prime[1]/(h*txnrate)))) - R_prime[1])^2 + ((h*(1 - 2^(-D_prime[1]/(h*txnrate)))) - R_prime[1])^2 }
starth = 0
if(sum(is.na(R_prime))==3){ return(hold.row) }
try(fit.hold <- optim(starth, f))
try(hold.row <- c(fit.hold$par, fit.hold$value))

output = hold.row # root, zerodev


print(paste(gsub("\n","",as.character(output[1]))," ",gsub("\n","",as.character(output[2]))))

