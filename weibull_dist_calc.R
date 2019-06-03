#!/usr/local/bin/Rscript
library(stats)





# Fragment and Read Function
get_reads <- function(df , eta_val = 200, insertsize = c(200, 300))
{
  # Select a fragment from each transcript, size select, and return the starting position of the
  # resulting reads relative to the length of the transcript
  # Inputs: 
  #   df - dataFrame of:
  #     splice types - 0,1,2,3
  #     lengths - the lengths of the transcript
  #   eta_val - the eta value input to the Weibull distribution
  #   insertsize - a list of length two with the lower and upper bound for fragments
  # Outputs:
  #   fragments - a data frame with the start positions of the filtered fragments and the index
  #               of the transcript that it came from (columns: transcript, start)

  # sample lengths from a weibull distribution for each transcript and transform them to the length 
  # of the transcripts
  
  lengths <- df[['lengths']]
  splice_types <- df[['splice_types']]
  print(splice_types)
  deltas = log10(lengths)
  ns_minus_1 = pmax(round(lengths/eta_val/gamma(1/deltas + 1)) - 1, 0)
  xis = lapply(ns_minus_1, function(n) {diff(sort(c(runif(n), 0, 1)))})
  xis_transformed = mapply(function(x, d) {x^(1/d)}, xis, deltas, SIMPLIFY = F)
  delta_is = mapply(function(len, x_t) {round(len*x_t/sum(x_t))}, lengths, xis_transformed, SIMPLIFY = F)

  # get all the start and end points of the fragments
  starts = lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(sample(min(insertsize[1], d[1]), 1), cumsum(d[1:(length(d)-1)]))
    } else{
      sample(min(insertsize[1], d), 1)
    }
  })
  ends = lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(cumsum(d[1:(length(d)-1)]), sum(d)-sample(min(insertsize[1], sum(d) - d[length(d)]), 1))
    } else{
      d
    }
  })

  # convert to a data frame of fragments and associated transcritp index
  #fragments = data.frame(transcript = rep(1:length(deltas), lengths(delta_is)),
  fragments = data.frame(transcript = rep(1:length(deltas), unlist(lapply(delta_is, length))),
                         start = unlist(starts),
                         end = unlist(ends))
  fragments$length = fragments$end - fragments$start

  # Filter fragments by length and return
  fragments = fragments[fragments$length >= insertsize[1] & fragments$length <= insertsize[2],]
  return(fragments[c('transcript', 'start')])
}



# read in info (as from command line)
## args: lengths, splice types (as 2 separate arrays)
args = (commandArgs(trailingOnly=TRUE))
print("running weibull dist")

# ####################
# # FOR TESTING: ##
# eta_val = 200
# insertsize = c(200, 300)
# args <-  "[[3800, 3406, 3800, 1303, 3800, 3800, 3286, 3800, 3800, 3800, 3800, 1310, 3800, 3800, 3800, 3800, 3179, 3800, 3800, 3800, 3559, 3800, 3800, 3800, 3800, 3692, 3800, 3800, 3800, 3800, 3800, 3523, 1067, 3800, 2364, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3417, 3800, 2557, 3800, 3800, 3800, 3800, 3800, 3800, 2227, 1886, 3800, 3800, 3800, 3800, 3800, 3008, 2228, 3800, 1099, 3800, 3800, 3800, 3800, 3800, 3800, 1170, 3800, 3800, 3800, 3800, 2412, 3800, 3800, 3800, 1098, 3579, 3800, 3800, 3800, 3800, 3800, 2694, 3800, 1212, 3800, 3800, 3800, 2233, 3800, 3800, 2590, 1202, 3800, 3004, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3348, 3800, 3513, 3800, 1008, 3800, 3800, 3613, 3800, 2223, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 1707, 1751, 3800, 3800, 3034, 3800, 3800, 1572, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3456, 3800, 3800, 3751, 3800, 3800, 1080, 3800, 3800, 3800, 3800, 3234, 3298, 2554, 3800, 3800, 3800, 3465, 1828, 3800, 3800, 3800, 2464, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 2499, 3800, 3800, 3800, 3800, 3800, 1428, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3134, 3800, 2878, 1331, 1641, 2263, 3470, 3800, 3800, 3800, 3800, 3800, 3800, 1916, 2166, 3800, 1833, 3800, 3800, 3800, 1001, 3800, 3800, 3800, 3718, 3800, 3800, 3340, 3800, 2952, 3800, 1007, 3800, 3800, 3800, 3800, 3800, 1532, 1295, 3800, 3800, 1028, 3800, 2214, 2603, 1450, 3800, 3800, 2591, 3800, 3800, 3800, 3800, 2978, 3800, 3643, 2457, 3800, 3800, 2176, 3027, 3800, 3800, 3800, 3800, 2647, 3507, 2018, 3436, 3800, 2579, 1927, 3011, 3800, 3800, 3800, 3800, 3800, 2275, 3800, 3800, 3800, 3800, 3800, 2974, 1588, 3800, 2142, 3800, 3382, 3800, 1601, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 1198, 3800, 3800, 3800, 3800, 2220, 2551, 3076, 3800, 1401, 3800, 3800, 3800, 3800, 2071, 3800, 2677, 2432, 1393, 2048, 3638, 3800, 3800, 3800, 3732, 2552, 3800, 3800, 3800, 2786, 3800, 3800, 1753, 3800, 3800, 3474, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 2168, 3532, 3264, 3800, 3800, 3800, 3800, 2416, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 2401, 3800, 3800, 3800, 3800, 3800, 3800, 2594, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 3800, 1159, 3800, 1550, 1006, 3800, 3800, 3800, 2503, 2139, 1864, 3800, 3800, 1769, 3800, 3315, 3058], 
#             [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0]]"
# ####################


# Passing in an array of lengths and splice types from python script
q <- gsub("\\[|\\]" ,"", unlist(strsplit(args, split="]"))) # splits input into separate lists

lengths <-as.numeric(unlist(strsplit(q[1], split=", ")))
splice_types <-(as.numeric(unlist(strsplit(q[2], split=", ")))) ## TODO: need to include these in analysis
splice_types <- splice_types[!is.na(splice_types)]
df <- data.frame(splice_types = splice_types,lengths=lengths)
start_pos = get_reads(df)

print(paste(gsub("\n","",as.character(start_pos[1]))," ",gsub("\n","",as.character(start_pos[2]))))
















##  OLD ********###############################
# read in info (as from command line)
## args: eta_val, length
#args = (commandArgs(trailingOnly=TRUE))
#eta_val = as.numeric(args[1])
# NOTE: for simplicity treated lengths  here as a vector, but really just passing in a single value from the .py script
#lengths = c(as.numeric(args[2]))





#deltas = log10(lengths)
#ns_minus_1 = pmax(round(lengths/eta_val/gamma(1/deltas + 1)) - 1, 0)
#xis = lapply(ns_minus_1, function(n) {diff(sort(c(runif(n), 0, 1)))})
#xis_transformed = mapply(function(x, d) {x^(1/d)}, xis, deltas, SIMPLIFY = F)
#delta_is = mapply(function(len, x_t) {round(len*x_t/sum(x_t))}, lengths, xis_transformed, SIMPLIFY = F)
#print(as.character(delta_is))

































