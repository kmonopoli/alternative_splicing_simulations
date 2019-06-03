#!/usr/local/bin/Rscript
library(stats)





# Fragment and Read Function
fragment_and_read <- function(df , eta_val = 200, insertsize = c(200, 300))
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
#args <- "[[5946, 1475, 500, 500, 5741, 500, 500, 500, 500, 500, 500, 500, 1602, 4574, 5270, 5723, 500, 4889, 3233, 5263, 500, 500, 500, 500, 500, 3351, 500, 500, 2474, 500, 500, 500, 500, 500, 500, 2804, 500, 500, 5620, 500, 1860, 5587, 500, 5449, 500, 2758, 4463, 500, 3214, 1011, 500, 500, 500, 500, 500, 500, 2186, 500, 1170, 500, 500, 500, 500, 5585, 500, 500, 500, 500, 1946, 1027, 6106, 5222, 500, 500, 2795, 1041, 500, 500, 500, 500, 500, 3680, 500, 500, 500, 500, 500, 5661, 500, 500, 500, 3543, 3985, 500, 500, 500, 500, 500, 3434, 500, 500, 500, 500, 500, 500, 500, 1942, 500, 500, 500, 500, 4007, 500, 500, 500, 500, 500, 5231, 500, 500, 500, 5903, 500, 500, 1413, 3458, 1532, 500, 500, 1051, 500, 500, 500, 5177, 500, 3572, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 4499, 500, 500, 500, 500, 500, 5960, 500, 500, 500, 500, 500, 500, 3221, 500, 4304, 1397, 500, 500, 2388, 500, 500, 500, 500, 500, 1488, 500, 500, 5917, 5209, 500, 4012, 500, 500, 2084, 1690, 500, 1877, 2263, 500, 5317, 500, 500, 2030, 500, 500, 500, 500, 500, 5031, 500, 3811, 2227, 500, 500, 500, 3787, 500, 500, 2121, 1422, 500, 500, 2141, 500, 500, 500, 1706, 500, 4424, 500, 500, 500, 2157, 5779, 500, 3861, 500, 500, 5847, 500, 500, 500, 5950, 500, 500, 500, 500, 500, 500, 3614, 500, 5188, 500, 500, 500, 500, 500, 500, 2165, 4589, 500, 5301, 500, 500, 500, 500, 3369, 500, 500, 1275, 1264, 500, 5566, 500, 500, 5051, 1052, 500, 500, 500, 500, 500, 1345, 500, 500, 4227, 500, 3205, 2351, 2667, 4936, 500, 500, 500, 5852, 500, 500, 3000, 1368, 5781, 1150, 500, 500, 1386, 500, 3838, 500, 500, 3027, 500, 500, 500, 4518, 3844, 500, 500, 500, 500, 3676, 500, 500, 1485, 500, 500, 1831, 3859, 500, 500, 500, 3151, 500, 500, 1687, 2889, 500, 5850, 5721, 4461, 1582, 500, 500, 4837, 2272, 500, 500, 500, 500, 500, 3837, 500, 2299, 500, 500, 1909, 500, 4583, 500, 4105, 2461, 500, 500, 500, 2230, 1022, 500, 5865, 500, 5033, 4626, 500, 500, 500, 500, 5273, 500, 500, 500, 1065, 500, 500, 500, 500, 500, 4846, 500, 500, 500, 1439, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 1871, 500, 500, 500, 500, 500, 500, 500, 2528, 1792, 3705], [0, 4, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2, 0, 2, 2, 4, 2, 2, 2, 2, 2, 2, 4, 2, 2, 0, 2, 1, 0, 2, 0, 2, 4, 0, 2, 0, 4, 2, 2, 2, 2, 2, 2, 4, 2, 4, 2, 2, 2, 2, 0, 2, 2, 2, 2, 4, 4, 0, 0, 2, 2, 4, 4, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 2, 2, 2, 0, 0, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 4, 0, 4, 2, 2, 4, 2, 2, 2, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 0, 2, 0, 4, 2, 2, 4, 2, 2, 2, 2, 2, 4, 2, 2, 0, 0, 2, 0, 2, 2, 4, 4, 2, 4, 4, 2, 0, 2, 2, 3, 2, 2, 2, 2, 2, 0, 2, 0, 4, 2, 2, 2, 0, 2, 2, 3, 4, 2, 2, 4, 2, 2, 2, 4, 2, 0, 2, 2, 2, 4, 0, 2, 0, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 0, 2, 0, 2, 2, 2, 2, 2, 2, 4, 0, 2, 0, 2, 2, 2, 2, 0, 2, 2, 4, 4, 2, 0, 2, 2, 0, 4, 2, 2, 2, 2, 2, 4, 2, 2, 0, 2, 0, 4, 4, 0, 2, 2, 2, 0, 2, 2, 0, 4, 0, 4, 2, 2, 4, 2, 0, 2, 2, 0, 2, 2, 2, 0, 0, 2, 2, 2, 2, 0, 2, 2, 4, 2, 2, 4, 0, 2, 2, 2, 0, 2, 2, 4, 4, 2, 0, 0, 0, 4, 2, 2, 0, 4, 2, 2, 2, 2, 2, 0, 2, 4, 2, 2, 4, 2, 0, 2, 0, 4, 2, 2, 2, 4, 4, 2, 0, 2, 0, 0, 2, 2, 2, 2, 0, 2, 2, 2, 3, 2, 2, 2, 2, 2, 0, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 4, 4, 0]]"
# ####################


# Passing in an array of lengths and splice types from python script
q <- gsub("\\[|\\]" ,"", unlist(strsplit(args, split="]"))) # splits input into separate lists

lengths <-as.numeric(unlist(strsplit(q[1], split=", ")))
splice_types <-(as.numeric(unlist(strsplit(q[2], split=", ")))) ## TODO: need to include these in analysis
splice_types <- splice_types[!is.na(splice_types)]
df <- data.frame(splice_types = splice_types,lengths=lengths)
start_pos = fragment_and_read(df)
print(start_pos)
print(paste(gsub("\n","",as.character(start_pos[1]))," ",gsub("\n","",as.character(start_pos[2]))))

