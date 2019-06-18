#!/usr/local/bin/Rscript
library(stats)





# Fragment and Read Function
fragment_and_read <- function(df , eta_val = 200, insertsize = c(200, 300))
{
  # Select a fragment from each transcript, size select, and return the starting position of the
  # resulting reads relative to the length of the transcript
  # Inputs: 
  #   df - dataFrame of:
  #     splice types - 1,2,3,4,5
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
#  eta_val = 200
#  insertsize = c(200, 300)
# args <- "[[22586, 16995, 12689, 17029, 18694, 11771, 11877, 25911, 10357, 10272, 11012, 25225, 21933, 16445, 20039, 20105, 16601, 18602, 19930, 12307, 10206, 27567, 22559, 15235, 11775, 20848, 11807, 29593, 31846, 10882, 10338, 13904, 17926, 10923, 11190, 16691, 12726, 13877, 11555, 29699, 11353, 24848, 16882, 11386, 22298, 21995, 25291, 18675, 12475, 15549, 16421, 16081, 12832, 25318, 20814, 19237, 19047, 13526, 14321, 24270, 25436, 18646, 37614, 13623, 28674, 17059, 19762, 16538, 21445, 20248, 27289, 11007, 26564, 13337, 26110, 11228, 11760, 13956, 22426, 24973, 13035, 19856, 17773, 12376, 18342, 26478, 20491, 20250, 24585, 15958, 12847, 30606, 12865, 27925, 14514, 22670, 13309, 13183, 18478, 13514, 23464, 11298, 23136, 23290, 12540, 17842, 22650, 19298, 25356, 23913, 24998, 18601, 31113, 22466, 16713, 25484, 15277, 25511, 10246, 19344, 13423, 18565, 13506, 25670, 21263, 18243, 10967, 11876, 14519, 19386, 20209, 25855, 20542, 11469, 28035, 19655, 30791, 17742, 10400, 28675, 22552, 22078, 14962, 17831, 15233, 20122, 26183, 22212, 23485, 18396, 20426, 19624, 20544, 12998, 15401, 26576, 15047, 19000, 15883, 21253, 23673, 20675, 21635, 11565, 11263, 26406, 24708, 19694, 17550, 14536, 21069, 10647, 16486, 27166, 23964, 19487, 12544, 17432, 27392, 22842, 22116, 23186, 24159, 18182, 24167, 16972, 17937, 30207, 12411, 21252, 14899, 13243, 18121, 31231, 25164, 10862, 22403, 13496, 17645, 25751, 20338, 30781, 22638, 15650, 12581, 22898, 21775, 13468, 14745, 23166, 20672, 22813, 15810, 11482, 17426, 14472, 13135, 17861, 11369, 21573, 26531, 18849, 14986, 12259, 23361, 24153, 10560, 16534, 13572, 15180, 26688, 17969, 14828, 10775, 21421, 16474, 13315, 13099, 15533, 29128, 17054, 19278, 13015, 24485, 13154, 19505, 24049, 16235, 16822, 15574, 23183, 24288, 18084, 19567, 12794, 22817, 30447, 21412, 28107, 22852, 19904, 11000, 18182, 23615, 23529, 19107, 15231, 30934, 19660, 26559, 18022, 12452, 10579, 23079, 23102, 15149, 25905, 11041, 22988, 17969, 16252, 24042, 13212, 24349, 23768, 11987, 24365, 19703, 11628, 11266, 24407, 22570, 16291, 17647, 14408, 20673, 21952, 11066, 21548, 14255, 12304, 15702, 16155, 11527, 15359, 10483, 23347, 22124, 29273, 10803, 26325, 22977, 15905, 15507, 20382, 31950, 30593, 22443, 10350, 15272, 18895, 23515, 11154, 17045, 19777, 24574, 16352, 22538, 16723, 13983, 12997, 15070, 23779, 27803, 29853, 12262, 23531, 10599, 22554, 16370, 23125, 19031, 16782, 13287, 25755, 14263, 26480, 17831, 13181, 15117, 24317, 19322, 22139, 24192, 24545, 28909, 31255, 14165, 10583, 12685, 14272, 13917, 10616, 12649, 23379, 24862, 24212, 13148, 19542, 24263, 10049, 23499, 29368, 12016, 26354, 19942, 19535, 13807, 17969, 15273, 13278, 17354, 19273, 23525, 11583, 26247, 11076, 17193, 16976, 24584, 21841, 14464, 23154, 15682, 25942, 24647, 11298, 20865, 10166, 10501], [1, 5, 5, 5, 5, 5, 5, 2, 5, 1, 5, 1, 4, 5, 5, 5, 5, 5, 1, 4, 5, 5, 3, 5, 4, 4, 5, 2, 5, 5, 5, 5, 3, 5, 5, 1, 1, 5, 1, 2, 5, 5, 3, 5, 5, 4, 5, 1, 5, 1, 1, 1, 5, 3, 5, 4, 1, 4, 5, 3, 2, 5, 5, 1, 5, 5, 3, 1, 4, 3, 1, 1, 1, 5, 3, 5, 5, 1, 5, 3, 5, 3, 5, 5, 5, 1, 5, 4, 5, 5, 1, 5, 5, 2, 5, 2, 4, 5, 5, 5, 3, 1, 5, 5, 5, 5, 3, 5, 1, 3, 2, 5, 5, 1, 1, 3, 5, 3, 5, 1, 5, 3, 1, 3, 5, 5, 1, 5, 4, 5, 1, 3, 5, 5, 2, 5, 1, 3, 5, 1, 3, 2, 4, 1, 5, 5, 3, 5, 3, 3, 3, 3, 2, 5, 1, 2, 4, 5, 3, 5, 5, 5, 2, 5, 5, 1, 5, 3, 1, 4, 3, 5, 5, 2, 5, 5, 1, 5, 2, 3, 4, 5, 3, 5, 5, 1, 3, 2, 5, 3, 5, 5, 5, 1, 3, 5, 3, 5, 4, 5, 5, 5, 4, 4, 1, 4, 1, 4, 5, 5, 3, 3, 5, 1, 1, 4, 5, 5, 5, 2, 3, 5, 4, 1, 3, 3, 5, 4, 1, 5, 3, 1, 5, 1, 5, 3, 1, 4, 3, 1, 1, 5, 5, 3, 5, 3, 3, 3, 1, 1, 3, 1, 1, 2, 5, 5, 2, 5, 1, 4, 4, 4, 1, 5, 5, 3, 3, 5, 4, 3, 5, 5, 1, 1, 1, 1, 3, 1, 5, 1, 5, 3, 5, 3, 3, 5, 2, 2, 5, 5, 5, 5, 3, 1, 5, 5, 5, 1, 5, 5, 4, 5, 4, 5, 3, 5, 3, 2, 1, 5, 1, 4, 5, 5, 3, 1, 1, 2, 5, 5, 5, 3, 4, 5, 3, 3, 1, 3, 1, 3, 1, 5, 3, 1, 2, 5, 3, 5, 5, 1, 3, 5, 3, 4, 1, 5, 1, 5, 5, 3, 3, 4, 2, 1, 2, 5, 2, 5, 4, 5, 5, 5, 5, 4, 5, 3, 2, 5, 5, 3, 5, 2, 2, 4, 1, 5, 4, 5, 1, 5, 5, 5, 5, 3, 5, 3, 5, 3, 5, 1, 5, 4, 3, 5, 2, 3, 5, 5, 5, 5]]"
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

