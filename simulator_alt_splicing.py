#!/anaconda2/bin/python
import random
import numpy as np
import pandas

## attributes TODO: make so read from command line
d_dists = [3] # d distances?
expr_lvls = [4]

## parameters
# constants
exon = 0.2
u_dist = 0.5
n_millions = 100 # total number of millions of transcripts to consider
transc_rate = 1.5 # rate of transcription
# variables (to simulate over)
labelings = [5,10,20,60]
#h_s = list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
h_s = [0.2]
introns =[40.0] #NOTE: need to make larger because alt spliced introns are larger (look up)
#introns =list(np.arange(0.04,0.09,0.02)+list(np.arange(0.1,1,0.1))+list(np.arange(1,50,2))

#d_dists = list(np.arange(0.5, 5, 0.5))
#expr_lvls = list(np.arange(1,100, 4))

# for alt splicing TODO: just for testing need to change to simulate over
#                   NOTE: also in future need to pick actual values because now these are determined by intron length and anything 
#                       that isn't divisible by 4 will give a non-integral response which won't work in the simulations
#h_s_alt_e = list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
h_s_alt_e = [0.2]
exon_se = introns[0]/4
intron_1 = (introns[0]-exon_se)/2
intron_2 = introns[0]-(exon_se+intron_1)
psi_se=[0.5]#list(np.arange(0.0,1,0.1)) #Psi of SE (skipped exon) ## TODO: need a function that calculates this

#
#       MODEL:
#
#                  ---------------------------intron---------------
#
#     [  u_dist   ]----intron_1----[/// exon_se ///]---intron_2----[  exon  (d_dist?)  AAAAAA]
#


## Simulation Function
## need to loop through some variables/params to run this
def simulate(intron, exon, u_dist, d_dist, labeling, h, expr_lvl, n_millions, transc_rate, intron_1, intron_2, exon_se, h_alt_e, psi_se):
    # simulate transcription of region - generate a list of end sites for a transcript that are distributed 
    # depending on where simulated transcription started
    #
    # Outputs: (in list form)
    #   [[1]] data frame of all reads with the start position, the matching pattern, and a 
    #         string with all the associated parameters (columns: start, match, name)
    #   [[2]] the number of spliced junction reads
    #   [[3]] the number of unspliced junction reads
    #   [[4]] whether or not each read comes from a spliced or unspliced transcript
    #
    intron = 1000*intron
    exon   = 1000*exon
    u_dist = 1000*u_dist
    d_dist = 1000*d_dist
    transc_rate = 1000*transc_rate
    intron_1 = 1000*intron_1
    intron_2 = 1000*intron_2
    exon_se = 1000*exon_se

    # Generate expression_level*n_millions transcripts uniformly from the labeled region 
    # end_sites = round(seq(from = 1, to = intron + exon + D_dist + labeling*transcription_rate,
    #                     length.out = expression_level*n_millions))
    end_sites =  random.sample(range(1,int(intron + exon + d_dist + labeling*transc_rate)), expr_lvl*n_millions)
    end_sites =[u_dist+x for x in end_sites]
    
    # Determine if the transcripts are spliced or not
    # And determine resulting lengths of transcripts that were spliced
    spliced = [splice(x, intron, intron_1, exon_se, u_dist, h, h_alt_e, transc_rate, psi_se) for x in end_sites]

     #TODO: Get the reads from the transcripts and map them to the gene


    return spliced

# Determine if (and how) a given transcript is spliced and their resulting lengths (possibly make a separate function)
#  ## TODO: update so uses exponential function (see slides from ap)
#
# Spliced only if:
#   1) Long enough:
#       a) For alt:
#           end_site > intron_1
#       b) For constitutive:
#           end_site > intron
#
#   2) Probability of splicing 
#       Depends on:
#           a) psi_se
#           b) transcription_rate
#           c) half_life (of the exon)
#       
#   0 --> unspliced
#   1 --> spliced
#   2 --> alternative spliced (retain exon_se)
#   3 --> alternative spliced short (retain exon_se) but exon not yet transcribed #NOTE: should I include this one?
#   4 --> unspliced (not long enough to splice)
#
def splice(end_site, intron, intron_1, exon_se, u_dist, h, h_alt_e, transc_rate, psi_se): 
    #print "****** "+str(end_site)+" *********"
    # too short
    if(end_site < intron_1+u_dist):
        #print "too short "+str(end_site)
        return [end_site,4]
    # alternatively spliced / not spliced SHORT #NOTE: should this even be here?
    elif((end_site >= intron_1+u_dist) and (end_site < intron+u_dist)): #TODO: need to determine >= or just >?
        # runs a probability based on psi_se, transcription_rate, and half_life # TODO: for now just depends on Psi
        n = np.random.choice([0, 3], p=[1-psi_se,psi_se])
        if n == 0:
            #print "Short unsplied "+str(end_site)
            return [end_site,0]
        else:
            #print "Short spliced "+str(end_site-intron_1)
            #return u_dist+exon_se
            return [end_site-intron_1,3] # is this correct since haven't fully transcribed exon_se yet?
    # alternatively spliced / constitutetively spliced # NOTE: should this also have a not spliced option?
    elif(end_site >= intron+u_dist):
        n = np.random.choice([1, 2], p=[1-psi_se,psi_se])
        if n == 1: # spliced normally
            #print "Spliced normally "+str(end_site-intron)
            #return u_dist+exon
            return [end_site - (intron), 1]
        else: # alternatively spliced
            #print "Spliced alt "+str(end_site-(intron_1+intron_2))
            #return u_dist + exon_se+exon
            return [end_site - (intron_1+intron_2), 2]
    else:
        print "\nERROR: end_site = ",end_site," quitting..."
        quit()

# Fragment and Read Function
#
# Takes a transcript and selects a fragment from it, then performs size selection, and returns the starting position of the
# resulting reads relative to the length of the transcript
#
# Inputs:
#   lengths - the length of a transcript
#
# Outputs:
#   fragments - a data frame with the start positions of the filtered fragments and the index
#               of the transcript that it came from (columns: transcript, start)
def get_reads(length):
    eta_val = 200 # the eta value input to the Weibull distribution
    insertsize = [200, 300] # a list of length two with the lower and upper bound for fragments

    # sample lengths from a weibull distribution for the transcript and transform them to the length
    # of the transcript
    #
    # NOTE: may need to call subprocess to do weibull stats
    cwd = os.getcwd()
    subprocess.call([cwd+"/weibull_dist_calc.R",eta_val, length ])

  deltas = log10(lengths)
  ns_minus_1 = pmax(round(lengths/eta_val/gamma(1/deltas + 1)) - 1, 0)
  xis = lapply(ns_minus_1, function(n) {diff(sort(c(runif(n), 0, 1)))})
  xis_transformed = mapply(function(x, d) {x^(1/d)}, xis, deltas, SIMPLIFY = F)
  delta_is = mapply(function(len, x_t) {round(len*x_t/sum(x_t))}, lengths, xis_transformed, SIMPLIFY = F)























## START 
#d_dist = d_dists[0]
#expr_lvl = expr_lvls[0]
#labeling = labelings[0]
#h = h_s[0]
#intron = introns[0]
#h_alt_e = h_s_alt_e[0]
intron_1 = intron_1
intron_2 = intron_2

ls = []
for labeling in labelings:
    for d_dist in d_dists:
        for expr_lvl in expr_lvls:
            for h in h_s:
                for h_alt_e in h_s_alt_e:
                    for intron in introns:
                        for psi in psi_se: # NOTE: should I iterate through these or just calculate?
                            e = simulate(intron, exon, u_dist, d_dist, labeling, h, expr_lvl, n_millions, transc_rate, intron_1, intron_2, exon_se, h_alt_e, psi)
                            
                            ls.append(e)

print e
df = pandas.DataFrame(ls,index = None, columns = None)
df = df.transpose()
df.to_csv("output.csv", sep=',')
