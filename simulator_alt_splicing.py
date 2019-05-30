#!/anaconda2/bin/python
import random
import numpy as np
import pandas
import os
import subprocess
from subprocess import Popen, PIPE

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
labelings = [5]#[5,10,20,60]
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
h_s_alt_e = [20]
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
    print spliced
    #TODO: Get the reads from the transcripts and map them to the gene
    reads = [get_reads(x[0]) for x in spliced[0:10]] ## TODO: TESTING remove [0:10] later
    return reads

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
    t = None # splice type
    final_len = 0 # length of transcript after splicing
    
    
    # too short to splice
    if(end_site < intron_1+u_dist):
        #print "too short "+str(end_site)
        t = 4
        final_len = u_dist+end_site
    
    
    # alternatively spliced / not spliced SHORT #NOTE: should this even be here?
    elif((end_site >= intron_1+u_dist) and (end_site < intron+u_dist)): #TODO: need to determine >= or just >?
        ## OLD (possibly wrong)
        ## runs a probability based on psi_se, transcription_rate, and half_life # TODO: for now just depends on Psi
        ## n = np.random.choice([0, 3], p=[1-psi_se,psi_se])
        
        # Updated
        n = np.random.uniform(0,1,1).tolist()[0]>2**(-(end_site-intron_1)/h_alt_e/transc_rate)
        
        if n:
            #print "Short unspliced "+str(end_site)
            t = 0
            final_len = u_dist+end_site
        else:
            #print "Short spliced "+str(end_site-intron_1)
            t = 3
            final_len = u_dist + min(end_site, intron_1+exon_se+d_dist) - intron_1 

    # alternatively spliced / constitutetively spliced # NOTE: should this also have a not spliced option?
    elif(end_site >= intron+u_dist):
        ## OLD (possibly wrong)
        # n = np.random.choice([1, 2], p=[1-psi_se,psi_se])

        # Updated
        n = np.random.uniform(0,1,1).tolist()[0]>2**(-(end_site-intron)/h/transc_rate)
        m = np.random.uniform(0,1,1).tolist()[0]>2**(-(end_site-intron)/h_alt_e/transc_rate)



        if m: # alternatively spliced
            #print "Spliced alt "+str(end_site-(intron_1+intron_2))
            t = 2
            final_len = u_dist + min(end_site, intron_1+exon_se+d_dist) - (intron_1+intron_2)
        if n: # spliced normally
            #print "Spliced normally "+str(end_site-intron)
            t = 1
            final_len = u_dist + min(end_site, intron+exon+d_dist) - intron
        else: #unspliced
            #print "unspliced"
            t = 0
            final_len = u_dist+ end_site 

    else:
        print "\nERROR: end_site = ",end_site," quitting..."
        quit()

    return [int(final_len),t]

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
    # call a subprocess and get results (which were printed to the output of the R script)
    cwd = os.getcwd()
    process = Popen([cwd+"/weibull_dist_calc.R",str(eta_val),str(length)], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    delta_is =  [int(x) for x in stdout.replace('[1] "c(','').replace(')"','').replace("\n","").replace("\\n","").split(", ")]
    starts = get_reads_helper_start(delta_is, insertsize)
    ends = get_reads_helper_end(delta_is, insertsize)
    lengths = [y-x for x, y in zip(starts, ends)]
    # Filter fragments by length and return
    starts = [starts[i] for i,v in enumerate(starts) if ((lengths[i] >=insertsize[0]) and (lengths[i]<=insertsize[1]))]
    #ends = [ends[i] for i,v in enumerate(ends) if ((lengths[i] >=insertsize[0]) and (lengths[i]<=insertsize[1]))]
    lengths = [x for x in lengths if ((x >=insertsize[0]) and (x<=insertsize[1]))]
    return [starts,lengths]

    
# get the start points of the fragments 
def get_reads_helper_start(delta_is, insertsize):
    if len(delta_is) > 1:
        d = range(min(insertsize[0], delta_is)+1)[1:]
        return random.sample(d,1)+ pandas.Series(delta_is[0:len(delta_is)-1]).cumsum().values.tolist()
    else:
        return [random.sample(min(insertsize[0], delta_is),1)]
    
# get the end points of the fragments 
def get_reads_helper_end(delta_is, insertsize):
    if len(delta_is)>1:
        return pandas.Series(delta_is[0:len(delta_is)-1]).cumsum().values.tolist() + [sum(delta_is) - random.sample(range(min(insertsize[0], (sum(delta_is)-delta_is[-1]))+1)[1:],1)[0]]
    else:
        return delta_is
















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

df = pandas.DataFrame(ls,index = None, columns = None)
df = df.transpose()
df.to_csv("output.csv", sep=',')
