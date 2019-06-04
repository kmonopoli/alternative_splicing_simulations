#!/anaconda2/bin/python
import random
import numpy as np
import pandas
import os
import subprocess
from subprocess import Popen, PIPE




## attributes TODO: make so read from command line
d_dists = [3] # list(np.arange(0.5, 5, 0.5))
expr_lvls = [4] #list(np.arange(1,100, 4))

## parameters
# constants
exon = 0.3
u_dist = 0.5
n_millions = 100 # total number of millions of transcripts to consider
transc_rate = 1.5 # rate of transcription


# variables (to simulate over)
labelings = [5]#[5,10,20,60]
h_s = [0.2]#[100]#[0.2] # list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
introns =[0.04]#[0.04]#[40.0] #list(np.arange(0.04,0.09,0.02)+list(np.arange(0.1,1,0.1))+list(np.arange(1,50,2)) #NOTE: need to make larger because alt spliced introns are larger (look up)

# for alt splicing 
#   TODO: just for testing need to change to simulate over
#   NOTE: also in future need to pick actual values because now these are determined by intron length and anything 
#   that isn't divisible by 4 will give a non-integral response which won't work in the simulations
h_s_alt_e = [40] # list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
exon_se = introns[0]/4
intron_1 = (introns[0]-exon_se)/2
intron_2 = introns[0]-(exon_se+intron_1)

psi_se_s=[0]#[0.5]#list(np.arange(0.0,1,0.1)) #Psi of SE (skipped exon) ## TODO: need a function that calculates this






## Simulation Function
## need to loop through some variables/params to run this
#
# Outputs: (in list form)
#   [[1]] TODO: data frame of all reads with the start position, the matching pattern, and a string with all the associated parameters (columns: start, match, name)
#   [[2]] the number of spliced junction reads
#   [[3]] the number of unspliced junction reads
#   [[4]] whether or not each read comes from a spliced or unspliced transcript
def simulate(intron, exon, u_dist, d_dist, labeling, h, expr_lvl, n_millions, transc_rate, intron_1, intron_2, exon_se, h_alt_e, psi_se):
    intron = 1000*intron
    exon   = 1000*exon
    u_dist = 1000*u_dist
    d_dist = 1000*d_dist
    transc_rate = 1000*transc_rate
    intron_1 = 1000*intron_1
    intron_2 = 1000*intron_2
    exon_se = 1000*exon_se
    rl = 50
    # Generate expression_level*n_millions transcripts uniformly from the labeled region 
    # end_sites = round(seq(from = 1, to = intron + exon + D_dist + labeling*transcription_rate,
    #                     length.out = expression_level*n_millions))
    end_sites =  random.sample(range(1,int(intron + exon + d_dist + labeling*transc_rate)), expr_lvl*n_millions)
    end_sites =[u_dist+x for x in end_sites]

   

    # Determine if the transcripts are spliced or not
    # And determine resulting lengths of transcripts that were spliced
    spliced = [splice(x, intron, intron_1, intron_2, exon_se, u_dist, d_dist, h, h_alt_e, transc_rate, psi_se) for x in end_sites]
    # Get the reads from the transcripts and map them to the gene
    start_reads = get_reads(spliced)
    start_pos = pandas.DataFrame({"transcript":start_reads[0],"start":start_reads[1]})
    print start_pos
    
    # Get splice types (numeric)
    start_pos.insert(0,"splice_type",[spliced[x-1][1] for x in start_pos['transcript']])
    
    
    # Get read start positions
    ## TODO: need to update to allow for alt splicing
    reads = pandas.DataFrame({"start":start_pos['start'] - u_dist + intron *(start_pos['splice_type']==1)*(start_pos['start'] > u_dist)})
    
    # Determine if read is junction read
    ## TODO: need to update for alt splicing
    reads.insert(1, "junction", ((reads['start'] > (-1*rl +10)) * (reads['start'] <= -9) * (start_pos['splice_type'] == 1)), True) 


    # Data to return
    spliced_num = sum(reads['junction'])
    unspliced_num = sum((start_pos['splice_type']==0) & (start_pos['start'] > (intron + u_dist-rl+10)) & (start_pos['start'] <= (intron + u_dist-9)))
    intron_num = sum((start_pos['splice_type']==0) & (start_pos['start'] >=1 ) & (start_pos['start'] < 50) )

    return [reads,intron_num,spliced_num,unspliced_num,start_pos['splice_type']]

# Determine if (and how) a given transcript is spliced and their resulting lengths (possibly make a separate function)
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
#    MODEL:
#                  ---------------------------intron---------------
#
#     [  u_dist   ]----intron_1----[/// exon_se ///]---intron_2----[  exon  (d_dist?)  AAAAAA]
#
#
#
def splice(end_site, intron, intron_1, intron_2, exon_se, u_dist, d_dist, h, h_alt_e, transc_rate, psi_se): 
    
    t = -1 # splice type
    final_len = 0 # length of transcript after splicing
    
    if(end_site < intron_1+u_dist): # too short to splice
        
        t = 4
        final_len = u_dist+end_site
    
    # short alternatively spliced / not spliced  
    elif((end_site >= intron_1+u_dist) and (end_site < intron+u_dist)):
        
        n = np.random.uniform(0,1,1).tolist()[0]>2**(-1*(end_site-intron_1)/h_alt_e/transc_rate)
        
        if n: # short alt. spliced
            t = 3
            final_len = u_dist + min(end_site, intron_1+exon_se+d_dist) - intron_1
            
        else: # short unspliced
            t = 4
            final_len = u_dist+end_site
  
    # alternatively spliced / constitutetively spliced / unspliced
    elif(end_site >= intron+u_dist):
        # n = np.random.choice([1, 2], p=[1-psi_se,psi_se])

        n = np.random.uniform(0,1,1).tolist()[0]>2**(-1*(end_site-intron)/h/transc_rate)
        m = np.random.uniform(0,1,1).tolist()[0]>2**(-1*(end_site-intron)/h_alt_e/transc_rate)

        if m: # alternatively spliced
            t = 2
            final_len = u_dist + min(end_site, intron_1+intron_2+exon_se+d_dist) - (intron_1+intron_2)
        elif n: # spliced normally
            t = 1
            final_len = u_dist + min(end_site, intron+exon+d_dist) - intron
        else: # unspliced
            t = 0
            final_len = u_dist+ end_site 
        
    else:
        print "\n ERROR: end_site = ",end_site," quitting..."
        quit()
#    print end_site, " <-- *"
#    print "intron_1 -> ",intron_1
#    print "exon_se -> ",exon_se
#    print "intron_2 -> ",intron_2
#    print "d_dist -> ",d_dist
#    print "u_dist -> ",u_dist
##    print "sum -> ",intron_1+exon_se+intron_2
##    print "intron -> ",intron
#    print "splice type --> ",t
#    print "final length --> ", final_len
#    print "\n"
    return [int(final_len),t]





# Fragment and Read Function
#
# Takes a transcript and selects a fragment from it, then performs size selection, and returns the starting position of the
# resulting reads relative to the length of the transcript
#
# Inputs:
#   data = [lengths, splice_type]
#       lengths - the length of a transcript
#       splice type
#
# Outputs:
#   fragments - an array of arrays ith the start positions of the filtered fragments and the index
#               of the transcript that it came from (columns: transcript, start)
def get_reads(data):
    lengths = [x[0] for x in data]
    if (sum(1 for l in lengths if l < 0)) >0 :
        print "WARNING: some lengths passed into get_reads function are negative, which will not work when calculating weibull distribution"
    
    splice_type = [x[1] for x in data]
    data = [lengths,splice_type]

    cwd = os.getcwd()
    process = Popen([cwd+"/weibull_dist_calc.R",str(data)], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    stdout = [map(int, q) for q in [x.replace(")","").replace("\n","").replace(" ","").replace('"','').split(",") for x in stdout.split("c(")][1:]]
    if( stdout == []):
        print "ERROR in get_reads function that calls weibull_dist_calc.R -  no output received"
        #quit()
    return stdout


    
    
    

# maps reads to sequence, returns start_sites based on the type of splicing that occured
def map_data(df, u_dist, intron, intron_1, intron_2):
    start_sites = []
    i = 0
    while i<df.shape[0]:
        s = (df['start_pos'].iloc[i]-int(u_dist))
        l = int(df['start_pos'].iloc[i]>u_dist) # 0 if too short, 1 if past u_dist
        n = 0 
        if  (df['splice_type'].iloc[i] == 0 or df['splice_type'].iloc[i] == 4): # 0 --> unspliced | 4 --> unspliced (not long enough to splice)
            n = 0
        elif(df['splice_type'].iloc[i] == 1): # 1 --> spliced
            n=-1*int(intron)
        elif(df['splice_type'].iloc[i] == 2): # 2 --> alternative spliced (retain exon_se)
            n=-1*(int(intron_1+intron_2))
        elif(df['splice_type'].iloc[i] == 3): # 3 --> alternative spliced short (retain exon_se) but exon not yet transcribed #NOTE: should I include this one?
            n=-1*(int(intron_1))
        else:
            print "ERROR splice type not supported: ",(df['splice_type'].iloc[i])
            quit()

        start_sites.append(s+l*n)
        i+=1

    # add reads to dataframe
    df.insert(0, "start_site", start_sites, True)





















## START 


ls = []
# TODO: need to run these for loops to go through all possibilities
#for labeling in labelings:
#    for d_dist in d_dists:
#        for expr_lvl in expr_lvls:
#            for h in h_s:
#                for h_alt_e in h_s_alt_e:
#                    for intron in introns:
#                        for psi in psi_se: # NOTE: should I iterate through these or just calculate?
#                            #e = simulate(intron, exon, u_dist, d_dist, labeling, h, expr_lvl, n_millions, transc_rate, intron_1, intron_2, exon_se, h_alt_e, psi)
##################################
## FOR TESTING: TODO: iterate through in real simulation
d_dist = d_dists[0]
expr_lvl = expr_lvls[0]
labeling = labelings[0]
h = h_s[0]
intron = introns[0]
h_alt_e = h_s_alt_e[0]
psi_se = psi_se_s[0]
##################################


intron = 1000*intron
exon   = 1000*exon
u_dist = 1000*u_dist
d_dist = 1000*d_dist
transc_rate = 1000*transc_rate
intron_1 = 1000*intron_1
intron_2 = 1000*intron_2
exon_se = 1000*exon_se
rl = 50
# Generate expression_level*n_millions transcripts uniformly from the labeled region 
# end_sites = round(seq(from = 1, to = intron + exon + D_dist + labeling*transcription_rate,
#                     length.out = expression_level*n_millions))
end_sites =  random.sample(range(1,int(intron + exon + d_dist + labeling*transc_rate)), expr_lvl*n_millions)
end_sites =[u_dist+x for x in end_sites]

   
# Determine if the transcripts are spliced or not
# And determine resulting lengths of transcripts that were spliced
spliced = [splice(x, intron, intron_1, intron_2, exon_se, u_dist, d_dist, h, h_alt_e, transc_rate, psi_se) for x in end_sites]

# Get the reads from the transcripts and map them to the gene
start_reads = get_reads(spliced)

start_pos = pandas.DataFrame({"transcript":start_reads[0],"start":start_reads[1]})

# Get splice types (numeric)
start_pos.insert(0,"splice_type",[spliced[x-1][1] for x in start_pos['transcript']])

# Get read start positions
## TODO: need to update to allow for alt splicing
reads = pandas.DataFrame({"start":start_pos['start'] - u_dist + intron *(start_pos['splice_type']==1)*(start_pos['start'] > u_dist)})

# Determine if read is junction read
## TODO: need to update for alt splicing
reads.insert(1, "junction", ((reads['start'] > (-1*rl +10)) * (reads['start'] <= -9) * (start_pos['splice_type'] == 1)), True) 


# Data to return
spliced_num = sum(reads['junction'])
unspliced_num = sum((start_pos['splice_type']==0) & (start_pos['start'] > (intron + u_dist-rl+10)) & (start_pos['start'] <= (intron + u_dist-9)))
intron_num = sum((start_pos['splice_type']==0) & (start_pos['start'] >=1 ) & (start_pos['start'] < 50) )




#df = pandas.DataFrame(e,index = None, columns = None)
#df0 = e[0]
#print e

e = start_pos
e.to_csv (os.getcwd()+'/export_dataframe.csv', index = None, header=True, sep=',')
#ls.append(e)

                            
                         
#df = pandas.DataFrame(ls,index = None, columns = None)
#df = df.transpose()
#df.to_csv("output.csv", sep=',')
