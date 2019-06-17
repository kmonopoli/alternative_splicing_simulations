#!/anaconda2/bin/python
import random
import numpy as np
import pandas
import os
import subprocess
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt





## attributes TODO: make so read from command line
d_dists = [0.5]#list(np.arange(0.5, 5, 0.5))
expr_lvls = [100]#[4] #list(np.arange(1,100, 4))

## parameters
# constants
exon = 0.3
u_dist = 0.5
n_millions = 100 # total number of millions of transcripts to consider
transc_rate = 1.5 # rate of transcription


# variables (to simulate over)
labelings = [5]#[5,10,20,60]
hs_intron_3 = [0.2]#[100]#[0.2] # list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
introns_3 =[2.0]#[40.0] #list(np.arange(0.04,0.09,0.02)+list(np.arange(0.1,1,0.1))+list(np.arange(1,50,2)) #NOTE: need to make larger because alt spliced introns are larger (look up)

# for alt splicing 
#   TODO: just for testing need to change to simulate over
#   NOTE: also in future need to pick actual values because now these are determined by intron length and anything 
#   that isn't divisible by 4 will give a non-integral response which won't work in the simulations
hs_intron_1 = [20.0]#[40.0] # list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
hs_intron_2 = [20.0]#[40.0] # list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))

exon_se = introns_3[0]/4
introns_1 = [(introns_3[0]-exon_se)/2]
introns_2 = [introns_3[0]-(exon_se+introns_1[0])]

psi_se_s=[1]#[0.5]#list(np.arange(0.0,1,0.1)) #Psi of SE (skipped exon) ## TODO: need a function that calculates this






## Simulation Function
## need to loop through some variables/params to run this
#
# Outputs: (in list form)
#   [[1]] TODO: data frame of all reads with the start position, the matching pattern, and a string with all the associated parameters (columns: start, match, name)
#   [[2]] the number of spliced junction reads
#   [[3]] the number of unspliced junction reads
#   [[4]] whether or not each read comes from a spliced or unspliced transcript
def simulate(intron_3, exon, u_dist, d_dist, labeling, h_intron_3, expr_lvl, n_millions, transc_rate, intron_1, intron_2, exon_se, h_intron_1, h_intron_2, psi_se):

    
    intron_3 = 1000*intron_3
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
#    end_sites =  random.sample(range(1,int(intron_3 + exon + d_dist + labeling*transc_rate)), expr_lvl*n_millions)
    end_sites =  np.random.choice(range(1,int(intron_3 + exon + d_dist + labeling*transc_rate)), expr_lvl*n_millions)

    end_sites =[u_dist+x for x in end_sites]

#    print len(end_sites)
    # Determine if the transcripts are spliced or not
    # And determine resulting lengths of transcripts that were spliced
    spliced = [splice(x, intron_3, intron_1, intron_2, exon_se, u_dist, d_dist, h_intron_3, h_intron_1, h_intron_2, transc_rate, psi_se) for x in end_sites]
    # transpose data and put in dataframe
    spliced_df = pandas.DataFrame({"lengths":[x[0] for x in spliced],"splice_type":[x[1] for x in spliced]})
    
    # Get the reads from the transcripts and map them to the gene
    start_reads = get_reads(spliced) # TODO : this should use rl?
    
    start_pos = pandas.DataFrame({"transcript":start_reads[0],"start":start_reads[1]})
    
    # Get splice types (numeric)
    start_pos.insert(0,"splice_type",[spliced[x-1][1] for x in start_pos['transcript']])
    
        
    
    # Get read start positions
    start_pos["read_start"] = np.nan
    #    1 --> intron_1 spliced out
    start_pos.read_start.loc[(start_pos['splice_type'] == 1)] = (start_pos['start'] - u_dist) + intron_1 * (start_pos['start'] > u_dist) # TODO: not sure how to do this
    #    2 --> intron_2 spliced out
    start_pos.read_start.loc[(start_pos['splice_type'] == 2)] = (start_pos['start'] - u_dist) + intron_2 * (start_pos['start'] > u_dist) # TODO: not sure how to do this
    #    3 --> intron_1 and intron_2 spliced out (exon_se included)
    start_pos.read_start.loc[(start_pos['splice_type'] == 3)] = (start_pos['start'] - u_dist) + intron_1 * (start_pos['start'] > u_dist) + intron_2 * (start_pos['start'] > u_dist+intron_1)# TODO: not sure how to do this
    #    4 --> intron_3 spliced out
    start_pos.read_start.loc[(start_pos['splice_type'] == 4)] = (start_pos['start'] - u_dist) + intron_3 * (start_pos['start'] > u_dist)
    #    5 --> unspliced
    start_pos.read_start.loc[(start_pos['splice_type'] == 5)] = (start_pos['start'] - u_dist) + intron_3 * (start_pos['start'] > u_dist)

    # Determine if read is exon-exon (ee) junction read (only for splice types 1,2,3)
    start_pos["junction"] = "not junction read"
    junct_read_overlap = 10
    start_pos.junction.loc[(start_pos['splice_type'] == 1) & # intron_1 excluded
                           (start_pos['read_start']  > ((-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <= (-1*junct_read_overlap))] = "exon_us - exon_se"
    start_pos.junction.loc[(start_pos['splice_type'] == 2) & # intron_2 excluded
                           (start_pos['read_start']  > ((intron_1 + exon_se) +(-1*rl) + junct_read_overlap)) &
                           (start_pos['read_start']  <=(( intron_1 + exon_se))+(-1*junct_read_overlap))] = "exon_se - exon"
    start_pos.junction.loc[(start_pos['splice_type'] == 3) &  # intron_1 and intron_2 excluded (exon_se included)
                           ((start_pos['read_start']  > ((-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <=(-1*junct_read_overlap)))] = "exon_us - exon_se"
    start_pos.junction.loc[(start_pos['splice_type'] == 3) &  # intron_1 and intron_2 excluded (exon_se included)
                           ((start_pos['read_start']  > ((intron_1 + exon_se) +(-1*rl) + junct_read_overlap)) &
                           (start_pos['read_start']  <=((intron_1 + exon_se))+(-1*junct_read_overlap)))] = "exon_se - exon"
    start_pos.junction.loc[(start_pos['splice_type'] == 4) & # intron_3 excluded
                           (start_pos['read_start']  > ((-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <=(-1*junct_read_overlap))] = "exon_us - exon"
    
    # TODO: Determine if read is intron-exon (ie) junction read 
    start_pos.junction.loc[(start_pos['splice_type'] == 1) & # intron_1 excluded
                           (start_pos['read_start']  > ((intron_1 + exon_se)+(-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <= ((intron_1 + exon_se)+(-1*junct_read_overlap)))] = "exon_se - intron_2"
    start_pos.junction.loc[(start_pos['splice_type'] == 1) & # intron_1 excluded
                           (start_pos['read_start']  > ((intron_1 + exon_se + intron_2)+(-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <= ((intron_1 + exon_se + intron_2)+(-1*junct_read_overlap)))] = "intron_2 - exon"
    start_pos.junction.loc[(start_pos['splice_type'] == 2) & # intron_2 excluded
                           (start_pos['read_start']  > ((-1*rl) + junct_read_overlap)) &
                           (start_pos['read_start']  <=(-1*junct_read_overlap))] = "exon_us - intron_1"
    start_pos.junction.loc[(start_pos['splice_type'] == 2) & # intron_2 excluded
                           (start_pos['read_start']  > (intron_1 +(-1*rl) + junct_read_overlap)) &
                           (start_pos['read_start']  <=(intron_1 +(-1*junct_read_overlap)))] = "intron_1 - exon_se"
                           
    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
                           (start_pos['read_start']  > ((-1*rl) + junct_read_overlap)) &
                           (start_pos['read_start']  <=(-1*junct_read_overlap))] = "exon_us - intron_1"
    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
                           (start_pos['read_start']  > (intron_1 +(-1*rl) + junct_read_overlap)) &
                           (start_pos['read_start']  <=(intron_1 +(-1*junct_read_overlap)))] = "intron_1 - exon_se"
    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
                           (start_pos['read_start']  > ((intron_1 + exon_se)+(-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <= ((intron_1 + exon_se)+(-1*junct_read_overlap)))] = "exon_se - intron_2"
    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
                           (start_pos['read_start']  > ((intron_1 + exon_se + intron_2)+(-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <= ((intron_1 + exon_se + intron_2)+(-1*junct_read_overlap)))] = "intron_2 - exon"
    return start_pos , spliced_df #e,f

# Determine if (and how) a given transcript is spliced and their resulting lengths 
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
#   OLD: 
#   0 --> unspliced
#   1 --> exclude intron_3
#   2 --> exclude intron_1 and intron_2
#   3 --> alternative spliced short (retain exon_se) but exon not yet transcribed #NOTE: should I include this one?
#   4 --> unspliced (not long enough to splice)
#
#   
#   NEW:
#    1 --> intron_1 excluded
#    2 --> intron_2 excluded
#    3 --> intron_1 and intron_2 excluded (exon_se included)
#    4 --> intron_3 excluded
#    5 --> unspliced
#
#    MODEL:
#                  ---------------------------intron_3---------------
#
#     [  u_dist   ]----intron_1----[/// exon_se ///]---intron_2----[  exon ][ (d_dist?)  ][AAAAAA (L*r)]
#
#
#
def splice(end_site, intron_3, intron_1, intron_2, exon_se, u_dist, d_dist, h_intron_3, h_intron_1, h_intron_2, transc_rate, psi_se): 
    # TODO: include psi here (see notes) ## n = np.random.choice([1, 2], p=[1-psi_se,psi_se])
    t = -1 # splice type
    final_len = 0 # length of transcript after splicing
    spl_path = np.random.choice([0, 1], p=[1-psi_se,psi_se]) # splice path - whether it proceed down path of alt splicing (include exon_se), or normal splicing (exclude exon_se)
    # spl_path will be 1 if alt splicing (include exon_se), 0 if normal splicing (exclude exon_se)
    
    # Alternative Splicing (include exon_se)
    if(spl_path == 1): # TODO be sure 1 and 0 are not flipped for Psi 
        if(end_site < intron_1+u_dist): # too short to splice
            t = 5 # unspliced
            final_len = u_dist+end_site
    
        # short - exclude intron_1 / unspliced
        elif((end_site >= intron_1+u_dist) and (end_site < intron_3+u_dist)):
            n = np.random.uniform(0,1,1).tolist()[0]>2**(-1*(end_site-intron_1)/h_intron_1/transc_rate)
            if n: # intron_1 spliced out
                t = 1
                #OLD:final_len = u_dist + min(end_site, intron_1+exon_se+d_dist) - intron_1
                final_len = u_dist+(end_site-intron_1)
            else: # unspliced
                t = 5
                final_len = u_dist+end_site
        # exclude intron_1 and intron_2 / exlude intron_1 / exclude intron_2 / unspliced
        elif(end_site >= intron_3+u_dist):
            # intron_1
            m = np.random.uniform(0,1,1).tolist()[0]>2**(-1*(end_site-intron_3)/h_intron_1/transc_rate)
            # intron_2
            n = np.random.uniform(0,1,1).tolist()[0]>2**(-1*(end_site-intron_3)/h_intron_2/transc_rate)

            if (n and m): # exclude intron_1 and intron_2
                t = 3
                #OLD: final_len = u_dist + min(end_site, intron_1+intron_2+exon_se+d_dist) - (intron_1+intron_2)
                final_len = u_dist + end_site-(intron_1+intron_2)
            elif m: # exclude intron_1 
                t = 1
                final_len = u_dist + end_site-(intron_1)
            elif n: # exclude intron_2
                t = 2
                final_len = u_dist + end_site-(intron_1)
            else: # unspliced
                t = 5
                final_len = u_dist+ end_site 
        else:
            print "\n ERROR: end_site = ",end_site," quitting..."
            quit()
            
    # Normal Splicing (exclude exon_se)
    else: 
        if(end_site < intron_1+u_dist): # too short to splice
            t = 5 # unspliced
            final_len = u_dist+end_site
    
        # short - unspliced
        elif((end_site >= intron_1+u_dist) and (end_site < intron_3+u_dist)):
            t = 5
            final_len = u_dist+end_site
        # exclude intron_3 / unspliced
        elif(end_site >= intron_3+u_dist):
            # intron_3
            n = np.random.uniform(0,1,1).tolist()[0]>2**(-1*(end_site-intron_3)/h_intron_3/transc_rate)
            if n: # exclude intron_3
                t = 4
#                final_len = u_dist + min(end_site, intron_3+exon+d_dist) - intron_3
                final_len = u_dist + end_site - intron_3
            else: # unspliced
                t = 5
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
    data = [[x[0] for x in data],[x[1] for x in data]]

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

################################## ################################## ##################################
## FOR TESTING: TODO: iterate through in real simulation
#d_dist = d_dists[0]
expr_lvl = expr_lvls[0]
labeling = labelings[0]
h_intron_3 = hs_intron_3[0]
intron_3 = introns_3[0]
h_intron_1 = hs_intron_1[0]
h_intron_2 = hs_intron_2[0]
psi_se = psi_se_s[0]
lebeling = labelings[0]
intron_1 = introns_1[0]
intron_2 = introns_2[0]

# TODO: need to run these for loops to go through all possibilities
file_num=1
#for labeling in labelings:
for d_dist in d_dists:
#        for expr_lvl in expr_lvls:
#            for intron_1 in introns_1: 
#                for intron_2 in introns_2:
#                    for intron_3 in introns_3:
#                        for h_intron_1 in hs_intron_1:
#                            for h_intron_2 in hs_intron_2:
#                                for h_intron_3 in hs_intron_3:
#                                    for psi_se in psi_se_s: 
    e,f = simulate(intron_3, exon, u_dist, d_dist, labeling, h_intron_3, expr_lvl, n_millions, transc_rate, intron_1, intron_2, exon_se, h_intron_1, h_intron_2, psi_se)
    #
    # For plotting splice diagram line data
    intron_3 = 1000*intron_3
    exon   = 1000*exon
    u_dist = 1000*u_dist
    d_dist = 1000*d_dist
    transc_rate = 1000*transc_rate
    intron_1 = 1000*intron_1
    intron_2 = 1000*intron_2
    exon_se = 1000*exon_se
    
    
    
    ## add data for test plots
    e.insert(0,"h_intron_1",len(e)*[h_intron_1])
    e.insert(0,"h_intron_2",len(e)*[h_intron_2])
    e.insert(0,"h_intron_3",len(e)*[h_intron_3]) # normal splicing
    e.insert(0,"transc_rate",len(e)*[transc_rate])
    e.insert(0,"intron_3",len(e)*[intron_3])
    e.insert(0,"intron_1",len(e)*[intron_1])
    e.insert(0,"intron_2",len(e)*[intron_2])
    e.insert(0,"u_dist",len(e)*[u_dist])
    e.insert(0,"exon",len(e)*[exon])
    e.insert(0,"exon_se",len(e)*[exon_se])
    e.insert(0,"psi_se",len(e)*[psi_se])
    ## Relabel splice types to readable words
    e.insert(0,"splice_type_read",e['splice_type'])
    e['splice_type_read'] = np.where(e['splice_type_read'] == 1, "intron_1 excluded", e['splice_type_read'])
    e['splice_type_read'] = np.where(e['splice_type_read'] == '2', "intron_2 excluded", e['splice_type_read'])
    e['splice_type_read'] = np.where(e['splice_type_read'] == '3', "intron_1 & intron_2 excluded", e['splice_type_read'])
    e['splice_type_read'] = np.where(e['splice_type_read'] == '4', "intron_3 excluded", e['splice_type_read'])
    e['splice_type_read'] = np.where(e['splice_type_read'] == '5', "unspliced", e['splice_type_read'])
    
    
    #Counts
    mn,mx =(round(e['read_start'].min()/100,0)*100, round( e['read_start'].max()/100,0)*100+100)
    
    e.loc[e['splice_type_read'] == 'intron_1 excluded','read_start'].value_counts()
    
    
    ct_df = pandas.DataFrame({"ct_intron_1_excluded":e.loc[e['splice_type_read'] == 'intron_1 excluded','read_start'].value_counts(),
                                "ct_intron_2_excluded":e.loc[e['splice_type_read'] == 'intron_2 excluded','read_start'].value_counts(),
                                "ct_intron_1_intron_2_excluded":e.loc[e['splice_type_read'] == 'intron_1 & intron_2 excluded','read_start'].value_counts(),
                                "ct_intron_3_excluded":e.loc[e['splice_type_read'] == 'intron_3 excluded','read_start'].value_counts(),
                                "ct_unspliced":e.loc[e['splice_type_read'] == 'unspliced','read_start'].value_counts()
    #                            "bin_start":(division1[:-1]),
                                })
    ct_df.insert(0,"index",ct_df.index)
    
    
    
    f.insert(0,"splice_type_read",f['splice_type'])
    f['splice_type_read'] = np.where(f['splice_type_read'] ==  1, "intron_1 excluded", f['splice_type_read'])
    f['splice_type_read'] = np.where(f['splice_type_read'] == '2', "intron_2 excluded", f['splice_type_read'])
    f['splice_type_read'] = np.where(f['splice_type_read'] == '3', "intron_1 & intron_2 excluded", f['splice_type_read'])
    f['splice_type_read'] = np.where(f['splice_type_read'] == '4', "intron_3 excluded", f['splice_type_read'])
    f['splice_type_read'] = np.where(f['splice_type_read'] == '5', "unspliced", f['splice_type_read'])
    
    ## get starts and ends for each
    f.insert(0,"end_3",None)
    f.insert(0,"start_3",None)
    f.insert(0,"end_2",None)
    f.insert(0,"start_2",None)
    f.insert(0,"end_1",None)
    f.insert(0,"start_1",0) # always starts at 0
    
    # start_1 
    #    (will be same for all ->0 )
    # end_1
    
    f.loc[f['splice_type'] == 1 , 'end_1'] = u_dist 
    f.loc[f['splice_type'] == 2 , 'end_1'] = u_dist + intron_1 + exon_se
    f.loc[f['splice_type'] == 3 , 'end_1'] = u_dist 
    f.loc[f['splice_type'] == 4 , 'end_1'] = u_dist
    f.loc[f['splice_type'] == 5 , 'end_1'] = f['lengths'] 
    
    
    # start_2
    f.loc[f['splice_type'] == 1 , 'start_2'] = u_dist+intron_1
    f.loc[f['splice_type'] == 2 , 'start_2'] = u_dist+intron_3
    f.loc[f['splice_type'] == 3 , 'start_2'] = u_dist+intron_1 
    f.loc[f['splice_type'] == 4 , 'start_2'] = u_dist+intron_1+exon_se+intron_2
    
    # end_2
    f.loc[f['splice_type'] == 1 , 'end_2'] = f['lengths']+intron_1
    f.loc[f['splice_type'] == 2 , 'end_2'] = f['lengths']+intron_2
    f.loc[f['splice_type'] == 3 , 'end_2'] = u_dist+intron_1+exon_se 
    f.loc[f['splice_type'] == 4 , 'end_2'] = f['lengths']+intron_3
    # start_3
    f.loc[f['splice_type'] == 3 , 'start_3'] = u_dist+intron_1+exon_se +intron_2
    # end_3
    f.loc[f['splice_type'] == 3 , 'end_3'] = f['lengths']+intron_1+intron_2#u_dist+intron_1+exon_se +intron_2+exon
    
    
    f.insert(0,"u_dist",len(f)*[u_dist])
    f.insert(0,"d_dist",len(f)*[d_dist])
    f.insert(0,"intron_3",len(f)*[intron_3])
    f.insert(0,"exon",len(f)*[exon])
    f.insert(0,"intron_1",len(f)*[intron_1])
    f.insert(0,"intron_2",len(f)*[intron_2])
    f.insert(0,"exon_se",len(f)*[exon_se])
    f.insert(0,"psi_se",len(f)*[psi_se])
    
    ## Get junction read data
    jnc_df = pandas.DataFrame({"start position":e['read_start'],
                               "junction read":e['junction'],
                               "read length":50,
                                })
    e.to_csv (os.getcwd()+'/sim_outputs'+'/export_dataframe_'+str(file_num)+'.csv', index = None, header=True, sep=',')
    f.to_csv (os.getcwd()+'/sim_outputs'+'/export_full_transcripts_'+str(file_num)+'.csv', index = None, header=True, sep=',')
    ct_df.to_csv (os.getcwd()+'/sim_outputs'+'/export_start_posns_'+str(file_num)+'.csv', index = None, header=True, sep=',')
    jnc_df.to_csv (os.getcwd()+'/sim_outputs'+'/export_junction_reads_'+str(file_num)+'.csv', index = None, header=True, sep=',')
#    TODO: clear memory of dataframes
#    del e
#    del f
#    del ct_df
#    del jnc_df
    file_num+=1

                            
#df = pandas.DataFrame(ls,index = None, columns = None)
#df = df.transpose()
#df.to_csv("output.csv", sep=',')

