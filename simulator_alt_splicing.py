#!/anaconda2/bin/python
#import random
#import subprocess
#import matplotlib.pyplot as plt
import numpy as np
import pandas
import os
import math
import gc
from subprocess import Popen, PIPE




## attributes TODO: make so read from command line
d_dists = [5.0]#[0.5]#list(np.arange(0.5, 5, 0.5))
expr_lvls = [100.0]#[100]#[4] #list(np.arange(1,100, 4))

## parameters
# constants
exon = 3.0#5.0#0.3
u_dist = 5.0#5.0#0.5
n_millions = 100.0 # total number of millions of transcripts to consider
transc_rate = 5.0#1.5 # rate of transcription


# variables (to simulate over)
labelings = [5.0]#,10,20]#[5,10,20,60]
introns_3 =[10.0]#[15.0]#[40.0] #list(np.arange(0.04,0.09,0.02)+list(np.arange(0.1,1,0.1))+list(np.arange(1,50,2)) #NOTE: need to make larger because alt spliced introns are larger (look up)

hs_intron_1 = [40.0]#[40.0] # list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
hs_intron_2 = [10.0]#[40.0] # list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
hs_intron_3 = [10.0]#[100]#[0.2] # list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))


# for alt splicing 
#   TODO: just for testing need to change to simulate over
#   NOTE: also in future need to pick actual values because now these are determined by intron length and anything 
#   that isn't divisible by 4 will give a non-integral response which won't work in the simulations
exon_se = introns_3[0]/4

introns_1 = [(introns_3[0]-exon_se)/2]
introns_2 = [introns_3[0]-(exon_se+introns_1[0])]

d_dists_1 = [introns_2[0] + exon + d_dists[0]]

psi_se_s=[0.5]#list(np.arange(0.0,1,0.1)) #Psi of SE (skipped exon) ## TODO: need a function that calculates this






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
#     [  u_dist   ]----intron_1----[/// exon_se ///]---intron_2----[  exon ][ (d_dist)  ][AAAAAA (L*r)]
#                                  <---------------- d_dist_1 -------------------------->
#
#
#
def splice(end_site, intron_3, intron_1, intron_2, exon_se, u_dist, d_dist, h_intron_3, h_intron_1, h_intron_2, transc_rate, psi_se, d_dist_1): 
    t = -1 # splice type
    final_len = 0 # length of transcript after splicing
    spl_path = np.random.choice([0, 1], p=[1-psi_se,psi_se]) # splice path - whether it proceed down path of alt splicing (include exon_se), or normal splicing (exclude exon_se)
    # spl_path will be 1 if alt splicing (include exon_se), 0 if normal splicing (exclude exon_se)
    
    # Alternative Splicing (include exon_se)
    if(spl_path == 1): 
        if(end_site < intron_1+u_dist): # too short to splice
            t = 5 # unspliced
            final_len = u_dist+end_site
    
        # short - exclude intron_1 / unspliced
        elif((end_site >= intron_1+u_dist) and (end_site < intron_3+u_dist)):
            n = np.random.uniform(0,1,1).tolist()[0]>2**(-1*(end_site-intron_1)/h_intron_1/transc_rate)
            if n: # intron_1 spliced out
                t = 1
                #OLD:
#                final_len = u_dist + min(end_site, intron_1+exon_se+d_dist_1) - intron_1
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
                #OLD: 
#                final_len = u_dist + min(end_site, intron_1+intron_2+exon_se+d_dist) - (intron_1+intron_2)
                final_len = u_dist + end_site-(intron_1+intron_2)
            elif m: # exclude intron_1 
                t = 1
                #OLD: 
#                final_len = u_dist + min(end_site, intron_1+intron_2+exon_se+d_dist_1) - (intron_1)
                final_len = u_dist + end_site-(intron_1)
            elif n: # exclude intron_2
                t = 2
                #OLD: 
#                final_len = u_dist + min(end_site, intron_1+intron_2+exon_se+d_dist) - (intron_1)
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
#   fragments - an array of arrays with the start positions of the filtered fragments and the index
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


def sum_square_eqn_solver(d_prime, r_prime, transcription_rate):
    data = [d_prime, r_prime, transcription_rate]
#    print data
    cwd = os.getcwd()
    process = Popen([cwd+"/sum_sq_equation_solver.R",str(data)], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    
    stdout = stdout.replace('[1] "running sum square equation solver"'+'\n[1] "','')[:-2].split("   ")
    stdout = [float(x) for x in stdout]
    if( stdout == []):
        print "ERROR in get_reads function that calls sum_sq_equation_solver.R -  no output received"
        #quit()
    return stdout









## START 


ls = []

################################## ################################## ##################################
## FOR TESTING: TODO: iterate through in real simulation
d_dist = d_dists[0]
d_dist_1 = d_dists_1[0]
expr_lvl = expr_lvls[0]
#labeling = labelings[0]
h_intron_3 = hs_intron_3[0]
intron_3 = introns_3[0]
h_intron_1 = hs_intron_1[0]
h_intron_2 = hs_intron_2[0]
psi_se = psi_se_s[0]
lebeling = labelings[0]
intron_1 = introns_1[0]
intron_2 = introns_2[0]


file_num=1

intron_3 = 1000*intron_3
exon   = 1000*exon
u_dist = 1000*u_dist

transc_rate = 1000*transc_rate
intron_1 = 1000*intron_1
intron_2 = 1000*intron_2
exon_se = 1000*exon_se

rl = 50

# TODO: need to run these for loops to go through all possibilities
for labeling in labelings:
#for d_dist in d_dists:
#        for expr_lvl in expr_lvls:
#            for intron_1 in introns_1: 
#                for intron_2 in introns_2:
#                    for intron_3 in introns_3:
#                        for h_intron_1 in hs_intron_1:
#                            for h_intron_2 in hs_intron_2:
#                                for h_intron_3 in hs_intron_3:
#                                    for psi_se in psi_se_s: 


    ######################### SIMULATION START #########################
    ### Simulation:
    ##  Outputs: 
    ##   [[1]] data frame of all reads with the start position, the matching pattern, and a string with all the associated parameters (columns: start, match, name)
    ##   [[2]] the number of spliced junction reads
    ##   [[3]] the number of unspliced junction reads
    ##   [[4]] whether or not each read comes from a spliced or unspliced transcript
    print "Labeling: ",labeling
    d_dist = 1000*d_dist
    d_dist_1 = 1000*d_dist_1
    # Generate expression_level*n_millions transcripts uniformly from the labeled region 
    end_sites =  np.random.choice(range(1,int(intron_3 + exon + d_dist + labeling*transc_rate)), int(expr_lvl*n_millions))
    end_sites =[u_dist+x for x in end_sites]

    # Determine if the transcripts are spliced or not
    # And determine resulting lengths of transcripts that were spliced
    spliced = [splice(x, intron_3, intron_1, intron_2, exon_se, u_dist, d_dist, h_intron_3, h_intron_1, h_intron_2, transc_rate, psi_se, d_dist_1) for x in end_sites]
    # transpose data and put in dataframe
    spliced_df = pandas.DataFrame({"lengths":[x[0] for x in spliced],"splice_type":[x[1] for x in spliced]})
    
    # Get the reads from the transcripts and map them to the gene
    start_reads = get_reads(spliced) 
    start_pos = pandas.DataFrame({"transcript":start_reads[0],"start":start_reads[1]})
    
    # Get splice types (numeric)
    start_pos.insert(0,"splice_type",[spliced[x-1][1] for x in start_pos['transcript']])

        

    # Get read start positions
    # TODO: there is a mistake here, it skips intron_1
    start_pos["read_start"] = np.nan
    #    1 --> intron_1 excluded
    start_pos.read_start.loc[(start_pos['splice_type'] == 1)] = (start_pos['start'] - u_dist) + intron_1 * (start_pos['start'] > u_dist) 
    #    2 --> intron_2 excluded
    start_pos.read_start.loc[(start_pos['splice_type'] == 2)] = (start_pos['start'] - u_dist) + 1 * (start_pos['start'] > u_dist) + intron_2 * (start_pos['start'] > u_dist+intron_1+exon_se) 
    #    3 --> intron_1 and intron_2 excluded (exon_se included)
    start_pos.read_start.loc[(start_pos['splice_type'] == 3)] = (start_pos['start'] - u_dist) + intron_1 * (start_pos['start'] > u_dist) + intron_2 * (start_pos['start'] > u_dist+exon_se)
    #    4 --> intron_3 excluded
    start_pos.read_start.loc[(start_pos['splice_type'] == 4)] = (start_pos['start'] - u_dist) + intron_3 * (start_pos['start'] > u_dist)
    #    5 --> unspliced
    start_pos.read_start.loc[(start_pos['splice_type'] == 5)] = (start_pos['start'] - u_dist) #+ intron_1 * (start_pos['start'] > u_dist)

    ######################### SIMULATION END #########################

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
    
#    start_pos.junction.loc[(start_pos['splice_type'] == 1) & # intron_1 excluded
#                           (start_pos['read_start']  > ((intron_1 + exon_se)+(-1*rl) +junct_read_overlap)) &
#                           (start_pos['read_start']  <= ((intron_1 + exon_se)+(-1*junct_read_overlap)))] = "exon_se - intron_2"
                           
    start_pos.junction.loc[(start_pos['splice_type'] == 1) & # intron_1 excluded
                           (start_pos['read_start']  > ((intron_1 + exon_se + intron_2)+(-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <= ((intron_1 + exon_se + intron_2)+(-1*junct_read_overlap)))] = "intron_2 - exon"
    
#    start_pos.junction.loc[(start_pos['splice_type'] == 2) & # intron_2 excluded
#                           (start_pos['read_start']  > ((-1*rl) + junct_read_overlap)) &
#                           (start_pos['read_start']  <=(-1*junct_read_overlap))] = "exon_us - intron_1"
    start_pos.junction.loc[(start_pos['splice_type'] == 2) & # intron_2 excluded
                           (start_pos['read_start']  > (intron_1 +(-1*rl) + junct_read_overlap)) &
                           (start_pos['read_start']  <=(intron_1 +(-1*junct_read_overlap)))] = "intron_1 - exon_se"
                           
#    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
#                           (start_pos['read_start']  > ((-1*rl) + junct_read_overlap)) &
#                           (start_pos['read_start']  <=(-1*junct_read_overlap))] = "exon_us - intron_3"
    
# introns 1 and 2 as well as exon_se do not exist for this splice type?

#    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
#                           (start_pos['read_start']  > ((intron_1 + exon_se)+(-1*rl) +junct_read_overlap)) &
#                           (start_pos['read_start']  <= ((intron_1 + exon_se)+(-1*junct_read_overlap)))] = "exon_se - intron_2"


#    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
#                           (start_pos['read_start']  > (intron_1 +(-1*rl) + junct_read_overlap)) &
#                           (start_pos['read_start']  <=(intron_1 +(-1*junct_read_overlap)))] = "intron_1 - exon_se"
#                           
    # OLD: whether or not the read is an intron_3 - exon junction read or an intron_2 - exon junction read depends on psi:
    # whether or not the read is an intron_3 - exon junction read or an intron_2 - exon junction read depends on ratio of splicing in intron_1 vs intron_3:


    spliced_num_int1 = float(len(start_pos[start_pos.junction == "exon_us - exon_se"]))
    spliced_num_int2 = float(len(start_pos[start_pos.junction == "exon_se - exon"])) # used later in half-life calculation
    spliced_num_int3 = float(len(start_pos[start_pos.junction == "exon_us - exon"]))
    
    
    if (spliced_num_int1 == 0 ):
        spl = 0.0
    elif (spliced_num_int3 == 0):
        spl = 1.0
    else:
        spl = spliced_num_int1/spliced_num_int3

    print "SPL --> ",spl
    # OLD (based on psi):spl_path = np.random.choice([0, 1], p=[1-psi_se,psi_se]) # splice path - whether it proceed down path of alt splicing (include exon_se), or normal splicing (exclude exon_se)
    spl_path = np.random.choice([0, 1], p=[1-spl,spl]) # splice path - whether it proceed down path of alt splicing (include exon_se), or normal splicing (exclude exon_se)

    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
                           (spl_path == 0) & 
                           (start_pos['read_start']  > ((intron_3)+(-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <= ((intron_3)+(-1*junct_read_overlap)))] = "intron_3 - exon"
    
                           
    start_pos.junction.loc[(start_pos['splice_type'] == 5) & # unspliced
                           (spl_path == 1) & 
                           (start_pos['read_start']  > ((intron_1 + exon_se + intron_2)+(-1*rl) +junct_read_overlap)) &
                           (start_pos['read_start']  <= ((intron_1 + exon_se + intron_2)+(-1*junct_read_overlap)))] = "intron_2 - exon"

    # used later in half-life calculation
    unspliced_num_int1 = float(len(start_pos[start_pos.junction == "intron_1 - exon_se"]))
    unspliced_num_int2 = float(len(start_pos[start_pos.junction == "intron_2 - exon"]))
    unspliced_num_int3 = float(len(start_pos[start_pos.junction == "intron_3 - exon"]))
    
                       
    # add junction read type (ie/ei and ee)
    start_pos["junc_read_type"] = "n/a"
    start_pos.junc_read_type.loc[(start_pos['junction'].str.contains("- intron"))] = "ei"
    start_pos.junc_read_type.loc[(start_pos['junction'].str.contains("intron_1 - |intron_2 - ", regex = True))] = "ie"
    start_pos.junc_read_type.loc[(start_pos['junc_read_type'] == "n/a")] = "ee"
    
    

    
    
    
    
    ## add simulation parameter data for test plots
    start_pos.insert(0,"h_intron_1",len(start_pos)*[h_intron_1])
    start_pos.insert(0,"h_intron_2",len(start_pos)*[h_intron_2])
    start_pos.insert(0,"h_intron_3",len(start_pos)*[h_intron_3]) # normal splicing
    start_pos.insert(0,"transc_rate",len(start_pos)*[transc_rate])
    start_pos.insert(0,"intron_3",len(start_pos)*[intron_3])
    start_pos.insert(0,"intron_1",len(start_pos)*[intron_1])
    start_pos.insert(0,"intron_2",len(start_pos)*[intron_2])
    start_pos.insert(0,"u_dist",len(start_pos)*[u_dist])
    start_pos.insert(0,"exon",len(start_pos)*[exon])
    start_pos.insert(0,"exon_se",len(start_pos)*[exon_se])
    start_pos.insert(0,"psi_se",len(start_pos)*[psi_se])
    ## Relabel splice types to readable words
    start_pos.insert(0,"splice_type_read",start_pos['splice_type'])
    start_pos['splice_type_read'] = np.where(start_pos['splice_type_read'] == 1, "intron_1 excluded", start_pos['splice_type_read'])
    start_pos['splice_type_read'] = np.where(start_pos['splice_type_read'] == '2', "intron_2 excluded", start_pos['splice_type_read'])
    start_pos['splice_type_read'] = np.where(start_pos['splice_type_read'] == '3', "intron_1 & intron_2 excluded", start_pos['splice_type_read'])
    start_pos['splice_type_read'] = np.where(start_pos['splice_type_read'] == '4', "intron_3 excluded", start_pos['splice_type_read'])
    start_pos['splice_type_read'] = np.where(start_pos['splice_type_read'] == '5', "unspliced", start_pos['splice_type_read'])
    
    
    # Counts
    mn,mx =(round(start_pos['read_start'].min()/100,0)*100, round( start_pos['read_start'].max()/100,0)*100+100)
    
    start_pos.loc[start_pos['splice_type_read'] == 'intron_1 excluded','read_start'].value_counts()
    
    
    ct_df = pandas.DataFrame({"ct_intron_1_excluded":start_pos.loc[start_pos['splice_type_read'] == 'intron_1 excluded','read_start'].value_counts(),
                                "ct_intron_2_excluded":start_pos.loc[start_pos['splice_type_read'] == 'intron_2 excluded','read_start'].value_counts(),
                                "ct_intron_1_intron_2_excluded":start_pos.loc[start_pos['splice_type_read'] == 'intron_1 & intron_2 excluded','read_start'].value_counts(),
                                "ct_intron_3_excluded":start_pos.loc[start_pos['splice_type_read'] == 'intron_3 excluded','read_start'].value_counts(),
                                "ct_unspliced":start_pos.loc[start_pos['splice_type_read'] == 'unspliced','read_start'].value_counts()
                                })
    ct_df.insert(0,"index",ct_df.index)
        

    
    
    spliced_df.insert(0,"splice_type_read",spliced_df['splice_type'])
    spliced_df['splice_type_read'] = np.where(spliced_df['splice_type_read'] ==  1, "intron_1 excluded", spliced_df['splice_type_read'])
    spliced_df['splice_type_read'] = np.where(spliced_df['splice_type_read'] == '2', "intron_2 excluded", spliced_df['splice_type_read'])
    spliced_df['splice_type_read'] = np.where(spliced_df['splice_type_read'] == '3', "intron_1 & intron_2 excluded", spliced_df['splice_type_read'])
    spliced_df['splice_type_read'] = np.where(spliced_df['splice_type_read'] == '4', "intron_3 excluded", spliced_df['splice_type_read'])
    spliced_df['splice_type_read'] = np.where(spliced_df['splice_type_read'] == '5', "unspliced", spliced_df['splice_type_read'])
    
    ## get starts and ends for each
    spliced_df.insert(0,"end_3",None)
    spliced_df.insert(0,"start_3",None)
    spliced_df.insert(0,"end_2",None)
    spliced_df.insert(0,"start_2",None)
    spliced_df.insert(0,"end_1",None)
    spliced_df.insert(0,"start_1",0) # always starts at 0
    # start_1 
    #    (will be same for all ->0 )
    # end_1
    spliced_df.loc[spliced_df['splice_type'] == 1 , 'end_1'] = u_dist 
    spliced_df.loc[spliced_df['splice_type'] == 2 , 'end_1'] = u_dist + intron_1 + exon_se
    spliced_df.loc[spliced_df['splice_type'] == 3 , 'end_1'] = u_dist 
    spliced_df.loc[spliced_df['splice_type'] == 4 , 'end_1'] = u_dist
    spliced_df.loc[spliced_df['splice_type'] == 5 , 'end_1'] = spliced_df['lengths'] 
    # start_2
    spliced_df.loc[spliced_df['splice_type'] == 1 , 'start_2'] = u_dist+intron_1
    spliced_df.loc[spliced_df['splice_type'] == 2 , 'start_2'] = u_dist+intron_3
    spliced_df.loc[spliced_df['splice_type'] == 3 , 'start_2'] = u_dist+intron_1 
    spliced_df.loc[spliced_df['splice_type'] == 4 , 'start_2'] = u_dist+intron_1+exon_se+intron_2
    # end_2
    spliced_df.loc[spliced_df['splice_type'] == 1 , 'end_2'] = spliced_df['lengths']+intron_1
    spliced_df.loc[spliced_df['splice_type'] == 2 , 'end_2'] = spliced_df['lengths']+intron_2
    spliced_df.loc[spliced_df['splice_type'] == 3 , 'end_2'] = u_dist+intron_1+exon_se 
    spliced_df.loc[spliced_df['splice_type'] == 4 , 'end_2'] = spliced_df['lengths']+intron_3
    # start_3
    spliced_df.loc[spliced_df['splice_type'] == 3 , 'start_3'] = u_dist+intron_1+exon_se +intron_2
    # end_3
    spliced_df.loc[spliced_df['splice_type'] == 3 , 'end_3'] = spliced_df['lengths']+intron_1+intron_2#u_dist+intron_1+exon_se +intron_2+exon
    
    # add simulation parameter info
    spliced_df.insert(0,"u_dist",len(spliced_df)*[u_dist])
    spliced_df.insert(0,"d_dist",len(spliced_df)*[d_dist])
    spliced_df.insert(0,"intron_3",len(spliced_df)*[intron_3])
    spliced_df.insert(0,"exon",len(spliced_df)*[exon])
    spliced_df.insert(0,"intron_1",len(spliced_df)*[intron_1])
    spliced_df.insert(0,"intron_2",len(spliced_df)*[intron_2])
    spliced_df.insert(0,"exon_se",len(spliced_df)*[exon_se])
    spliced_df.insert(0,"psi_se",len(spliced_df)*[psi_se])
    spliced_df.insert(0,"d_dist_1",len(spliced_df)*[d_dist_1])
    

    # Get junction read data
    jnc_df = pandas.DataFrame({"start position":start_pos['read_start'],
                               "junction_read":start_pos['junction'],
                               "junction_read_type":start_pos['junc_read_type'],
                               "read length":50,
                                })
    # remove exon-intron junction reads (not useful for analysis) 
    jnc_df = jnc_df.drop(jnc_df[jnc_df.junction_read_type == "ei"].index)
    
    
    

    # Estimate half-lives from Junction Reads
    # TODO put calculations after loop (so can use multiple calculations to calculate)
    # initialize estimated halflives  (DO NOT REMOVE)
    est_h_int1 = "N/A"
    est_h_int2 = "N/A"
    est_h_int3 = "N/A"
    
    
    # intron 1
    if(spliced_num_int1 == 0 or unspliced_num_int1 == 0):
        print "\nWARNING: spliced_num_int1 or  is zero, estimated halflife for intron 1 will be incorrect"
    else:
        # compute d_prime_int1
        #   d_prime_int1 can be two values depending on splicing of intron_2:
        #    ratio_int2N >= 1 means more unspliced than spliced so intron_2 will be more likely to be INCLUDED in d_prime_int1 calculation
        #    ratio_int2N <  1 means more spliced than unspliced so intron_2 will be more likely to be EXCLUDED in d_prime_int1 calculation

#        if (spliced_num_int2 == 0 and unspliced_num_int2 ==0): # no intron_2 reads - treats as all unspliced
#            ratio_int2N = 1.0
#        elif (spliced_num_int2 == 0): # all unspliced
#            ratio_int2N = 1.0
#        elif  (unspliced_num_int2 ==0) : # all spliced
#            ratio_int2N = 0.0
#        else:
#            ratio_int2N = unspliced_num_int2/spliced_num_int2
#            # make into a ratio
#            if ratio_int2N > 1:
#                ratio_int2N = 1.0-(1.0/ratio_int2N)
#        spl_int_2 = np.random.choice([0, 1], p=[ratio_int2N,1-ratio_int2N]) # 0 -> intron_2 included (unspliced), 1 -> intron_2 excluded (spliced)
#        # TODO be sure the 1 and 0 aren't flipped!
#        # intron_2 included (unspliced)
#        if spl_int_2 == 1:
#            d_prime_int1 =  (spliced_df.exon_se[0] + spliced_df.intron_2[0] + spliced_df.exon[0] + spliced_df.d_dist[0])+(float(labeling)*float(transc_rate))
#        # intron_2 excluded (spliced)
#        elif spl_int_2 == 0:
#            d_prime_int1 =  (spliced_df.exon_se[0] + spliced_df.d_dist_1[0])+(float(labeling)*float(transc_rate))

        d_prime_int1 =  (spliced_df.exon_se[0] + spliced_df.d_dist_1[0])+(float(labeling)*float(transc_rate))
        ratio_int1 = unspliced_num_int1/spliced_num_int1
        r_prime_int1 = (d_prime_int1*math.log(2))/(float(transc_rate)*((1.0/ratio_int1)+1.0))
        est_h_int1 = sum_square_eqn_solver(d_prime_int1,r_prime_int1,transc_rate)[0]
        
    # intron 2
    
    if(spliced_num_int2 == 0 or unspliced_num_int2 ==0 ):
        print "\nWARNING: spliced_num_int2 or  is zero, estimated halflife for intron 2 will be incorrect"
    else:
        ratio_int2 = unspliced_num_int2/spliced_num_int2
        d_prime_int2 =  (spliced_df.exon[0] + spliced_df.d_dist[0])+(float(labeling)*float(transc_rate))
        r_prime_int2 = (d_prime_int2*math.log(2))/(float(transc_rate)*((1.0/ratio_int2)+1.0))
        est_h_int2 = sum_square_eqn_solver(d_prime_int2,r_prime_int2,transc_rate)[0]

    # intron 3

    if(spliced_num_int3 == 0 or unspliced_num_int3 ==0):
        print "\nWARNING: spliced_num_int3 or  is zero, estimated halflife for intron 3 will be incorrect"

    else:
        ratio_int3 = unspliced_num_int3/spliced_num_int3
        d_prime_int3 =  (spliced_df.exon[0] + spliced_df.d_dist[0])+(float(labeling)*float(transc_rate))
        r_prime_int3 = (d_prime_int3*math.log(2))/(float(transc_rate)*((1.0/ratio_int3)+1.0))
        est_h_int3 = sum_square_eqn_solver(d_prime_int3,r_prime_int3,transc_rate)[0]
    print spliced_num_int3
    print unspliced_num_int3
#    # corrections for introns 1 and 2
#    corr_int1 = 1/((unspliced_num_int2+spliced_num_int2) / (unspliced_num_int1+spliced_num_int1+unspliced_num_int2+spliced_num_int2))
#    corr_int2 = (unspliced_num_int1+spliced_num_int1) / (unspliced_num_int2+spliced_num_int2+unspliced_num_int1+spliced_num_int1)
#    
#    print "before corr 1: ", est_h_int1
#    print "before corr 2: ", est_h_int2
#
#    
#    est_h_int1 = est_h_int1*corr_int1
#    est_h_int2 = est_h_int2*corr_int2
        
    h_est_df = pandas.DataFrame({"intron":["intron_1","intron_2","intron_3"],"h_actual":[h_intron_1,h_intron_2,h_intron_3],"h_estimated":[est_h_int1,est_h_int2,est_h_int3]})
    
    

    # TODO need to name these files using parameter values, rather than file_num
    start_pos.to_csv (os.getcwd()+'/sim_outputs'+'/export_dataframe_'+str(file_num)+'.csv', index = None, header=True, sep=',')
    spliced_df.to_csv (os.getcwd()+'/sim_outputs'+'/export_full_transcripts_'+str(file_num)+'.csv', index = None, header=True, sep=',')
    ct_df.to_csv (os.getcwd()+'/sim_outputs'+'/export_start_posns_'+str(file_num)+'.csv', index = None, header=True, sep=',')
    jnc_df.to_csv (os.getcwd()+'/sim_outputs'+'/export_junction_reads_'+str(file_num)+'.csv', index = None, header=True, sep=',')
    h_est_df.to_csv (os.getcwd()+'/sim_outputs'+'/export_half-lives_'+str(file_num)+'.csv', index = None, header=True, sep=',')
    print h_est_df
#    # clear dataframes from memory
#    del ct_df
#    del jnc_df
#    del spliced_df
#    del start_pos
#
#    gc.collect()

    file_num+=1
                            
#df = pandas.DataFrame(ls,index = None, columns = None)
#df = df.transpose()
#df.to_csv("output.csv", sep=',')






















