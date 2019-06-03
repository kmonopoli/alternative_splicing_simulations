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
exon = 0.3
u_dist = 0.5
n_millions = 100 # total number of millions of transcripts to consider
transc_rate = 1.5 # rate of transcription
# variables (to simulate over)
labelings = [5]#[5,10,20,60]
#h_s = list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
h_s = [0.2]#[100]#[0.2]
introns =[0.04] #[5]#[0.04]#[40.0] #NOTE: need to make larger because alt spliced introns are larger (look up)
#introns =list(np.arange(0.04,0.09,0.02)+list(np.arange(0.1,1,0.1))+list(np.arange(1,50,2))

#d_dists = list(np.arange(0.5, 5, 0.5))
#expr_lvls = list(np.arange(1,100, 4))

# for alt splicing TODO: just for testing need to change to simulate over
#                   NOTE: also in future need to pick actual values because now these are determined by intron length and anything 
#                       that isn't divisible by 4 will give a non-integral response which won't work in the simulations
#h_s_alt_e = list(np.arange(0.2,0.9,0.1))+list(np.arange(1,10,0.75))+list(np.arange(11,100,2))
h_s_alt_e = [0.2]
exon_se = 0# introns[0]/4
intron_1 = 0#(introns[0]-exon_se)/2
intron_2 = 0#introns[0]-(exon_se+intron_1)
psi_se=[0]#[0.5]#list(np.arange(0.0,1,0.1)) #Psi of SE (skipped exon) ## TODO: need a function that calculates this



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
    df = pandas.DataFrame({'end_sites':end_sites})

   

    # Determine if the transcripts are spliced or not
    # And determine resulting lengths of transcripts that were spliced
    spliced = [splice(x, intron, intron_1, exon_se, u_dist, h, h_alt_e, transc_rate, psi_se) for x in end_sites]
    #spliced = spliced[0:10] ## TODO: TESTING remove [0:10] later
    df.insert(0, "spliced",spliced, True)

    
    
    # Get the reads from the transcripts 
    reads = [get_reads(x[0]) for x in spliced] 
    df.insert(0, "reads",reads, True)
    return df
'''
    ## add splice type to all elements with same splice type
    splice_types =[ x[1] for x in spliced]
    splice_counts =  [len(x[0]) for x in reads]
    sp = []
    for x, y in zip(splice_types, splice_counts):
        sp.append([x]*y)
    sp = [item for sublist in [x for x in sp] for item in sublist]
    starts = [item for sublist in [x[0] for x in reads] for item in sublist]
    lengths = [item for sublist in [x[1] for x in reads] for item in sublist]

    # create dataframe with all of the data
    df = pandas.DataFrame({'start_pos':starts,'length':lengths,'splice_type':sp})
    # Map reads to the gene
    map_data(df, u_dist, intron, intron_1, intron_2)

    # Identify junction reads
    rl = 50
    get_junction_reads(df,rl)
    return df
'''

# determines if a read is a junction read (over intron/exon boundary or exon/exon boundary)
# Adds info to dataframe: True if junction read, False otherwise
def get_junction_reads(df,rl):
    i = 0
    junction_rds = []
    while i<df.shape[0]:
        s = df['start_site'].iloc[i]
        spl_t = df['splice_type'].iloc[i] # splice type
        j_read = False
        if (spl_t == 0 or spl_t == 4): # unspliced
            # can get 4 types of junction reads:
            #   exon(u_dist)/intron_1
            #   intron_1/exon_se
            #   exon_se/intron_2
            #   intron_2/exon
            ## TODO: this is complicated and will depend on the length
            j_read = False

        elif(spl_t == 1): # spliced normally
            # only can get junction read of exon(u_dist)/exon
            j_read = (((s+rl) > 0) and ((s+rl) < rl))
            

        elif(spl_t == 2): # spliced alternatively
            # can get 2 types of junction reads:
            #   exon(u_dist)/exon_se 
            #   exon_se/exon
            j_read = False

        
        elif(spl_t == 3): # spliced alternatively but not done transcribing intron_2
            # can get 2 types of junction reads:
            #   exon(u_dist)/exon_se
            #   exon_se/intron_2
            j_read = False

        else:
            print "ERROR splice type not supported: ",spl_t
            quit()
        junction_rds.append(j_read)
        i+=1
    # add junction read info to dataframe
    df.insert(0, "junction_read",junction_rds, True)
    
    
    

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
#   data = [lengths, splice_type]
#       lengths - the length of a transcript
#       splice type
#
# Outputs:
#   fragments - an array of arrays ith the start positions of the filtered fragments and the index
#               of the transcript that it came from (columns: transcript, start)
def get_reads(data):
   

    lengths = [x[0] for x in data]
    splice_type = [x[1] for x in data]
    data = [lengths,splice_type]
    cwd = os.getcwd()
    process = Popen([cwd+"/weibull_dist_calc.R",str(data)], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    stdout = [map(int, q) for q in [x.replace(")","").replace("\n","").replace(" ","").replace('"','').split(",") for x in stdout.split("c(")][1:]]
    return stdout
















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
                            # = simulate(intron, exon, u_dist, d_dist, labeling, h, expr_lvl, n_millions, transc_rate, intron_1, intron_2, exon_se, h_alt_e, psi)
                            #e.to_csv (os.getcwd()+'/export_dataframe.csv', index = None, header=True, sep=',')
                            #ls.append(e)
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
                            
                            # Get the reads from the transcripts 
                            reads = get_reads(spliced)
                            start_pos = pandas.DataFrame({"transcript":reads[0],"reads":reads[1]})
                            
                            # Get splice types (numeric)
                            start_pos.insert(0,"splice_type",[spliced[x-1][1] for x in start_pos['transcript']])
                            
                            
                            
                            '''        
  reads = data.frame(start = start_pos$start - U_dist +
     intron*start_pos$spliced*(start_pos$start > U_dist))
  reads$name = paste('i', as.character(intron), 'e', as.character(exon),
                     'u', as.character(U_dist), 'd', as.character(D_dist),
                     'L', as.character(labeling), 'hl', as.character(half_life),
                     'X', as.character(expression_level),
                     'M', as.character(n_millions), sep=':')
  reads_formatted = formatreads(reads, intron, start_pos, U_dist, rl)
  # reads$match = '5M'
  # #junction = (reads$start > -rl + 1) & (reads$start <= 0) & start_pos$spliced
  reads$junction = (reads$start > -rl + 10) & (reads$start <= -9) & start_pos$spliced
  # reads$match[junction] = 
  #   paste(as.character(-reads$start[junction]), 'M', as.character(intron), 'N', 
  #         as.character(rl+reads$start[junction]), 'M', sep = '')
  # reads$name = paste('intron', as.character(intron), 'exon', as.character(exon), 
  #                    'U_dist', as.character(U_dist), 'D_dist', as.character(D_dist),
  #                    'labeling', as.character(labeling), 'half_life', as.character(half_life), 
  #                    'expression_level', as.character(expression_level), 
  #                    'n_millions', as.character(n_millions), sep=':')
  spliced_num = sum(reads$junction)
  unspliced_num = sum(!start_pos$spliced &
                        (start_pos$start > intron + U_dist-rl+10) &
                        (start_pos$start <= intron + U_dist-9))
  intron_num = sum(!start_pos$spliced &
                     (start_pos$start >= 1) & (start_pos$start < 50))
                              '''


                            
                         
df = pandas.DataFrame(ls,index = None, columns = None)
df = df.transpose()
df.to_csv("output.csv", sep=',')
