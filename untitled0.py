#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 06:21:16 2019

@author: kmonopoli
"""

#remove PCR duplicates and anything <18nt
#for i in $(find . -name "*.q20p100"); do /home/km39w/scripts/ms_short_64G "awk '(NR%4==1){name=\$1}(NR%4==2){total++;if((substr(\$1,length(\$1)-11,3)==\"GTC\")&&(substr(\$1,length(\$1)-5,3)==\"TAG\")&&(((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\"))||((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\")))){umis++;if(a[\$1]!=1){nondup++;a[\$1]=1;if(length(\$1)>47){longer18++; print name; print substr (\$1,16,length (\$1)-30);getline; print; getline; print substr (\$1,16,length (\$1)-30)}}}}END{print FILENAME\"\\t\"total\"\\t\"umis\"\\t\"nondup\"\\t\"longer18 > FILENAME\".dup\"}' $i > $i.deUMI.dedup.fq" ; done

dr= '/home/km39w/scripts/ms_short_64G'
ls = [x for x in os.listdir('.') if x[-8:] == ".q20p100"]
wk = "\ "[0]
ls = [ls[0]]
for fl in ls:
    bashCommand = /home/km39w/scripts/ms_short_64G "awk '(NR%4==1){name=\$1}(NR%4==2){total++;if((substr(\$1,length(\$1)-11,3)==\"GTC\")&&(substr(\$1,length(\$1)-5,3)==\"TAG\")&&(((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\"))||((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\")))){umis++;if(a[\$1]!=1){nondup++;a[\$1]=1;if(length(\$1)>47){longer18++; print name; print substr (\$1,16,length (\$1)-30);getline; print; getline; print substr (\$1,16,length (\$1)-30)}}}}END{print FILENAME\"\\t\"total\"\\t\"umis\"\\t\"nondup\"\\t\"longer18 > FILENAME\".dup\"}' $i > $i.deUMI.dedup.fq"
    
    
    
    
    
    
    print bashCommand