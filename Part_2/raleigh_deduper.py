#!/usr/bin/env python

import argparse
import re
from unittest import skip

def get_args():
    parser = argparse.ArgumentParser(description="This script will remove PCR duplicates and UMI errors")
    parser.add_argument("-f", help="file to run", type=str, required = True)
    parser.add_argument("-o", help="output file", type=str, required = True)
    parser.add_argument("-u", help="File containing UMI's", type=str, required = True)
    return parser.parse_args()

parameters = get_args()
f = parameters.f
o = parameters.o
u = parameters.u

umis = []
# a list to store the UMIs
chrom_dict = {}
# a dictionary to store the reads for each chromosome
dict_key = ()
# the variable that hold the key for each line in the chrom dictionary

dup = 0
# holds the number of duplicates found
mu = 0
# holds the number of missmatched UMI's

def position (x):                                   #x indicates a "" split line from a sam file
      bit_flag = int(x[4])                          #locates bit flag in a SAM file

      if bit_flag & 16 == 0:                        #identifies the posative strand
            a = re.findall('^[0-9]+S', x[5])        #identifies soft clipping and stores clipping value
            if a != []:
                  s = int(a[0][0:-1])
                  p = ("+",int(x[3]) - s)           #saves position with clipping value subtracted from the leftmost start position
            else:
                  p = ("+",x[3])                    #saves position as leftmost starting mopsition if no softclipping is detected and tags the strand as +

      else:                                         #identifies the negative strand
            a = re.findall('[0-9]+[DNMS]',x[5])     #finds all DNMS numbers in CIGAR string
            b = 0                                   #sets a variable to hold the value needed to add to the leftmost starting position
            if "S" in a[0]:                         #identifies if there is soft clipping at the beginning of the CIGAR string
                  for n,cig in enumerate(a):
                        if n !=0:
                              b += int(cig[0:-1])   #adds all the values in the CIGAR string that are not in the first S position
            else:
                  for cig in a:
                        b += int(cig[0:-1])         #adds all values that are in the CIGAR string
            p = ("-", int(x[3]) + b)                #adds the value gathered from the CIGAR strings to the leftmost start position and tags the strand as -
      return p                                      #returns the strand and position

def umi (x):
      h = x[0].split(":")[-1]
      return h                     #returns the UMI from the Qname


with open (u, "r") as uin:
    for line in (uin):
        line = line.strip()
        umis.append(line)       #saves all the UMI's to a list


with open (f, "r") as fin, open (o,"w") as oout:
    while True:                                                 #the while loop will allow me to close the loop wherever in the script I need to
        line = fin.readline().strip()
        if line == "":                                          #finds the end of the file
            for line in chrom_dict:
                print (chrom_dict[line],file = oout)            #writes out the lines saved in the chrom dictionary
            break                                               #ends the loop
        elif line[0] == "@":
            print (line,file = oout)                            #writes out the end file
        else:                                                   #this partof the loop will run through most of the file
            line_list = line.split()                            #this splits the line into a list at blank space
            p = position (line_list)                            #saves the position and strand of the line as p
            u = umi (line_list)                                 #saves the UMIof the line
            c = line_list[2]                                    #saves the chromosome number of the line 
            if dict_key != ():                                  #skips the first lap on each chromosome
                if c != dict_key[2]:                            #triggers at the change from one chromosome to the another
                    for line in chrom_dict:
                        print (chrom_dict[line],file = oout)    #prints out the lines saved for the dictionary
                    dict_key = ()                               #empties the dictionary key
                    chrom_dict = {}                             #empties the chromosome dictionary
            dict_key = (p[0],p[1],c,u)                          #populates the dictionary key with the strand, position, chromosome and UMI
            if u in umis:                                       #checks if the UMI is valid
                if dict_key in chrom_dict:                      #checks if the dictionary key has already been used
                    dup += 1                                    #enumerates the duplicate counter
                else:
                    chrom_dict[dict_key] = line                 #adds unique lines from the SAM to the chrom dictionary
            else:
                mu += 1                                         #enumerates the mismatched UMI's count
            


