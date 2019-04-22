##########################################################################
##                      DUBLIN CONVERTER SCRIPT                         ##
##                                                                      ##
##      SCRIPT TO CONVERT THE OSCILLOSCOPE .dat FILES IN TTRE           ##
##   Works converting and merging all the file inside the input folder  ##
##########################################################################

import time
import glob
#import matplotlib.pyplot as plt # https://matplotlib.org/
import ROOT
import numpy as np # http://www.numpy.org/
import os,sys
import ctypes
import optparse
import argparse

from ROOT import TFile, TTree
from array import array
from root_numpy import array2tree
from root_numpy import array2root
import re
import os
# from ROOT import TFile, TTree



# The code requires the complete input folder directory to run
parser = argparse.ArgumentParser(description='Run info.')
parser.add_argument('--path', metavar='path', type=str, help='Path to the input file',required=True)
args = parser.parse_args()

# Declare Input folder and list the filenumber 
in_folder = str(args.path)
filenumber=(len(glob.glob1(in_folder,"*.dat")))

# Open input folder
os.chdir(in_folder)

# Declare the output file name and path
for filename in os.listdir(in_folder):
        # Parsing the info of the input file's name
    if filename.endswith(".dat"):
        header, name, energy, HV, position, Aaxis_tmp, Baxis_tmp,Caxis_tmp = filename.split("_")

    out_folder = "tree_produced"
    out_name = name + ("_") + energy + ("_") + HV + ("_") + ("_v1")
    outputfile = in_folder + '/%s/%s.root'%(out_folder,out_name)

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

Caxis = array('i')
# Caxis_array = array( 'i', [filenumber] )
f = TFile( outputfile, 'recreate' )
t = TTree( 'pulse', 'Test beam samples' )


l = 0 

for filename in os.listdir(in_folder):
    TrigNumber_tmp=1

    # Parsing the info of the input file's name
    if filename.endswith(".dat"):
        header, name, energy, HV, position, Aaxis_tmp, Baxis_tmp,Caxis_tmp = filename.split("_")
        
        # Split the headers and extract the useful info
        fh = open(filename)
        line = fh.readline()
        totTrig, scrap1, scrap2, scrap3 = line.split(" ")
        line = fh.readline()
        sampleNum, scrap4, scrap5, scrap6, scrap7, scrap8 = line.split(" ")
        line = fh.readline()
        TimeDiv, scrap9, scrap10, scrap11, scrap12 = line.split(" ")
        continue 
                
TotData = int(filenumber)*int(sampleNum)*int(totTrig)
print(TotData)

TrigNumber = [] 
t.Branch( 'Caxis', Caxis, 'Caxis/i' )


for filename in os.listdir(in_folder):   
    if filename.endswith(".dat"):
        header, name, energy, HV, position, Aaxis_tmp, Baxis_tmp,Caxis_tmp = filename.split("_")
        fh = open(filename)
        Caxis_tmp=Caxis_tmp.replace('.dat','')
        Caxis_tmp = int(Caxis_tmp)        
        
        # Split the headers and extract the usual info
        line = fh.readline()
        line = fh.readline()
        line = fh.readline()
               
        time = 0.0
        TimeDiv = float(TimeDiv)
        sampleNum = int(sampleNum)
        totTrig = int(totTrig)
        TrigNumber_tmp = 0.0
        Time=[]
        
        while True:
            # read line
            line = fh.readline()
            # check if line is not empty
            if not line:
                break

            if not any(c.isalpha() for c in line):
#                 Voltage.append(line)
                Caxis.append(Caxis_tmp)
                time = time + TimeDiv
                Time.append(time)
                TrigNumber.append(TrigNumber_tmp)
                # my_vector.push_back(Caxis_tmp)
#                 my_vector.append(Caxis_tmp)
#                 print(l)
                l = l + 1 
                t.Fill()



            else:
                for j in range (1):
                    fh.readline()

        fh.close()

    else:
        continue




# Caxis_array = np.asarray(Caxis)
# Voltage = [float(i) for i in Voltage]
# TrigNumber = [int(i) for i in TrigNumber]
# Time = [float(i) for i in Time]
# Aaxis_tmp = [int(i) for i in Aaxis]
# Baxis_tmp = [int(i) for i in Baxis]
# Caxis_array= [int(i) for i in Caxis]


# my_vector = vec
# print(my_vector)
print(Caxis)

# for i in my_vector:
f.Write()
f.Close()

# print(TrigNumber_tmp)
