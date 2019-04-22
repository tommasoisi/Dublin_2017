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

# Declare Input folder
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
        
 

# my_vector = ROOT.vector('int')()
# my_vector = ROOT.vector('int')(2)


Aaxis=[]
Baxis=[]
Caxis=[]
# my_vector = ROOT.vector('int')()


f = TFile( outputfile, 'recreate' )
t = TTree( 'pulse', 'Test beam samples' )

maxn = 10
n = array( 'i', [ filenumber ] )
d = array( 'f', maxn*[ 0. ] )
# t.Branch( 'Caxis', Caxis, 'Caxis/I')
# t.Branch( 'mynum', n, 'mynum/I' )
t.Branch( 'myval', d, 'myval[mynum]/F' )

for filename in os.listdir(in_folder):
    TrigNumber_tmp=1

        # Parsing the info of the input file's name
    if filename.endswith(".dat"):
        header, name, energy, HV, position, Aaxis_tmp, Baxis_tmp,Caxis_tmp = filename.split("_")

        fh = open(filename)


        Caxis_tmp=Caxis_tmp.replace('.dat','')
        Aaxis.append(Aaxis_tmp)
        # t.Branch( 'Aaxis', Aaxis, 'A' )
        Baxis.append(Baxis_tmp)
        Caxis.append(Caxis_tmp)

        Caxis_tmp = int(Caxis_tmp)


        line = fh.readline()
        totTrig, scrap1, scrap2, scrap3 = line.split(" ")
        line = fh.readline()
        sampleNum, scrap4, scrap5, scrap6, scrap7, scrap8 = line.split(" ")
        line = fh.readline()
        TimeDiv, scrap9, scrap10, scrap11, scrap12 = line.split(" ")

        time = 0.0
        TimeDiv = float(TimeDiv)
        sampleNum = int(sampleNum)
        totTrig = int(totTrig)
        Voltage = []
        TrigNumber_tmp = 0.0

        # my_vector.push_back(Caxis_tmp)
        Time=[]
        print(Caxis_tmp)
        TrigNumber=[]

        while True:
            # read line
            line = fh.readline()
            trigUpdate = False
            # check if line is not empty
            if not line:
                break

            if not any(c.isalpha() for c in line):
                Voltage.append(line)
                time = time + TimeDiv
                Time.append(time)
                TrigNumber.append(TrigNumber_tmp)
                # my_vector.push_back(Caxis_tmp)
                Caxis.append(Caxis_tmp)



            else:
                for j in range (1):
                    fh.readline()
                # if trigUpdate == False:
                #     TrigNumber_tmp = TrigNumber_tmp + 1
                #     time = 0.0
                #     trigUpdate = True


        # print(len(my_vector))
        print(len(TrigNumber))
        print(len(Voltage))
        print(len(Time))

        t.Fill()

        # print(Time)

        # t.Branch( 'channel', channel, 'channel[4][1000]/F' )
        # t.Branch( 'time', time, 'time[1][1000]/F')

        # 25000 Number of samples per trigger
        # 2.5000000000e-010 Time division in seconds
        fh.close()

    else:
        continue




Voltage = [float(i) for i in Voltage]
TrigNumber = [int(i) for i in TrigNumber]
Time = [float(i) for i in Time]
Aaxis_tmp = [int(i) for i in Aaxis]
Baxis_tmp = [int(i) for i in Baxis]
Caxis = [int(i) for i in Caxis]
t.Branch( 'mynum', n, 'mynum/I' )

# print(Caxis)

# print(my_vector)

#vec = np.asarray(Caxis)
#ciccio = array2root(vec, name ='ciccio')
#array2root(vec,'test.root','ciccio')

#
# pulse = array2tree(array, name='pulse')
# array2root(array,'test3.root','pulse')
#

#
# while i < len(Caxis):
#     my_vector.push_back(Caxis[i])
#     i = i+1
# Caxis_vec = [int(i) for i in Caxis]
# my_vector = ROOT.vector('int')(len(Caxis))

# for i in Caxis:
#     my_vector.push_back(Caxis[i])

# while i < len(Caxis) :
#     my_vector.push_back(Caxis[i])
#     i += 1


# my_vector = vec
# print(my_vector)



# for i in my_vector:

f.Write()
f.Close()

# print(TrigNumber_tmp)
