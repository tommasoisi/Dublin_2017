import time
#import matplotlib.pyplot as plt # https://matplotlib.org/
import numpy as np # http://www.numpy.org/
import os,sys
import ctypes
import optparse
import argparse

#from ROOT import TFile, TTree
from array import array
import re
import ROOT
import os
# from ROOT import TFile, TTree




parser = argparse.ArgumentParser(description='Run info.')

parser.add_argument('--path', metavar='path', type=str, help='Path to the input file',required=True)
args = parser.parse_args()

in_folder = str(args.path)

os.chdir(in_folder)

for filename in os.listdir(in_folder):
    TrigNumber_tmp=1
    if filename.endswith(".dat"):
        header, name, energy, HV, position, Aaxis_tmp, Baxis_tmp,Caxis_tmp = filename.split("_")

out_folder = "tree_produced"

outputfile = in_folder + '/%s/%s.root'%(out_folder,unName)


# # if len(sys.argv)==3:
# #     path = sys.argv[1]
# #     outputfile = sys.argv[2]
# #
# # if len(sys.argv)==2:
# #     path = sys.argv[1]
# #     runName = path.split('/')[-1].strip('/')
# #     if runName=='':
# #         runName = path.split('/')[-2].strip('/')
# #     # print(runName)
# #     print("Using default output: %tree_produceds"%(outputfile))
# #
# #
#
# if len(sys.argv)<2 or len(sys.argv)>3:
# raise Exception("Usage: python converter.py inputDir [outputFile]")

# f = TFile( outputfile, 'recreate' )
# t = TTree( 'pulse', 'Test beam samples' )


fh = open(filename)

Aaxis=[0]
Baxis=[0]
Caxis=[0]



Caxis_tmp=Caxis_tmp.replace('.dat','')
Aaxis.append(Aaxis_tmp)
# t.Branch( 'Aaxis', Aaxis, 'A' )
Baxis.append(Baxis_tmp)
Caxis.append(Caxis_tmp)


for i in range (1):
    line = fh.readline()
    totTrig, scrap1, scrap2, scrap3 = line.split(" ")

for i in range (1):
    line = fh.readline()
    sampleNum, scrap4, scrap5, scrap6, scrap7, scrap8 = line.split(" ")


for i in range (1):
    line = fh.readline()
    TimeDiv, scrap9, scrap10, scrap11, scrap12 = line.split(" ")

time = 0.0

TimeDiv = float(TimeDiv)
sampleNum = int(sampleNum)
totTrig = int(totTrig)
Voltage = []
TrigNumber_tmp = 0.0
n=0

Time=[]

TrigNumber=[]

while True:
    # read line

    line = fh.readline()
    trigUpdate = False
    my_vector = ROOT.vector('int')()

    # check if line is not empty
    if not line:
        break

    if not any(c.isalpha() for c in line):
        Voltage.append(line)
        time = time + TimeDiv
        Time.append(time)
        TrigNumber.append(TrigNumber_tmp)
        my_vector.push_back(TrigNumber_tmp)

    else:
        for j in range (1):
            next(fh)
        if trigUpdate == False:
            TrigNumber_tmp = TrigNumber_tmp + 1
            time = 0.0
            trigUpdate = True



Voltage = [float(i) for i in Voltage]
TrigNumber = [int(i) for i in TrigNumber]
Time = [float(i) for i in Time]
Aaxis_tmp = [int(i) for i in Aaxis]
Baxis_tmp = [int(i) for i in Baxis]
Caxis_tmp = [int(i) for i in Caxis]



print(len(TrigNumber))
print(len(Voltage))
print(len(Time))


print(Time)

# t.Branch( 'TrigNumber', TrigNumber 'TrigNumber/i' )
# t.Branch( 'channel', channel, 'channel[4][1000]/F' )
# t.Branch( 'time', time, 'time[1][1000]/F')


# 25000 Number of samples per trigger
# 2.5000000000e-010 Time division in seconds
fh.close()



# print(TrigNumber_tmp)
