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

TrigNumber = np.zeros(1,dtype=np.dtype("u4"))
Aaxis = np.zeros(1,dtype=np.dtype("u4"))
Baxis = np.zeros(1,dtype=np.dtype("u4"))
Caxis = np.zeros(1,dtype=np.dtype("u4"))
SampleNum = 0

f = TFile( outputfile, 'recreate' )
t = TTree( 'ttree', 'Test beam samples' )
t.Branch( 'Aaxis', Aaxis, 'Aaxis/i' )
t.Branch( 'Baxis', Baxis, 'Baxis/i' )
t.Branch( 'Caxis', Caxis, 'Caxis/i' )
t.Branch( 'TrigNumber', TrigNumber, 'TrigNumber/i' )

SamplesNumber = 25000
Voltage = np.zeros(SamplesNumber,dtype=np.Float_t)
Time = np.zeros(SamplesNumber,dtype=np.float32)
t.Branch( 'Time', Time, 'time[{}]/F'.format(SamplesNumber))
t.Branch( 'Voltage', Voltage, 'voltage[{}]/Float_t'.format(SamplesNumber))


for filename in os.listdir(in_folder):
    TrigNumber_tmp=1

    # Parsing the info of the input file's name
    if filename.endswith(".dat"):
        header, name, energy, HV, position, Aaxis_tmp, Baxis_tmp,Caxis_tmp = filename.split("_")

        # Split the headers and extract the useful info
        fh = open(filename)

        Caxis_tmp=Caxis_tmp.replace('.dat','')
        Caxis = int(Caxis_tmp)

        line = fh.readline()
        totTrig, scrap1, scrap2, scrap3 = line.split(" ")
        line = fh.readline()
        numberOfSamples, scrap4, scrap5, scrap6, scrap7, scrap8 = line.split(" ")
        line = fh.readline()
        TimeDiv, scrap9, scrap10, scrap11, scrap12 = line.split(" ")

        time = 0.0
        TimeDiv = float(TimeDiv)
        totTrig = int(totTrig)
        for i,tt in enumerate(np.linspace(0,(SamplesNumber-1)*TimeDiv,SamplesNumber)):
            Time[i] = tt

        if int(numberOfSamples) != SamplesNumber:
            print('Strange numberOfSamples: {}'.format(numberOfSamples))

        while True:
            # read line
            line = fh.readline()
            # check if line is not empty
            if not line:
                break

            #if the line doesen't contain text
            if not any(c.isalpha() for c in line):
                if SampleNum < Voltage.size:
                    Voltage[SampleNum] = float(line)
                    SampleNum += 1
            else:
                if SampleNum == SamplesNumber:
                    SampleNum = 0
                    t.Fill()
                    TrigNumber += 1
                elif SampleNum > 0:
                    print('Strange SampleNum: {}'.format(SampleNum))
                fh.readline()


        t.Fill()
        fh.close()

    else:
        continue


# Write in the output root file
f.Write()
f.Close()
