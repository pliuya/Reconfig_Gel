#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 12:09:09 2019
calculte stree and strain for each node
@author: yaliu
"""
import sys
if 'StressFunc' in sys.modules:  
    del sys.modules["StressFunc"]

import numpy as np
import pandas as pd
import StressFunc




datapath ='/Users/yaliu/Desktop/Project/Code/1-Reconfig-Gel/6-Stress-Strain/phi/'#'/Users/yaliu/Desktop/Project/Code/1-Reconfig-Gel/2-Ramp/'

inputinitname=  datapath + 'fullsample.gel.vtk'#'T34.late1000.gel.vtk'
inputgelname = datapath + '0000000155.gel.vtk'
inputstressname = datapath + '0000000155.stress.data'#'flatsample.gel.vtk'#'fullsample.gel.vtk'


outputstressname = datapath + 'stress.txt'
outputstress = open(outputstressname,'w+')
#outputstrain = open(outputstrainname,'w+')


inputinit = open(inputinitname,'r')
inputgel = open(inputgelname,'r')
inputstress = open(inputstressname,'r')

init()
#record initial configuration
#get index for top layer
print('0a: Start Reading fullsamplegel vtk')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
read_data(inputinit)
clear_array()
print('0b: End Reading fullsamplegel vtk')   

print('1a: Start Reading Current gel.vtk File')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
GelCount = 0
read_data(inputgel)

print('1b:End of ReadData')

print('2a: Start calculating strain')  
GelCount = 0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
while True:
    line = inputstress.readline()
    if line == '\n':
        continue
    if len(line) == 0:
        break
    fx = float(line.split()[0])
    fy = float(line.split()[1])
    fz = float(line.split()[2])
    tg = int(line.split()[3])
    fcomp = [fx,fy,fz]
    if GelCount in topgelinfo:
       stressinfo.append(fcomp) 
    GelCount += 1
print('End of Read Stress')

totinfo = np.hstack((gelcontour,stressinfo))
stress_analysis(totinfo)
for ti in range(MaxXLabel):
    outputstress.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(ti,profilez[ti],normalstress[ti],
                       normalvectorx[ti],normalvectory[ti]))

print('End of Calculating Stress')
