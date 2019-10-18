#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 10:34:11 2019
#calculate Hydrostatic stress and generate vtk file
@author: yaliu
"""
import numpy as np
import pandas as pd

datapath ='/Users/yaliu/Desktop/Project/Code/1-Reconfig-Gel/6-Stress-Strain/New-1/'#'/Users/yaliu/Desktop/Project/Code/1-Reconfig-Gel/2-Ramp/'

#input file inform
filelabel = '0000000100'
inputnameX = datapath + filelabel + '.x.vtk'
inputnameY = datapath + filelabel + '.y.vtk'
inputnameZ = datapath + filelabel + '.z.vtk'

inputX = open(inputnameX,'r')
inputY = open(inputnameY,'r')
inputZ = open(inputnameZ,'r')

xStress = []
yStress = []
zStress = []

#output file inform
outputname = datapath + filelabel + '.stress.vtk'
outputstress = open(outputname,'w+')

for tempfile in [inputX,inputY,inputZ]:
    isData = False
    while(True):
        line = tempfile.readline()
        if line == '\n':
            continue
        if len(line) == 0:
            break
        if line.split()[0] == 'LOOKUP_TABLE':
            isData = True
            continue
        if isData == True:
            if tempfile == inputX:
                xStress.append(float(line.split()[0]))
            elif tempfile == inputY:
                yStress.append(float(line.split()[0]))
            else:
                zStress.append(float(line.split()[0]))
#output stress
inputNew = open(inputnameZ,'r')       
isData = False     
while(True):
        linenew = inputNew.readline()
        if linenew == '\n':
            continue
        if len(linenew) == 0:
            break            
        if linenew.split()[0] == 'LOOKUP_TABLE':
            isData = True
            outputstress.write(linenew)
        if isData == False:
            outputstress.write(linenew)
            continue
        if isData == True:
            for i in range(len(xStress)):
                a =  (xStress[i] + yStress[i] + zStress[i])/3.
                outputstress.write('{}\n'.format(a))
        break
        
        
        
        
        
        
        
        