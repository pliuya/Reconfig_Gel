#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 12:09:09 2019
# read vtk file and extract information given coordinates
@author: yaliu
"""
import numpy as np
import pandas as pd

writeconfig = True


datapath ='/Users/yaliu/Desktop/Project/Code/1-Reconfig-Gel/2-Ramp/Fast36-2/'
inputname=  datapath + '0000000089.gel.vtk'
outputname =  datapath + 'result.txt'


spotheight = 29

inputfile = open(inputname,'r')
outputfile = open(outputname,'w+')

GelCoord = []
CellData = []
PhiData  = []


is_Coord = False
is_Cell  = False
is_Cell_Type = False
is_Phi   = False
is_Calc = False  
      
GelCount = 0
CellCount = 0
PhiCount = 0

print('I: Start Reading File')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
while True:
    line = inputfile.readline()
    if line == '\n':
        continue
    if len(line) == 0:
        print('a1: No input file/empty input file/finish input')
        break
    if line.split()[0] == 'POINTS' and  line.split()[1] == '35552':
        print('II:Start Coordinate')
        is_Coord = True
        is_Cell  = False
        is_Cell_Type = False
        is_Phi   = False
        continue
    if line.split()[0] == 'CELLS' and line.split()[1] == '31000':
        print('III:Start Cell Data')
        is_Coord = False
        is_Cell  = True
        is_Cell_Type = False
        is_Phi   = False
        continue
    if line.split()[0] == 'CELL_TYPES':
        is_Coord = False
        is_Cell  = False
        is_Cell_Type = True
        is_Phi   = False
        continue
    if line.split()[0] == 'LOOKUP_TABLE' and  line.split()[1] == 'default':
        print('IV:Start Phi')
        is_Coord = False
        is_Cell  = False
        is_Cell_Type = False
        is_Phi   = True
        continue  
    if is_Coord == True:
        mono_x = float(line.split()[0])
        mono_y = float(line.split()[1])
        mono_z = float(line.split()[2])
        mono_coord = [GelCount,mono_x,mono_y,mono_z]
        if mono_z>=spotheight and mono_z<spotheight+1:
            GelCoord.append(mono_coord)
        GelCount += 1
    if is_Cell  == True:
        label1 = int(line.split()[1])
        label2 = int(line.split()[2])
        label3 = int(line.split()[3])
        label4 = int(line.split()[4])
        label5 = int(line.split()[5])
        label6 = int(line.split()[6])
        label7 = int(line.split()[7])
        label8 = int(line.split()[8])
        cellinfo = [CellCount,label1,label2,label3,label4,label5,label6,label7,label8]
        CellData.append(cellinfo) 
        CellCount += 1
    if is_Phi == True:
        phivalue = float(line.split()[0])
        PhiData.append(phivalue)
      
#extract phi for given gel height
"""
for tempx in GelCoord:
    label = tempx[0]
    for tempcell in CellData:
        if label in tempcell[1:]:
            philabel = tempcell[0]
            #outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
            if tempx[2] < 0:
               outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))

for tempx in GelCoord:
    label = tempx[0]
    for tempcell in CellData:
        if label in tempcell[1:]:
            philabel = tempcell[0]
            #outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
            if tempx[2] > 11:
               outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))

for tempx in GelCoord:
    label = tempx[0]
    for tempcell in CellData:
        if label in tempcell[1:]:
            philabel = tempcell[0]
            #outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
            if tempx[2] < 6.5 and tempx[2] > 5.5:
               outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))


for tempx in GelCoord:
    label = tempx[0]
    for tempcell in CellData:
        if label in tempcell[1:]:
            philabel = tempcell[0]
            #outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
            if tempx[2] < 9.5 and tempx[2] > 9:
               outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))

for tempx in GelCoord:
    label = tempx[0]
    for tempcell in CellData:
        if label in tempcell[1:]:
            philabel = tempcell[0]
            #outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
            if tempx[2] > 2 and tempx[2] < 2.8:
               outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
"""
for tempx in GelCoord:
    label = tempx[0]
    for tempcell in CellData:
        if label in tempcell[1:]:
            philabel = tempcell[0]
            #outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
            if tempx[2] < 6.5 and tempx[2] > 5:
               outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))

           
print('End of ReadData')