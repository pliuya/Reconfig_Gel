#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:37:06 2019

@author: yaliu
"""

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


datapath ='/Users/yaliu/Desktop/Project/Code/1-Reconfig-Gel/7-ThickGel/Shink/'#'/Users/yaliu/Desktop/Project/Code/1-Reconfig-Gel/2-Ramp/'
inputname=  datapath + '0000002000.gel.vtk'#'T34.late1000.gel.vtk'

outputname =  datapath + 'result.txt'
outputphiname = datapath + 'phi.txt'

configname = datapath + 'late.gel.vtk'

outputconfig = open(configname,'w+')
outputphi = open(outputphiname,'w+')
outputfile = open(outputname,'w+')

spotheight = 29

inputfile = open(inputname,'r')
inputfullfile = open(inputfullname,'r')


GelType  = []
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

minphi = 1
maxphi = 0        
print('I: Start Reading File')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
while True:
    line = inputfile.readline()
    if line == '\n':
        outputconfig.write(line)
        continue
    if len(line) == 0:
        print('a1: No input file/empty input file/finish input')
        break
    if line.split()[0] == 'LOOKUP_TABLE':
        print('IV:Start Phi')
        outputconfig.write(line)
        is_Phi   = True
        continue  
    if is_Phi == False:
        outputconfig.write(line)
    if is_Phi == True:
        phivalue = float(line.split()[0])
        if phivalue == 1.:
            outputconfig.write('1.0\t0\n') 
        else:
            outputconfig.write('{0:f}\t{1:d}\n'.format(phivalue,1))
            outputphi.writelines('{0:f}\n'.format(phivalue))
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

for tempx in GelCoord:
    label = tempx[0]
    for tempcell in CellData:
        if label in tempcell[1:]:
            philabel = tempcell[0]
            #outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
            if tempx[2] < 6.5 and tempx[2] > 5:
               outputfile.write('{0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\n'.format(tempx[1],tempx[2],tempx[3],PhiData[philabel]))
"""         
print('Max Phi is: {0:f}'.format(maxphi))
print('Min Phi is: {0:f}'.format(minphi))
print('End of ReadData')