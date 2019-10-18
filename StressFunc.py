#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:56:22 2019

@author: yaliu
"""

def init():
    global GelType, GelCoord, CellData, PhiData, GelInitCoord
    global CellCount, GelCount, PhiCount
    global topgelinfo, topgelcoord, gelcontour, stressinfo
    global is_Coord, is_Cell, is_Cell_Type, is_Calc, is_Phi
    global PHIIni,phi0,LamdaIni,unit_space,layer, randdisp
    global MinXLabel, MaxXLabel
    global profilex, profilez,normalstress
    global normalvectorx, normalvectory
    global AnalysisMethod #=1: top gel layers; =2: skin layer\
    
    normalvectorx = []
    normalvectory = []
    topgelinfo = []
    topgelcoord = []
    gelcontour = []
    stressinfo = []
    profilex = []
    profilez = []
    normalstress =[]

    is_Coord = False
    is_Cell  = False
    is_Cell_Type = False
    is_Phi   = False
    is_Calc = False 

    GelType  = []
    GelCoord = []
    CellData = []
    PhiData  = []
    GelInitCoord = []
    CellCount = 0
    GelCount = 0 
    PhiCount = 0
    layer = 32
    PHIIni = 0.0894219
    phi0 = 0.1286
    LamdaIni = pow((phi0 / PHIIni), (1.0 / 3.0))
    unit_space = LamdaIni
    randdisp = 0.02
    MinXLabel = 1000
    MaxXLabel = 0
    
def clear_array():
    #topgelinfo = []
    #topgelcoord = []
    #gelcontour = []
    #stressinfo = []

    is_Coord = False
    is_Cell  = False
    is_Cell_Type = False
    is_Phi   = False
    is_Calc = False 

    GelType  = []
    GelCoord = []
    CellData = []
    PhiData  = []
    GelInitCoord = []
    CellCount = 0
    GelCount = 0 
    PhiCount = 0   

def label_func(vx):
    return int(vx/1.2)
#int(vx/unit_space+0.02)#int(vx/1.2) 
    #int((mono_x+randdisp*1.1)/unit_space)


def read_data(filename):
    global GelType, GelCoord, CellData, PhiData, GelInitCoord
    global CellCount, GelCount, PhiCount
    global topgelinfo, topgelinfo, gelcontour, stressinfo
    global is_Coord, is_Cell, is_Cell_Type, is_Calc, is_Phi
    global PHIIni,phi0,LamdaIni,unit_space,layer, randdisp
    global MinXLabel, MaxXLabel
    AnalysisMethod = 1
    while True:
        line = filename.readline()
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
            if 'fullsample' in filename.name:
                tx = float(line.split()[0])
                ty = float(line.split()[1])
                tz = float(line.split()[2])
                t  = float(line.split()[3])
                isCalNode = False
                if AnalysisMethod == 1:
                    zlayer = int((tz+randdisp*1.1)/unit_space)
                    if zlayer in range(layer-2,layer):
                        isCalNode = True   
                        
                xlayer = label_func(tx)
                if MinXLabel > xlayer:
                    MinXLabel = xlayer
                if MaxXLabel < xlayer:
                    MaxXLabel = xlayer
                    
                if isCalNode == True:
                    tempcord = [tx,ty,tz]
                    topgelinfo.append(GelCount)
                    topgelcoord.append(tempcord)
            else:
                if GelCount in topgelinfo:
                    mono_x = float(line.split()[0])
                    mono_y = float(line.split()[1])
                    mono_z = float(line.split()[2])
                    mono_coord = [mono_x,mono_y,mono_z]
                    gelcontour.append(mono_coord)
            GelCount += 1

def stress_analysis(rawdata):
    global profilex,profilez,normalstress,normalvectorx,normalvectory
    tanget = []
    count =[]
    ftot_x = []
    ftot_y = []
    ftot_z = []
    profilex = []
    profiley = []
    normalstress = []
    
    # Extract profile (xz plan)
    for i in range(MaxXLabel+1):
        profilex.append(0)
        profilez.append(0)
        count.append(0)
        ftot_x.append(0)
        ftot_y.append(0)
        ftot_z.append(0)
    
    for i in rawdata:
        mono_x = i[0]
        mono_y = i[1]
        mono_z = i[2]
        f_x = i[3]
        f_y = i[4]
        f_z = i[5]
        xlabel = label_func(mono_x)
        profilez[xlabel] +=  mono_z
        profilex[xlabel] +=  mono_x
        count[xlabel] += 1
        ftot_x[xlabel] += f_x
        ftot_y[xlabel] += f_y
        ftot_z[xlabel] += f_z
    
    for i in range(MaxXLabel+1):
        if count[i] != 0:
            profilex[i] /= count[i]
            profilez[i] /= count[i]
            ftot_x[i] /= count[i]
            ftot_y[i] /= count[i]
            ftot_z[i] /= count[i]
        else:
            print('No data point at {}\n'.format(i))

    for i in range(MaxXLabel+1):
        if count[i] == 0:
            if i == MaxXLabel:
                profilex[i] = (profilex[i-1]+profilex[0])*0.5
                profilez[i] = (profilez[i-1]+profilez[0])*0.5
                ftot_x[i]=(ftot_x[i-1]+ftot_x[0])*0.5
                ftot_y[i]=(ftot_y[i-1]+ftot_y[0])*0.5
                ftot_z[i]=(ftot_z[i-1]+ftot_z[0])*0.5
            if i == 0:
                profilex[i] = (profilex[MaxXLabel]+profilex[1])*0.5
                profilez[i] = (profilez[MaxXLabel]+profilez[1])*0.5
                ftot_x[i]=(ftot_x[MaxXLabel]+ftot_x[1])*0.5
                ftot_y[i]=(ftot_y[MaxXLabel]+ftot_y[1])*0.5
                ftot_z[i]=(ftot_z[MaxXLabel]+ftot_z[1])*0.5

    normalstress = []
    normalvectorx = []
    normalvectory = []
    normalstress.append(0)
    normalvectorx.append(0)
    normalvectory.append(0)
    
    for i in range(1,MaxXLabel):
        a = profilez[i]-profilez[i-1]
        ax = 1.2
        b = profilez[i+1]-profilez[i]
        bx = 1.2
        u1 =[ax/(ax**2+a**2)**0.5,a/(ax**2+a**2)**0.5]
        u2 =[bx/(bx**2+b**2)**0.5,b/(bx**2+b**2)**0.5]
        c = u1[0]+u2[0]
        d = u1[1]+u2[1]
        tangentvector = [c/(c**2+d**2)**0.5,d/(c**2+d**2)**0.5]
        tempnormalvector = [d/(c**2+d**2)**0.5,-c/(c**2+d**2)**0.5]
        normalstress.append(ftot_x[i]*tempnormalvector[0]+ftot_y[i]*tempnormalvector[1]) 
        normalvectorx.append(tangentvector[0])
        normalvectory.append(tangentvector[1])
        #normalvectorx.append(ftot_x[i]*tempnormalvector[0])#(ftot_x[i])#(ftot_x[i]*tempnormalvector[0])
        #normalvectory.append(ftot_y[i]*tempnormalvector[1])#(ftot_y[i])#(ftot_y[i]*tempnormalvector[1])
    normalstress.append(0)
    normalvectorx.append(0)
    normalvectory.append(0)
 