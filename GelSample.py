import argparse
import random
from datetime import datetime

#parser = argparse.ArgumentParser()
#parser.add_argument("-p", "--perturbation", action="store", type=float, help="perturbation in percentage of unit_space")
#parser.add_argument("-f", "--foldname", action="store", type=str, help="fold where output vtk file")
#args = parser.parse_args()
perturbation = 0.02
shift = 0# total is 80*period-2*period*shift=128 
includebump = True 

foldname = '/Users/yaliu/Desktop/Project/Code/1-Reconfig-Gel'
fibers = []
filename = foldname+"/thickconf"
inputfile=open(filename,'r')
while True:
    line=inputfile.readline()
    if len(line)==0:
        print('no input file')
        break
    if line.split()[0] == '#py#' and line.split()[1] == 'N':
        Lx = int(line.split()[2])
        Ly = int(line.split()[3])
        Lz = int(line.split()[4])
        break
"""
    if line.split()[0] == '#py#' and line.split()[1] == 'fiber':
        one_fiber = []
        for i in range(6):
            one_fiber.append(int(line.split()[i+2]))
        fibers.append(one_fiber)
"""
inputfile.close()
#Ya: bump geometrical parameter
bump_top_thick = 10
period = 1 # 2: two bumps
bump_height = Lz - bump_top_thick
origin_Lx = Lx/period+2*shift
bump_radius = ((0.5*origin_Lx)**2 + bump_height**2)/(2*bump_height)
bump_centerL = [-1*shift,bump_height-bump_radius]
bump_centerR = [Lx+shift,bump_height-bump_radius]
bool_bump = False

deltaX = 1.0
#phi0 = 0.1286
#PHIIni=0.0428939 (t=20),PHIIni=0.0445666 (t=21),PHIIni=0.0464266(t=22),PHIIni=0.0485128 (t=23)
#PHIIni=0.0508769 (t=24),PHIIni=0.0535892 (t=25),PHIIni=0.0567486 (t=26),PHIIni=0.0605002 (t=27)
#PHIIni=0.0650679 (t=28),PHIIni=0.0708237 (t=29),PHIIni=0.0784509 (t=30),PHIIni=0.0894219 (t=31)
#PHIIni=0.108011 (t=32),PHIIni=0.167191 (t=33),PHIIni=0.352235 (t=34),PHIIni=0.41007 (t=35)
#PHIIni=0.448942 (t=36)

PHIIni = 0.0894219
phi0 = 0.1286
LamdaIni = pow((phi0 / PHIIni), (1.0 / 3.0))
unit_space = deltaX * LamdaIni


unit_space_Bump = 1

EPS = 1e-6

random.seed(datetime.now())

# set up nodes positions
# nodes_ijk [i,j,k,0/1]: 0:bump 1:gel 
nodes_position=[]
nodes_ijk=[]
bool_ijk=[]
num_bump = 0

periodsize = Lx/period
for i in range(Lx):
    for j in range(Ly):
        for k in range(Lz):
            x = (i + (2.0 * random.random() - 1.0) * perturbation) * unit_space
            y = (j + (2.0 * random.random() - 1.0) * perturbation) * unit_space
            z = (k + (2.0 * random.random() - 1.0) * perturbation) * unit_space
            #bool_bump = 0: bump; = 1: boundray; =2: gel
            el = 0
            er = 0
            if i<0.5*periodsize:
                if includebump == False:
                    el = 0
                else:
                    el = (bump_radius**2-(i+shift)**2)**0.5-(bump_radius-bump_height)
            elif i>=periodsize and i<1.5*periodsize:
                if includebump == False:
                    el2 = 0
                else:
                    el2 = (bump_radius**2-(i+shift-periodsize)**2)**0.5-(bump_radius-bump_height)
            elif i>=0.5*periodsize and i<periodsize:
                if includebump == False:
                    er = 0
                else:
                    er = (bump_radius**2-(periodsize+shift-1-i)**2)**0.5-(bump_radius-bump_height)
            elif i>=1.5*periodsize:
                if includebump == False:
                    er2 = 0
                else:
                    er2 = (bump_radius**2-(periodsize*2+shift-1-i)**2)**0.5-(bump_radius-bump_height)
            if i<0.5*periodsize and k < el:
                bool_bump = 0
            elif i>=periodsize and i<1.5*periodsize  and k < el2:
                bool_bump = 0
            elif i>=0.5*periodsize and i<periodsize and k <er :
                bool_bump = 0
            elif i>=1.5*periodsize and k <er2 :
                bool_bump = 0
            elif (i==0 or i==Lx-1) and (j==0 or j==Ly-1):
                bool_bump = 5
            elif i==0 or i==Lx-1:
                bool_bump = 4
            elif j==0 or j==Ly-1:
                bool_bump = 1
            elif k == Lz-1:
                bool_bump = 2
            else:
                bool_bump = 3    
            nodes_position.append([x,y,z,bool_bump])
            nodes_ijk.append([i,j,k])
            bool_ijk.append([x,y,z,bool_bump])

i_max = 0
j_max = 0
k_max = 0
for node in nodes_ijk:
    i = node[0]
    j = node[1]
    k = node[2]
    if i > i_max:
        i_max = i
    if j > j_max:
        j_max = j
    if k > k_max:
        k_max = k

# record if position i,j,k is occupied by a node, if yes, store the node index+1
occupation_ijk = [[[0 for k in range(k_max+1)] for j in range(j_max+1)] for i in range(i_max+1)]

for l in range(len(nodes_ijk)):
    occupation_ijk[nodes_ijk[l][0]][nodes_ijk[l][1]][nodes_ijk[l][2]] = l + 1

# find elements
elements=[]
for loop in nodes_ijk:
    i = loop[0]
    j = loop[1]
    k = loop[2]
    if i < i_max and j < j_max and k < k_max:
        node_1 = occupation_ijk[i][j][k] - 1
        node_2 = occupation_ijk[i][j][k+1] - 1
        node_3 = occupation_ijk[i+1][j][k+1] - 1
        node_4 = occupation_ijk[i+1][j][k] - 1
        node_5 = occupation_ijk[i][j+1][k] - 1
        node_6 = occupation_ijk[i][j+1][k+1] - 1
        node_7 = occupation_ijk[i+1][j+1][k+1] - 1
        node_8 = occupation_ijk[i+1][j+1][k] - 1
        if node_2 >= 0 and node_3 >= 0 \
            and node_4 >= 0 and node_5 >= 0 \
            and node_6 >= 0 and node_7 >= 0 \
            and node_8 >= 0:
                elements.append([node_1,node_2,node_3,node_4,node_5,node_6,node_7,node_8])


# output sample.gel.vtk file
filename = foldname+"/thicksample.gel.vtk"

outputfile=open(filename,'w+')

outputfile.write("# vtk DataFile Version 2.0\n")
outputfile.write("Unstructured Grid\n")
outputfile.write("ASCII\n")
outputfile.write("DATASET UNSTRUCTURED_GRID\n")
outputfile.write("POINTS %d double\n"%(len(nodes_position)))
for loop in nodes_position:
    outputfile.write("%12.5e %12.5e %12.5e\n"%(loop[0],loop[1],loop[2]))

outputfile.write("\n")
outputfile.write("CELLS %d %d\n"%(len(elements), 9*len(elements)))
for loop in elements:
    outputfile.write("%-6d"%8)
    for i in range(8):
        outputfile.write(" %-6d"%loop[i])
    outputfile.write("\n")
outputfile.write("\n")

outputfile.write("CELL_TYPES %d\n"%(len(elements)))
for loop in elements:
    outputfile.write("12\n")
outputfile.write("\n")

outputfile.write("CELL_DATA %d\n"%(len(elements)))
outputfile.write("SCALARS cellData double 1\n")
outputfile.write("LOOKUP_TABLE default\n")
for loop in elements:
    if nodes_position[loop[0]][3]==0 and nodes_position[loop[1]][3]==0 and nodes_position[loop[2]][3]==0\
         and nodes_position[loop[3]][3]==0 and nodes_position[loop[4]][3]==0 and nodes_position[loop[5]][3]==0\
         and nodes_position[loop[6]][3]==0 and nodes_position[loop[7]][3]==0:
             outputfile.write("%12.5e\n"%1)
    else:
        outputfile.write("%12.5e\n"%PHIIni)
    
      
             
outputfile.write("\n")

outputfile.close()

# output fullsample.gel.vtk file
fullfilename = foldname+"/fullthicksample.gel.vtk"

outputfullfile=open(fullfilename,'w+')

outputfullfile.write("# vtk DataFile Version 2.0\n")
outputfullfile.write("Unstructured Grid\n")
outputfullfile.write("ASCII\n")
outputfullfile.write("DATASET UNSTRUCTURED_GRID\n")
outputfullfile.write("POINTS %d double\n"%(len(nodes_position)))
for loop in nodes_position:
    outputfullfile.write("%12.5e %12.5e %12.5e %d\n"%(loop[0],loop[1],loop[2],loop[3]))

outputfullfile.write("\n")
outputfullfile.write("CELLS %d %d\n"%(len(elements), 9*len(elements)))
for loop in elements:
    outputfullfile.write("%-6d"%8)
    for i in range(8):
        outputfullfile.write(" %-6d"%loop[i])
    outputfullfile.write("\n")
outputfullfile.write("\n")

outputfullfile.write("CELL_TYPES %d\n"%(len(elements)))
for loop in elements:
    outputfullfile.write("12\n")
outputfullfile.write("\n")

outputfullfile.write("CELL_DATA %d\n"%(len(elements)))
outputfullfile.write("SCALARS cellData double 1\n")
outputfullfile.write("LOOKUP_TABLE default\n")
for loop in elements:
    outputfullfile.write("%12.5e\n"%PHIIni)
outputfullfile.write("\n")

outputfullfile.close()
"""
# output sample.fibers.vtk
outputfile=foldname+"/sample.fibers.vtk"
outputfile=open(filename,'w+')
outputfile.write("# vtk DataFile Version 2.0\n")
outputfile.write("polydata\n")
outputfile.write("ASCII\n")
outputfile.write("DATASET POLYDATA\n")
outputfile.write("POINTS %d double\n"%(len(nodes_position)))
for loop in nodes_position:
        outputfile.write("%12.5e %12.5e %12.5e\n"%(loop[0],loop[1],loop[2]))

outputfile.write("\n")
num_nodes_on_fiber = 0
fibers_length = []
for loop in fibers:
    count = 0
    for i in range(loop[0], loop[3] + 1):
        for j in range(loop[1], loop[4] + 1):
            for k in range(loop[2], loop[5] + 1):
                num_nodes_on_fiber = num_nodes_on_fiber + 1
                count = count + 1
    fibers_length.append(count)
outputfile.write("LINES %d %d\n"%(len(fibers), len(fibers)+num_nodes_on_fiber))
for l in range(len(fibers)):
    loop = fibers[l]
    outputfile.write("%-6d"%fibers_length[l])
    for i in range(loop[0], loop[3] + 1):
        for j in range(loop[1], loop[4] + 1):
            for k in range(loop[2], loop[5] + 1):
                outputfile.write(" %-6d"%(occupation_ijk[i][j][k] - 1))
    outputfile.write("\n")
outputfile.write("\n")

outputfile.close()
"""
'''
#output position

outputfile=file(args.foldname+"/sample.bump.vtk",'w+')
outputfile.write("# vtk DataFile Version 2.0\n")
outputfile.write("Unstructured Grid\n")
outputfile.write("ASCII\n")
outputfile.write("DATASET UNSTRUCTURED_GRID\n")
outputfile.write("POINTS %d double\n"%(num_bump))
for loop in nodes_position:
    if bool_ijk[loop[0]][loop[0]][loop[2]][3] == 0:
        outputfile.write("%12.5e %12.5e %12.5e\n"%(loop[0],loop[1],loop[2]))

outputfile.write("\n")
#count number of cell 



outputfile.write("CELLS %d %d\n"%(len(elements), 9*len(elements)))
for loop in elements:
    outputfile.write("%-6d"%8)
    for i in range(8):
        outputfile.write(" %-6d"%loop[i])
    outputfile.write("\n")
outputfile.write("\n")

outputfile.write("CELL_TYPES %d\n"%(len(elements)))
for loop in elements:
    outputfile.write("12\n")
outputfile.write("\n")

outputfile.write("CELL_DATA %d\n"%(len(elements)))
outputfile.write("SCALARS cellData double 1\n")
outputfile.write("LOOKUP_TABLE default\n")
for loop in elements:
    outputfile.write("%12.5e\n"%PHIIni)
outputfile.write("\n")

outputfile.close()
'''
fileinfoname = foldname+"/info.txt"
infooutputfile=open(fileinfoname,'w+')
for loop in bool_ijk:
    infooutputfile.write("%12.5e %12.5e %12.5e %d\n"%(loop[0],loop[1],loop[2],loop[3]))
infooutputfile.close()
# for loop in elements:
#     print loop
