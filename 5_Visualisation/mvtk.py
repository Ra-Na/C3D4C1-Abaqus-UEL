# This python file generates VTK files that can be viewed with Paraview. 
# Both python ver. 2 and 3 work, but python 3 is much faster

# For a better quality output the field values are interpolted at the edge 
# centers such that a better cell model can be used in VTK/Paraview

# You have to specify the input files. 

inpfile="mesh0.inp"
inpfile2="IPVALS56_01.txt"
inpfile3="UNVALS56_01.txt"

inpfile="mesh3.inp"
inpfile2="IPVALS56_mesh3_01.txt"
inpfile3="UNVALS56_mesh3_01.txt"

outfile="Output1.vtk"

import string
from string import whitespace
import numpy as np
import sys
import time
start_time = time.time()

# initialize lists
nodes=list()
els=list()
intpunkte=list()

print("Read input files...")

# read node list from inp file
ifile=open(inpfile, "r")
read =str(ifile.readline())
while read[0:3] != '*No':
    read = str(ifile.readline())
read = str(ifile.readline())
while read[0]!="*":
    data=read.split(",")
    nodes.append([int(data[0]),float(data[1]),float(data[2]),float(data[3]),[0],[0],[0],[0],[0],[0]])
    knotenanzahl=int(data[0])
    read = str(ifile.readline())

# number of original nodes
unnodes=int(len(nodes))

# read element list - four corners, put placeholder for six edge centers 
while read[0:2] != '*E':
    read = str(ifile.readline())

read = str(ifile.readline())
while read[0]!="*":
    data=read.split(",")
    els.append([int(data[0]),int(data[1]),int(data[2]),int(data[3]),int(data[4]),0,0,0,0,0,0])
    elementnummer=int(data[0])
    read = str(ifile.readline())
ifile.close()

# Append displacmeents to node list
ifile3=open(inpfile3, "r")
for i in range(len(els)):
    read1 =str(ifile3.readline()).split()
    for j in range(int(read1[1])):
        read2 =str(ifile3.readline()).split()
        nodes[els[int(read1[0])-1][j+1]-1][4]=[0,(float(read2[0]))]
        nodes[els[int(read1[0])-1][j+1]-1][5]=[0,(float(read2[1]))]
        nodes[els[int(read1[0])-1][j+1]-1][6]=[0,(float(read2[2]))]

# Calculate the edge nodes, checking whether the node is already in the list before appending
# The order is prescribed by the VTK format CELLS 211 2321
# node 5:  between corner nodes 1 and 2
# node 6:  between corner nodes 2 and 3
# node 7:  between corner nodes 3 and 1
# node 8:  between corner nodes 1 and 4
# node 9:  between corner nodes 2 and 4
# node 10: between corner nodes 3 and 4
print(str(len(nodes))+" nodes, "+str(len(els))+" elements.")
print("Creating edge nodes...")

alreadycreated={};

for i in range(len(els)):
    print(str(float(i)/float(len(els))*100.0)+"%")
    sys.stdout.write("\033[F\033[K")
# node 5
    current=str(sorted([nodes[els[i][1]-1][0],nodes[els[i][2]-1][0]]))
    try:
        Nummer_Knoten_5=alreadycreated[current]
    except KeyError:
        aa5=[int(len(nodes)+1),float(0.5*(nodes[els[i][1]-1][1]+nodes[els[i][2]-1][1])),float(0.5*(nodes[els[i][1]-1][2]+nodes[els[i][2]-1][2])),float(0.5*(nodes[els[i][1]-1][3]+nodes[els[i][2]-1][3])),[0],[0],[0],[0],[0],[0]]
        nodes.append(aa5)
        Nummer_Knoten_5=int(len(nodes))
        alreadycreated[current]=int(len(nodes))

# node 6
    current=str(sorted([nodes[els[i][2]-1][0],nodes[els[i][3]-1][0]]))
    try:
        Nummer_Knoten_6=alreadycreated[current]
    except KeyError:
        aa6=[int(len(nodes)+1),float(0.5*(nodes[els[i][2]-1][1]+nodes[els[i][3]-1][1])),float(0.5*(nodes[els[i][2]-1][2]+nodes[els[i][3]-1][2])),float(0.5*(nodes[els[i][2]-1][3]+nodes[els[i][3]-1][3])),[0],[0],[0],[0],[0],[0]]
        nodes.append(aa6)
        Nummer_Knoten_6=int(len(nodes))
        alreadycreated[current]=int(len(nodes))

# node 7
    current=str(sorted([nodes[els[i][3]-1][0],nodes[els[i][1]-1][0]]))
    try:
        Nummer_Knoten_7=alreadycreated[current]
    except KeyError:
        aa7=[int(len(nodes)+1),float(0.5*(nodes[els[i][3]-1][1]+nodes[els[i][1]-1][1])),float(0.5*(nodes[els[i][3]-1][2]+nodes[els[i][1]-1][2])),float(0.5*(nodes[els[i][3]-1][3]+nodes[els[i][1]-1][3])),[0],[0],[0],[0],[0],[0]]
        nodes.append(aa7)
        Nummer_Knoten_7=int(len(nodes))
        alreadycreated[current]=int(len(nodes))

# node 8
    current=str(sorted([nodes[els[i][1]-1][0],nodes[els[i][4]-1][0]]))
    try:
        Nummer_Knoten_8=alreadycreated[current]
    except KeyError:
        aa8=[int(len(nodes)+1),float(0.5*(nodes[els[i][1]-1][1]+nodes[els[i][4]-1][1])),float(0.5*(nodes[els[i][1]-1][2]+nodes[els[i][4]-1][2])),float(0.5*(nodes[els[i][1]-1][3]+nodes[els[i][4]-1][3])),[0],[0],[0],[0],[0],[0]]
        nodes.append(aa8)
        Nummer_Knoten_8=int(len(nodes))
        alreadycreated[current]=int(len(nodes))

# node 9
    current=str(sorted([nodes[els[i][2]-1][0],nodes[els[i][4]-1][0]]))
    try:
        Nummer_Knoten_9=alreadycreated[current]
    except KeyError:
        aa9=[int(len(nodes)+1),float(0.5*(nodes[els[i][2]-1][1]+nodes[els[i][4]-1][1])),float(0.5*(nodes[els[i][2]-1][2]+nodes[els[i][4]-1][2])),float(0.5*(nodes[els[i][2]-1][3]+nodes[els[i][4]-1][3])),[0],[0],[0],[0],[0],[0]]
        nodes.append(aa9)
        Nummer_Knoten_9=int(len(nodes))
        alreadycreated[current]=int(len(nodes))

# node 10
    current=str(sorted([nodes[els[i][3]-1][0],nodes[els[i][4]-1][0]]))
    try:
        Nummer_Knoten_10=alreadycreated[current]
    except KeyError:
        aa10=[int(len(nodes)+1),float(0.5*(nodes[els[i][3]-1][1]+nodes[els[i][4]-1][1])),float(0.5*(nodes[els[i][3]-1][2]+nodes[els[i][4]-1][2])),float(0.5*(nodes[els[i][3]-1][3]+nodes[els[i][4]-1][3])),[0],[0],[0],[0],[0],[0]]
        nodes.append(aa10)
        Nummer_Knoten_10=int(len(nodes))
        alreadycreated[current]=int(len(nodes))

    els[i][5]=Nummer_Knoten_5
    els[i][6]=Nummer_Knoten_6
    els[i][7]=Nummer_Knoten_7
    els[i][8]=Nummer_Knoten_8
    els[i][9]=Nummer_Knoten_9
    els[i][10]=Nummer_Knoten_10

print("Reading integration point values...")
# read integration point values

ifile2=open(inpfile2, "r")

for j in range(len(els)):
    print(str(float(j)/float(len(els))*100.0)+"%")
    sys.stdout.write("\033[F\033[K")
    i = 1
    intpunkte.append([])
    Zeile = ifile2.readline()
    data = Zeile.split()
    elmnr = str(data[0])
    anzahlzeilenfuerelement = int(data[1])

    while i <= anzahlzeilenfuerelement:
        Datenzeile = ifile2.readline()
        data2 = Datenzeile.split()
        intpunkte[j].append([int(data[0]),float(data2[0]),float(data2[1]),float(data2[2]),float(data2[3]),float(data2[4]),float(data2[5]),float(data2[6]), float(data2[7]), float(data2[8]), i])
        i += 1

ifile2.close

# We extrapolate the 29 results values from the integration points to the corner and edge nodes. 
# Firstly, the 10 parameters of a cubic interpolation polynomial are adopted to the 29 values
# by square error minimization. 

# Set up and solve eq. sys.
# b, c, d for stress values (pass(1,2,3) in UEL) 
# e, f, g for displacements 

A = np.zeros((29, 10))
b = np.zeros((29, 1))
c = np.zeros((29, 1))
d = np.zeros((29, 1))
e = np.zeros((29, 1))
f = np.zeros((29, 1))
g = np.zeros((29, 1))

# List of monomials
intkoord=list()

# List of field values at 29 points
loes=list()
print("Extrapolating values at corners and edge centers from integration points...")
for i in range(len(els)):
    print(str(float(i)/float(len(els))*100.0)+"%")
    sys.stdout.write("\033[F\033[K")
    zeilen=29
    for j in range(zeilen):
        # monomials 1, x, y, z, xx, yy, zz, xy, xz, yz
        intkoord.append([1.0, intpunkte[i][j][1], intpunkte[i][j][2], intpunkte[i][j][3], intpunkte[i][j][1]*intpunkte[i][j][1], intpunkte[i][j][2]*intpunkte[i][j][2], intpunkte[i][j][3]*intpunkte[i][j][3], intpunkte[i][j][1]*intpunkte[i][j][2], intpunkte[i][j][1]*intpunkte[i][j][3], intpunkte[i][j][2]*intpunkte[i][j][3]])
        loes.append([intpunkte[i][j][4], intpunkte[i][j][5], intpunkte[i][j][6], intpunkte[i][j][7], intpunkte[i][j][8], intpunkte[i][j][9]])

    for k in range(zeilen):
        A[k] = np.array(intkoord[k])
        b[k] = np.array(loes[k][0])
        c[k] = np.array(loes[k][1])
        d[k] = np.array(loes[k][2])
        e[k] = np.array(loes[k][3])
        f[k] = np.array(loes[k][4])
        g[k] = np.array(loes[k][5])

    B=np.transpose(A)
    C=np.dot(B,A)
    b1=np.dot(B,b)
    c1=np.dot(B,c)
    d1=np.dot(B,d)
    e1=np.dot(B,e)
    f1=np.dot(B,f)
    g1=np.dot(B,g)

    x1 = np.linalg.solve(C, b1)
    x2 = np.linalg.solve(C, c1)
    x3 = np.linalg.solve(C, d1)
    x4 = np.linalg.solve(C, e1)
    x5 = np.linalg.solve(C, f1)
    x6 = np.linalg.solve(C, g1)
    
    for m in range(10):
        x_node = nodes[els[i][m+1]-1][1]
        y_node = nodes[els[i][m+1]-1][2]
        z_node = nodes[els[i][m+1]-1][3]

# append solutions
        nodes[els[i][m+1]-1][7].append(1*float(x1[0])+x_node*float(x1[1])+y_node*float(x1[2])+z_node*float(x1[3])+x_node*x_node*float(x1[4])+y_node*y_node*float(x1[5])+z_node*z_node*float(x1[6])+x_node*y_node*float(x1[7])+x_node*z_node*float(x1[8])+y_node*z_node*float(x1[9]))
        nodes[els[i][m+1]-1][8].append(1*float(x2[0])+x_node*float(x2[1])+y_node*float(x2[2])+z_node*float(x2[3])+x_node*x_node*float(x2[4])+y_node*y_node*float(x2[5])+z_node*z_node*float(x2[6])+x_node*y_node*float(x2[7])+x_node*z_node*float(x2[8])+y_node*z_node*float(x2[9]))
        nodes[els[i][m+1]-1][9].append(1*float(x3[0])+x_node*float(x3[1])+y_node*float(x3[2])+z_node*float(x3[3])+x_node*x_node*float(x3[4])+y_node*y_node*float(x3[5])+z_node*z_node*float(x3[6])+x_node*y_node*float(x3[7])+x_node*z_node*float(x3[8])+y_node*z_node*float(x3[9]))
        nodes[els[i][m+1]-1][4].append(1*float(x4[0])+x_node*float(x4[1])+y_node*float(x4[2])+z_node*float(x4[3])+x_node*x_node*float(x4[4])+y_node*y_node*float(x4[5])+z_node*z_node*float(x4[6])+x_node*y_node*float(x4[7])+x_node*z_node*float(x4[8])+y_node*z_node*float(x4[9]))
        nodes[els[i][m+1]-1][5].append(1*float(x5[0])+x_node*float(x5[1])+y_node*float(x5[2])+z_node*float(x5[3])+x_node*x_node*float(x5[4])+y_node*y_node*float(x5[5])+z_node*z_node*float(x5[6])+x_node*y_node*float(x5[7])+x_node*z_node*float(x5[8])+y_node*z_node*float(x5[9]))
        nodes[els[i][m+1]-1][6].append(1*float(x6[0])+x_node*float(x6[1])+y_node*float(x6[2])+z_node*float(x6[3])+x_node*x_node*float(x6[4])+y_node*y_node*float(x6[5])+z_node*z_node*float(x6[6])+x_node*y_node*float(x6[7])+x_node*z_node*float(x6[8])+y_node*z_node*float(x6[9]))


    del intkoord[:]
    del loes[:]

for i in range(len(nodes)):
    del nodes[i][4][0]
    del nodes[i][5][0]
    del nodes[i][6][0]
    del nodes[i][7][0]
    del nodes[i][8][0]
    del nodes[i][9][0]

# Since most nodes are shared by several elements, we have as many extrapolations
# as adjacent elements. Replace value by average.

print("Averaging results from adjacent elements...")
for n in range(len(nodes)):
    nodes[n][4]=np.mean(nodes[n][4])
    nodes[n][5]=np.mean(nodes[n][5])
    nodes[n][6]=np.mean(nodes[n][6])
    nodes[n][7]=np.mean(nodes[n][7])
    nodes[n][8]=np.mean(nodes[n][8])
    nodes[n][9]=np.mean(nodes[n][9])

print("Writing VTK file for Paraview...")

# Write results to vtk file
ofile=open(outfile, "w")
ofile.write("# vtk DataFile Version 2.0"+"\n"+"Titel der Datei"+"\n"+"ASCII"+"\n"+"DATASET UNSTRUCTURED_GRID"+"\n")
ofile.write("\n" + "POINTS "+str(len(nodes))+" float"+"\n")

for i in range(len(nodes)):
    ofile.write(    str(nodes[i][1])+"   " +     str(nodes[i][2])+ "   " +     str(nodes[i][3])+"\n")

ofile.write("\n" + "CELLS " + str(len(els)) + " " + str(len(els)*11) + "\n")

# VTK indexing starts at zero!
for i in range(len(els)):
    ofile.write("10 " + str(els[i][1]-1) + " " + str(els[i][2]-1) + " " + str(els[i][3]-1) + " " + str(els[i][4]-1) + " " + str(els[i][5]-1) + " " + str(els[i][6]-1) + " " + str(els[i][7]-1) + " " + str(els[i][8]-1) + " " + str(els[i][9]-1) + " " + str(els[i][10]-1) + "\n")

# all cells are of type 24
ofile.write("\n" + "Cell_Types " + str(len(els)) + "\n")
for i in range(len(els)):
    ofile.write("24" + "\n")

# finally write node values to file
ofile.write("\n" + "POINT_DATA " + str(len(nodes)) + "\n" + "VECTORS Verschiebungen float" + "\n")
for i in range(len(nodes)):
    ofile.write(str(nodes[i][4]) + " " + str(nodes[i][5]) + " " + str(nodes[i][6]) + "\n")

ofile.write("\n" + "VECTORS Spannungen float" + "\n")

for i in range(len(nodes)):
    ofile.write(str(nodes[i][7]) + " " + str(nodes[i][8]) + " " + str(nodes[i][9]) + "\n")

ofile.close()

print("Elasped time: " + str(time.time() - start_time))
