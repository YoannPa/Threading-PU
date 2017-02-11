#!/usr/bin/python
# -*- coding: utf-8 -*-

##Threading PU - Double Dynamic Programing

__author__ = "Yoann PAGEAUD"
__version__  = "1.0.0"
__copyright__ = "copyleft"
__date__ = "2017/01"

'''This script aims at align a protein sequence on a 3D structure. Only alpha 
carbons (CA) are considered.'''


#IMPORTS
import argparse #generate a parser for arguments of the program 
import re #for the use of regular expressions
import numpy as np #for the handling of matrices
import matplotlib.pyplot as plt #for the window of visualization
import mpl_toolkits.mplot3d #for 3d dynamic rendering


#PARSER
parser= argparse.ArgumentParser() #initiate the parser
#Necessary arguments for the program
parser.add_argument("FASTA_File", type=argparse.FileType('r'),\
help = "Path to Protein Sequence.")
parser.add_argument("PDB_File", type=argparse.FileType('r'),\
help = "Path to 3D template.")
parser.add_argument("DOPE_File", type=argparse.FileType('r'),\
help = "Path to potential energies associated to distances\
between atoms of amino acid.")
parser.add_argument('Output_File',type=str,\
help='Name of the output file')
#Facultative Argument for the program 
parser.add_argument('-g','--gap_score',type=float, default = 0,\
help='Value to apply a gap score for the alignment (has to be an integer or a float > 0)')
'''By default, the Gap score is set to 0.'''

args=parser.parse_args()

#Create and Open output file
output=open(args.Output_File,'w') #Will contain all results

#REGEX
'''Regular Expression used to extract informations necessary for the 3d
structure from the PDb file.'''
motifPdb = \
r"^ATOM\s+\d+\s+CA\s+(\w{3})\s(\w+)\s+(\w+)\s+(-*\d+.\d+)\s+(-*\d+.\d+)\s+(-*\d+.\d+)"
motifPdb_comp = re.compile(motifPdb)
'''Regular Expression used to extract interaction energy levels between CAs of
amino acids from the DOPE file.'''  
motifDope = r"^(\w{3})\sCA\s(\w{3})\sCA(\s.*)" 
motifDope_comp = re.compile(motifDope)


#LISTS & DICTS

#For the Protein Sequences
aaSeq=[] #List of amino acids extracted from the FASTA file

#For the 3D Structure 
aaList=[] #List of amino acids extracted from the PDB file
aaNum=[] #List of amino acids numbers extracted from the PDB file
aaXcoord=[] #List of amino acids X coordinates
aaYcoord=[] #List of amino acids Y coordinates
aaZcoord=[] #List of amino acids Z coordinates
dicoDope={} #Dictionary of energy levels associated to a tuple of 2 amino acids
distList=[] #List of distances calculated between CAs of the 3D structure
dicoDist={} #Dictionary of distances associated to a tuple of 2 CA positions

#Parse FASTA file:
#Convert amino acids from 1 letter code to 3 leters code with dicoAA
dicoAA = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
     'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

for line in args.FASTA_File:
    if line[0]!= '>': #FASTA Header
        for aa in line:
            if aa!='\n': #FASTA Sequence
                if aa in dicoAA:
                    aaSeq.append(dicoAA[aa])
args.FASTA_File.close() #Close FASTA file

#Parse PDB file
print >> output,'#STRUCTURE_SUMMARY:'

for line in args.PDB_File:
    retourPdb = motifPdb_comp.match(line)
    if retourPdb != None : #return lines matching the regular expression
        '''print a table resuming important informations extracted from the PDB
        file :
        group(1) = Amino acid name; 
        group(2) = Protein chain name;
        group(3) = Amino acid number;
        group(4) to group(6) = XY & Z coordinates'''
        print >> output, retourPdb.group(1) + "\t" + retourPdb.group(2) + "\t"\
        + retourPdb.group(3) + "\t" + retourPdb.group(4) + "\t"\
        + retourPdb.group(5) + "\t" + retourPdb.group(6)
         
        aaList.append(retourPdb.group(1)) #add amino acid names to aaList
        aaNum.append(retourPdb.group(3)) #add amino acid numbers to aaNum
        aaXcoord.append(float(retourPdb.group(4))) #add X coordinates
        aaYcoord.append(float(retourPdb.group(5))) #add Y coordinates
        aaZcoord.append(float(retourPdb.group(6))) #add Z coordinates
args.PDB_File.close() #Close PDB file

print >> output,' '

matCoord=np.column_stack([aaXcoord,aaYcoord,aaZcoord])
'''Stack of X,Y & Z coordinates of CAs by column into a matrix.'''

#Parse DOPE file
for line in args.DOPE_File:
    retourDope = motifDope_comp.match(line)
    if retourDope != None : #return lines matching the regular expression
        '''Fill the dicoDope dictionary with two amino acids as key
        and 30 energy levels as a value.'''
        dicoDope[retourDope.group(1),retourDope.group(2)]=\
        map(float,retourDope.group(3).split()) 
        '''split() function allow to separate each energy level. map() function
        allow to convert energy levels of string format to float format.'''
args.DOPE_File.close()

#CA distances calculus
for ca in range(len(matCoord[:-3,0])):
    '''The iteration stops 3 steps before because the minimum distance
    calculated is between 2 amino acids separatedby 3 positions.'''
    distList.append(np.ndarray.tolist(np.sqrt(\
    (matCoord[ca,0]-matCoord[ca+3:,0])**2\
    + (matCoord[ca,1]-matCoord[ca+3:,1])**2\
    + (matCoord[ca,2]-matCoord[ca+3:,2])**2)))
    '''Distances are calculated thanks to the square root of the matricial sum
    of differences between each coordinates (x1-x2 ; y1-y2 ; z1-z2) each to the
    power of two.'''
    
#Convert list of lists to dictionary of distances
posOne=0
for distances in distList:
    posOne+=1
    for distance in range(len(distances)):
        dicoDist[posOne,posOne+distance+3]=distances[distance] 
        '''Go through the list of distances between CAs, and fill a dictionary
        of distances where positions of amino acids are taken as key, and the
        value is the corresponding distance calculated between those 2 amino
        acid's CAs.''' 

#Low Level Matrix
def low_level_matrix(fixed_pos, fixed_AA, AA_sequence, total_positions,\
CA_distances, DOPE_values, gap_score):
    mat = np.zeros((total_positions + 1, len(AA_sequence) + 1))
    '''initiate a matrix with as many rows as there are positions on the 3D 
    structure and as many columns as there are amino acids in the protein 
    sequence, plus one row and one column full of 0 corresponding to the 
    initialization of the alignment matrix.'''
    for i in xrange(1, len(AA_sequence) + 1):
        for j in xrange(1, total_positions + 1):
            
            left = mat[j, i - 1] + gap_score
            '''Add gap score to the value coming from the left cell.'''
            up = mat[j - 1, i] + gap_score
            '''Add gap score to the value coming from the upper cell.'''
            
            if j <= fixed_pos - 4 or j >= fixed_pos + 4 and i != fixed_AA:
                if CA_distances.has_key((fixed_pos, j)):
                    d = CA_distances[fixed_pos, j]
                else:
                    d = CA_distances[j, fixed_pos] 
                    '''Switch between the amino acids to find the corresponding
                    distance if the first combination doesn't match any key in
                    the distances dictionary'''
                    
                if d < 15.25:
                    diag = mat[j-1, i-1] + \
                    DOPE_values[AA_sequence[fixed_AA - 1], AA_sequence[i - 1]]\
                    [(int(round(d)*2)-1)] 
                    '''formula to get the corresponding energy level in the DOPE
                    dictionary values.'''
                else:
                    diag = mat[j-1, i-1] + 0
                    '''the energy level is supposed to be 0 when amino acids are
                    to distant.'''
                
                mat[j, i] = min(left, up, diag) 
                '''minimum value is taken to fill the next matrix cell.'''
                
            else:
                mat[j, i] = min(left, up, mat[j-1, i-1])
                '''minimum value is taken to fill the next matrix cell.'''
            
    np.set_printoptions(linewidth=200)
    '''Allow the matrix to be displayed without backspaces.'''
    
    return mat[-1, -1] 
    '''return the last value at the bottom right corner of the low level matrix
    (will be used to fill the high level matrix).'''

#High Level Matrix
def high_level_matrix(AA_sequence, total_positions, CA_distances, DOPE_values,\
gap_score):
    
    #2 matrices are created : a matrix with scores and another with directions
    
    #Define Scores Matrix
    matscore = np.zeros((total_positions + 1, len(AA_sequence) + 1))
    '''initiate a matrix with as many rows as there are positions on the 3D 
    structure and as many columns as there are amino acids in the protein 
    sequence, plus one row and one column full of 0 corresponding to the 
    initialization of the alignment matrix.'''
    
    #Define Directions Matrix
    matdir = np.chararray((total_positions + 1, len(AA_sequence) + 1))
    matdir[:] = ''
    '''initiate a matrix of strings with as many rows as there are positions on
    the 3D structure and as many columns as there are amino acids in the protein
    sequence, plus one row and one column full of empty strings corresponding to
    the initialization of the alignment matrix.'''
    
    directions = ['l', 'u', 'd'] #Direction options : left, up & diagonal
    
    for i in xrange(1, len(AA_sequence) + 1):
        matdir[0,i] = 'l' #first row is filled with 'l' for left origin
    for j in xrange(1, total_positions + 1):
        matdir[j,0] = 'u' #first column is filled with 'u' for up origin
        
    for i in xrange(1,len(AA_sequence)+1):
        for j in xrange(1,total_positions+1):
            left = matscore[j,i-1] + gap_score
            '''Add gap score to the value coming from the left cell.'''
            up = matscore[j-1,i] + gap_score
            '''Add gap score to the value coming from the upper cell.'''
            diag = matscore[j-1,i-1] + low_level_matrix(j, i, AA_sequence,\
            total_positions, CA_distances, DOPE_values, gap_score)
            '''Add value returned by the low level matrix to the value comind
            from the upper left cell, at the same position of the fixed amino
            acid used for the low level matrix.'''
            
            matscore[j,i] = min(left, up, diag)
            '''minimum value is taken to fill the next matrix cell.'''
            
            if left == matscore[j,i]:
                matdir[j,i] = matdir[j,i] + 'l'
            if up == matscore[j,i]:
                matdir[j,i] = matdir[j,i] + 'u'
            if diag == matscore[j,i]:
                matdir[j,i] = matdir[j,i] + 'd'
            '''multiple if are used to test every possible directions and to
            concatenate possible directions in a string'''
    
    #print both high level matrix representations and the alignment score
    print >> output,'#ALIGNMENT_DIRECTION_MATRIX:'
    print >> output,matdir
    print >> output,' '
    print >> output,'#ALIGNMENT_SCORE_MATRIX:'
    print >> output,matscore
    print >> output,' '
    print >> output,'#ALIGNMENT_SCORE: ', matscore[-1,-1]
    print >> output,' '
    
    return matdir #return the directions matrix for the backtracking function

#Backtracking
def backtracking(matrix_direction, total_positions, AA_sequence,\
matrix_coordinates):
    '''i and j are initialized at there hightest value to go throught the
    direction matrix by starting at the bottom right cell.'''
    i = len(AA_sequence)
    j = total_positions
    
    dicoValid={} 
    '''Initialization of a dictionary that will contain amino acid positions as
    keys, and there corresponding 3 letters code as values.'''
    
    while i > 0 or j > 0: #0 represents the end of the iterations
        if 'd' in matrix_direction[j,i]:
            dicoValid[j]=AA_sequence[i-1]
            i = i-1
            j = j-1
            '''i & j starting at their hightest values, they have to decrease at
            each iteration.'''
            
        elif matrix_direction[j,i][0] == 'u':
            j = j-1 #go to the cell above
        
        elif matrix_direction[j,i][0] == 'l':
            i = i-1 #go to the cell at the left
            
    print >> output,'#AMINO_ACID_POSITIONS:'
    print >> output,dicoValid 
    '''print the dictionary of positions where an amino acid has been placed.'''
    print >> output,' '
    
    #Occupied Coordinates Matrix
    occupCoord = np.empty((0,3), float)
    '''Create an empty matrix of floats with 3 columns and no rows that will
    contain 3D coordinates of positions occupied by an amino acid.'''
    
    for i in range(len(matrix_coordinates[:,0])): #number of positions
        if dicoValid.has_key(i+1): #if the position is occupied by an amino acid
            occupCoord = np.append(occupCoord,\
            np.array([matrix_coordinates[i,:]]), axis=0)
            '''3D coordinates of the position are added to matrix of occupied
            coordinates.'''
    print >> output,'#OCCUPIED_COORDINATES:'
    print >> output,occupCoord 
    '''print the matrix containing all occupied coordinates.'''
    
    return occupCoord #return the validated coordinates matrix


#PRINTERS
print >> output,'#SEQUENCE:'
print >> output,aaSeq #List of amino acids extracted from the FASTA file
print >> output,' '
print >> output,'#STRUCTURE:'
print >> output,aaList #List of amino acids extracted from the PDB file
print >> output,' '
print >> output,'#POSITIONS:'
print >> output,aaNum #List of amino acids numbers extracted from the PDB file
print >> output,' '
print >> output,'#COORDINATES:'
print >> output,matCoord 
'''Matrix of 3D coordinates of each positions extracted from the PDB file.'''
print >> output,' '

#Creation of a High Level Matrix
matAlign = high_level_matrix(aaSeq, len(aaList), dicoDist, dicoDope,\
args.gap_score)

#Creation of the Occupied Coordinates Matrix
occpdCoordinates=backtracking(matAlign, len(aaList), aaSeq, matCoord)

output.close() #close the output file

print 'DONE!' #printed when the alignment is done

#3D RENDERING
fig = plt.figure(figsize=plt.figaspect(0.5)) #creation of the figure

#First Subplot
ax = fig.add_subplot(1, 2, 1, projection='3d') 
ax.set_aspect("auto")
ax.set_autoscale_on(True) #Auto-scaling turned on
#Set axis labels
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
#Plot positions as small red spheres
for coord in range(len(matCoord)):
    ax.scatter(matCoord[coord,0],matCoord[coord,1],
    matCoord[coord,2], c='r', marker='o',
    depthshade=True,s=50)
#Plot red lines between positions
for coord in range(len(matCoord)-1):
    ax.plot(\
    [matCoord[coord,0],matCoord[coord+1,0]],\
    [matCoord[coord,1],matCoord[coord+1,1]],\
    [matCoord[coord,2],matCoord[coord+1,2]], c='r', linewidth=3)

#Second Subplot
ax = fig.add_subplot(1, 2, 2, projection='3d') 
ax.set_aspect("auto")
ax.set_autoscale_on(True) #Auto-scaling turned on
#Set axis labels
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
#Plot positions as small red spheres
for coord in range(len(matCoord)):
    ax.scatter(matCoord[coord,0],matCoord[coord,1],
    matCoord[coord,2], c='r', marker='o',
    depthshade=True,s=40)
#Plot red dashed lines between positions
for coord in range(len(matCoord)-1):
    ax.plot(\
    [matCoord[coord,0],matCoord[coord+1,0]],\
    [matCoord[coord,1],matCoord[coord+1,1]],\
    [matCoord[coord,2],matCoord[coord+1,2]],\
    c='r', linestyle='dashed', linewidth=2)
#Plot CAs as small green spheres
for coord in range(len(occpdCoordinates)):
    ax.scatter(occpdCoordinates[coord,0],occpdCoordinates[coord,1],
    occpdCoordinates[coord,2], c='g', marker='o',
    depthshade=True,s=50)
#Plot green lines between CAs
for coord in range(len(occpdCoordinates)-1):
    ax.plot(\
    [occpdCoordinates[coord,0],occpdCoordinates[coord+1,0]],\
    [occpdCoordinates[coord,1],occpdCoordinates[coord+1,1]],\
    [occpdCoordinates[coord,2],occpdCoordinates[coord+1,2]],\
    c='g', linewidth=3)
'''red positions are masked by green CAs if positions are occupied by a CAs.'''

plt.show() #Display the matplotlib window