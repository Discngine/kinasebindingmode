import math
import sys
import argparse
from os import path
from Bio import PDB


# Start with some argument parsing and checks
parser = argparse.ArgumentParser()
parser.add_argument("d", help="D of DFG pattern (String), formatted as MODEL:CHAIN:RESNUM", type=str)
parser.add_argument("e", help="Conserved E in the alphaC helix (String), formatted as MODEL:CHAIN:RESNUM", type=str) 
parser.add_argument("PDB_File", help="Full file path to the PDB structure file to process (String)", type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
args = parser.parse_args()







if args.verbose:
    print("verbosity turned on")
    print("D from DFG is : "+ args.d)
    print("E from helix AlphaC is : "+ args.e)
    print("PDB file Path : "+ args.PDB_File)
    if path.isfile(args.PDB_File):
        print("File is found")

# StringRepresentsInt
# This small helper function allow a check on a string and return True if the string can be converted into an integer
# @param s  the string to test
#
# return a boolean
#
def StringRepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

# if pdb file is not found, exit
if not path.isfile(args.PDB_File):
    sys.exit('Error! file not found!')
# basic check on first residue definition
d_split=args.d.split(':')
if not StringRepresentsInt(d_split[0]):
     sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
if len(d_split) != 3:
    sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM')
# basic check on second residue definition
e_split=args.e.split(':')
if not StringRepresentsInt(e_split[0]):
     sys.exit('Error with argument e format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
if len(e_split) != 3:
    sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM')



# getAtomDstance
# This function calculates the distance between 2 atoms based on the 3D coordinate vector of each atom
# @param atomVector1   the 3D coordinate vector of the first atom (in Angstrom)
# @param atomVector2   the 3D coordinate vector of the second atom (in Angstrom)
#
# @return   the distance in Angstrom between the 2 atoms
#
def getAtomDstance(atomVect1, atomVect2):
    distSquare = (atomVect2[0] - atomVect1[0])**2 + (atomVect2[1] - atomVect1[1])**2 + (atomVect2[2] - atomVect1[2])**2
    res = math.sqrt(distSquare)
    return res

# getAlphHelixType
# This function defines the the alpha helix type of the kinase based on the distance between 2 atoms of its atoms
# @param distance   the distance in Angstrom between 2 atoms
#
# @return   the type of the Alpha Helix (in, out, out-like)
#
def getAlphHelixType(distance):
    if distance > 4 and distance <= 7.2:
        type = 'alphaC-in'
    elif distance > 7.2 and distance <= 9.3:
        type = 'alphaC-out-like'
    elif distance > 9.3 and distance <= 14:
        type = 'alphaC-out'
    else:
        type = 'Not Applicable'
    return type


parser = PDB.PDBParser()
#structure = parser.get_structure('AAK1', '../../prepared_data/structures/5l4q.pdb')
structure = parser.get_structure('myPDB', args.PDB_File)
#models = list(structure.get_models())
#model = models[0]
#chains = list(model.get_chains())

#i = 0
#for atom in structure.get_atoms():
#    if atom.get_serial_number() == 5:
#        atomVect1 = atom.get_vector()
#    if atom.get_serial_number() == 50:
#        atomVect2 = atom.get_vector()
dm=int(d_split[0])-1
dc=d_split[1]
dr=int(d_split[2])
em=int(e_split[0])-1
ec=e_split[1]
er=int(e_split[2])

atom1 = structure[dm][dc][dr]["CA"]
atom2 = structure[em][ec][er]["CA"]

atomVect1 = atom1.get_vector()
atomVect2 = atom2.get_vector()

distance = getAtomDstance(atomVect1, atomVect2)
if args.verbose:
    print('vect1')
    print(atomVect1)
    print('vect2')
    print(atomVect2)
    print(distance)
print(getAlphHelixType(distance))
