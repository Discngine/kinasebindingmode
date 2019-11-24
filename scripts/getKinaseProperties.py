from Bio.PDB import *
import math

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


parser = PDBParser()
structure = parser.get_structure('AAK1', '../prepared_data/structures/5l4q.pdb')
i = 0
for atom in structure.get_atoms():
    if atom.get_serial_number() == 5:
        atomVect1 = atom.get_vector()
    if atom.get_serial_number() == 50:
        atomVect2 = atom.get_vector()

distance = getAtomDstance(atomVect1, atomVect2)
print('vect1')
print(atomVect1)
print('vect2')
print(atomVect2)
print(distance)
print(getAlphHelixType(distance))