from Bio.PDB import *
import math
import sys

# getCzetaFromDfgMotif
# This function calculates get the position from the Zeta carbon from the DFG motif
# @param chain         the chain from the pdb used for the calculation
# @param phe_number    the number of the Phenylalanine residue from the Kinase DFG motif.
#
# @return   the position of the Zeta Carbon of the phenylalamide
#
def getCzetaFromDfgMotif(chain, phe_number):
    for residue in chain.get_residues():
        if residue.get_id()[1] == phe_number and residue.get_resname()  == 'PHE':
            for atom in residue:
                if atom.get_name() == "CZ":
                    position = atom.get_vector()
        elif residue.get_id()[1] == phe_number and residue.get_resname()  != 'PHE':
            raise Exception('The residue mentionned is not a phenylalamide. The Calculation can not continue.')
    return position

# getCalphaFromAlphaHelix
# This function calculates get the position from the alpha carbon from the alpha helix
# @param chain         the chain from the pdb used for the calculation
# @param gly_number    the number of the glutamine residue from the kinase alpha helix.
#
# @return   the position of the alpha Carbon of the alpha helix
#
def getCalphaFromAlphaHelix(chain, glu_number):
    res_position = glu_number + 4 
    for residue in chain.get_residues():
        if residue.get_id()[1] == res_position:
            for atom in residue:
                if atom.get_name() == "CA":
                    position = atom.get_vector()
    return position

# getCalphaFromBeta3
# This function calculates get the position from the alpha carbon from the Beta-sheet 3
# @param chain         the chain from the pdb used for the calculation
# @param lys_number    the number of the Lysine residue from the Kinase Beta-sheet 3.
#
# @return   the position of the alpha Carbon of the lysine
#
def getCalphaFromBeta3(chain, lys_number):
    for residue in chain.get_residues():
        if residue.get_id()[1] == lys_number and residue.get_resname()  == 'LYS':
            for atom in residue:
                if atom.get_name() == "CA":
                    position = atom.get_vector()
        elif residue.get_id()[1] == lys_number and residue.get_resname()  != 'LYS':
            raise Exception('The residue mentionned is not a lysine. The Calculation can not continue.')
    return position

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
def getDFGType(D1, D2):
    if D1 <= 11 and D2 <= 11:
        type = 'DFGinter'
    elif D1 > 11 and D2 <= 14:
        type = 'DFGout'
    else:
        type = 'DFGin'
    return type    

def main():
    parser = PDBParser()
    argt = sys.argv
    if len(argt) == 5:
        pdbFile = argt[1]
        phe_res_num = int(argt[2])
        glu_res_num = int(argt[3])
        lys_res_num = int(argt[4])

        structure = parser.get_structure('KinaseName', pdbFile)
        model = structure[0]
        chain = model["A"]
        
        #getting the position of each atoms
        czeta_pos = getCzetaFromDfgMotif(chain, phe_res_num)
        helix_calpha_pos = getCalphaFromAlphaHelix(chain, glu_res_num)
        sheet_calpha_pos = getCalphaFromBeta3(chain, lys_res_num)

        D1 = getAtomDstance(czeta_pos,helix_calpha_pos)
        D2 = getAtomDstance(czeta_pos,sheet_calpha_pos)

        dfgType = getDFGType(D1,D2)
        print(D1)
        print(D2)
        
        print(dfgType)
    else:
        raise Exception('The number of argument is incorect')

    

main()  