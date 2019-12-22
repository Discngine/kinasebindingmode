import math
import sys
import argparse
from os import path
from os import mkdir
from Bio import PDB
import gzip
import shutil

# ex :  python getKinaseProperties.py  1:B:194 1:B:90 1:B:74 ../../prepared_data/structures/5l4q.pdb -v

# Start with some argument parsing and checks
parser = argparse.ArgumentParser()
parser.add_argument("d", help="D of DFG pattern (String), formatted as MODEL:CHAIN:RESNUM", type=str)
parser.add_argument("f", help="F of DFG pattern (String), formatted as MODEL:CHAIN:RESNUM", type=str)
parser.add_argument("e", help="Conserved E in the alphaC helix (String), formatted as MODEL:CHAIN:RESNUM", type=str) 
parser.add_argument("k", help="Front Cleft Beta Sheet LYS (String), formatted as MODEL:CHAIN:RESNUM", type=str) 
parser.add_argument("PDB_File", help="Full file path to the PDB structure file to process (String)", type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-hac", "--helixAlphaCDist", help="add the DFG alphaC helix distance to the main output line", action="store_true")
parser.add_argument("-d", "--download", help="allow the download of missing pdb file in a pdb subfolder. Ensure that PDB_File argument is a simple 4 letter pdb code instead a woole file name", action="store_true")
#parser.add_argument("-r", "--repository", help="suport the use of a pdb repository containing gzip pdb in subfolder (ex: wr/pdb4wrc.ent.gz)", action="store")
args = parser.parse_args()


if args.verbose:
    print("verbosity turned on")
    print("D from DFG is : "+ args.d)
    print("F from DFG is : "+ args.f)
    print("E from helix AlphaC is : "+ args.e)
    print("K from front clft Bsheet : "+ args.k)
    print("PDB file Path : -"+ args.PDB_File+"-")
    if path.isfile(args.PDB_File):
        print("File is found")
    if args.download and path.isdir("pdb") and path.isfile("pdb/pdb"+args.PDB_File+'.ent'):
        print("File "+args.PDB_File+".ent is found in the  pdb/ subfolder")
    elif args.download :
        print("File "+args.PDB_File+" will be downloaded in the pdb/ subfolder")
    #if args.repository :
    #    print("Structure repositoy used : " + args.repository)


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

def CheckParameters():
    if args.download :
        if not path.isdir("pdb"):
            if args.verbose:
                print("Creating pdb subfolder")
            mkdir("pdb")
        if(path.isfile("pdb/pdb"+args.PDB_File+'.ent')):
            args.PDB_File = "pdb/pdb"+args.PDB_File+'.ent';
        else:
            if args.verbose:
                print("Ready to download structure")
            pdbl = PDB.PDBList()
            pdbl.download_pdb_files([args.PDB_File.upper()], pdir='pdb', file_format='pdb')
            args.PDB_File = "pdb/pdb"+args.PDB_File+'.ent';
     
    # if pdb file is not found, exit
    if not path.isfile(args.PDB_File):
        sys.exit('Error! file '+args.PDB_File+' not found!')
    # basic check on first residue definition
    d_split=args.d.split(':')
    if not StringRepresentsInt(d_split[0]):
        sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
    if len(d_split) != 3:
        sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM')
    # basic check on secondresidue definition
    f_split=args.f.split(':')
    if not StringRepresentsInt(f_split[0]):
        sys.exit('Error with argument f format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
    if len(f_split) != 3:
        sys.exit('Error with argument f format. format should be MODEL:CHAIN:RESNUM')
    # basic check on third residue definition
    e_split=args.e.split(':')
    if not StringRepresentsInt(e_split[0]):
        sys.exit('Error with argument e format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
    if len(e_split) != 3:
        sys.exit('Error with argument e format. format should be MODEL:CHAIN:RESNUM')
    # basic check on fourth residue definition
    k_split=args.k.split(':')
    if not StringRepresentsInt(k_split[0]):
        sys.exit('Error with argument k format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
    if len(k_split) != 3:
        sys.exit('Error with argument k format. format should be MODEL:CHAIN:RESNUM')
    if args.PDB_File.endswith('.gz') :
        if args.verbose:
            print("Ready to unzip the pdb file")
        with gzip.open(args.PDB_File, 'r') as f_in, open('_localPdbFile.pdb', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        args.PDB_File = '_localPdbFile.pdb'


    
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
def getAtomDistance(atomVect1, atomVect2):
    distSquare = (atomVect2[0] - atomVect1[0])**2 + (atomVect2[1] - atomVect1[1])**2 + (atomVect2[2] - atomVect1[2])**2
    res = math.sqrt(distSquare)
    return res

# getAlphHelixType
# This function defines the the alpha helix type of the kinase based on the distance between 2 atoms of its atoms
# @param distance   the distance in Angstrom between 2 atoms
#
# @return   the type of the Alpha Helix (in, out, out-like)
#
def getAlphaCHelixType_fromDistance(distance):
    if distance > 4 and distance <= 7.2:
        type = 'HaC:in'
    elif distance > 7.2 and distance <= 9.3:
        type = 'HaC:out-like'
    elif distance > 9.3 and distance <= 14:
        type = 'HaC:out'
    else:
        type = 'HaC:na'
    return type


def getResidueAtom(strct, stringDefOfResidue, AtomName):
    splitStr=stringDefOfResidue.split(':')
    if not StringRepresentsInt(splitStr[0]):
        sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
    if len(splitStr) != 3:
        sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM')

    model=int(splitStr[0])-1
    chain=splitStr[1]
    residue=int(splitStr[2])
    try:
        atom = strct[model][chain][residue][AtomName]
    except KeyError:
        if args.verbose:
            print('failed to get atom '+stringDefOfResidue+' with name '+AtomName)
        atom = None
    return atom

def getAtomPosition(atom):
    return atom.get_vector()

def getResidueName(strct, stringDefOfResidue):
    splitStr=stringDefOfResidue.split(':')
    if not StringRepresentsInt(splitStr[0]):
        sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
    if len(splitStr) != 3:
        sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM')

    model=int(splitStr[0])-1
    chain=splitStr[1]
    residue=int(splitStr[2])
    res = strct[model][chain][residue]
    return res.get_resname()

def checkResidueName(strct, stringDefOfResidue, residueName):
    structResName = getResidueName(strct, stringDefOfResidue)
    if structResName == residueName:
        pass
    else:
        raise Exception('Expected residue seems wrong. ' + stringDefOfResidue + ' should be a '+ residueName + ' but a '+structResName +' was found')


def getAlphaCHelix_DFG_dist(strct, d, e):
    atom1 = getResidueAtom(strct, d, "CA")
    atom2 = getResidueAtom(strct, e, "CA")

    if atom1 is None or atom2 is None :
        return 'na'
    atomVect1 = getAtomPosition(atom1)
    atomVect2 = getAtomPosition(atom2)
    distance = getAtomDistance(atomVect1, atomVect2)
    return distance


def getAlphaCHelixType(strct, d, e):
    atom1 = getResidueAtom(strct, d, "CA")
    atom2 = getResidueAtom(strct, e, "CA")

    if atom1 is None or atom2 is None :
        return 'HaC:na'


    if args.verbose:
        print('ASP atom :')
        print(atom1)
        print('GLU atom :')
        print(atom2)

    atomVect1 = getAtomPosition(atom1)
    atomVect2 = getAtomPosition(atom2)

    if args.verbose:
        print('ASP of dfg position :')
        print(atomVect1)
        print('GLU of alpha C helix position :')
        print(atomVect2)
    distance = getAtomDistance(atomVect1, atomVect2)
    if args.verbose:
        print('Distance : ')
        print(distance)
        print('Helix type :')
        print(getAlphaCHelixType_fromDistance(distance))
    return getAlphaCHelixType_fromDistance(distance)

def getDFGType(strct,f,e,k):
    atom_F = getResidueAtom(strct, f, "CZ")
    atom_E = getResidueAtom(strct, e, "CA")
    atom_K = getResidueAtom(strct, k, "CA")

    if atom_F is None or atom_F is None or atom_F is None :
        return 'DFG:na'

    atomVect_F = getAtomPosition(atom_F)
    atomVect_E = getAtomPosition(atom_E)
    atomVect_K = getAtomPosition(atom_K)
    
    D1 = getAtomDistance(atomVect_F, atomVect_E)
    D2 = getAtomDistance(atomVect_F, atomVect_K)

    if D1 <= 11 and D2 <= 11:
        type = 'DFG:inter'
    elif D1 > 11 and D2 <= 14:
        type = 'DFG:out'
    else:
        type = 'DFG:in'
    if args.verbose:
        print('Distance 1 : '+str(D1))
        print('Distacnce2 : '+str(D2))
        print('DFG type : '+type)
    return type    

def main():
    CheckParameters()
    if args.verbose:
        parser = PDB.PDBParser()    
    else:
        parser = PDB.PDBParser(QUIET=True)    

    structure = parser.get_structure('myPDB', args.PDB_File)
    checkResidueName(structure, args.f, "PHE")
    checkResidueName(structure, args.d, "ASP")
    checkResidueName(structure, args.e, "GLU")
    checkResidueName(structure, args.k, "LYS")


    helixType = getAlphaCHelixType(structure, args.d, args.e)
    dfgType =  getDFGType(structure, args.f, args.e, args.k)
    if args.helixAlphaCDist:
        dfg_helixAc_dist = getAlphaCHelix_DFG_dist(structure, args.d, args.e)

    if args.verbose:
        print('______________________')
        print('pdb: '+args.PDB_File)
        print('aC helix: ' + helixType)
        print('DFG: '+dfgType)
        print('______________________')
    elif args.helixAlphaCDist:
        print (dfgType +' '+helixType + ' '+ str( dfg_helixAc_dist)   )
    else :
        print (dfgType +' '+helixType);
main()  
