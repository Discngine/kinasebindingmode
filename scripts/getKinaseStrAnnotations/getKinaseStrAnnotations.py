import math
import sys
import argparse
from os import path
from os import mkdir
from Bio import PDB
import gzip
import shutil
import csv
import pandas

# Start with some argument parsing and checks
parser = argparse.ArgumentParser()
# parser.add_argument("d", help="D of DFG pattern (String), formatted as MODEL:CHAIN:RESNUM", type=str)
parser.add_argument("f", help="F of DFG pattern (String), formatted as MODEL:CHAIN:RESNUM", type=str)
parser.add_argument("e", help="Conserved E in the alphaC helix (String), formatted as MODEL:CHAIN:RESNUM", type=str) 
parser.add_argument("k", help="Front Cleft Beta Sheet LYS (String), formatted as MODEL:CHAIN:RESNUM", type=str) 
parser.add_argument("PDB_File", help="Full file path to the PDB structure file to process (String)", type=str)
# optional args
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-hac", "--helixAlphaCDist", help="add the DFG alphaC helix distance to the main output line", action="store_true")
parser.add_argument("-FS", "--fieldSeparator", help="replace the default filed separator (in non verbose mode). Default is a space", nargs='?', const=' ', type=str,  action="store", default=' ')
parser.add_argument("-d", "--download", help="allow the download of missing pdb file in a pdb subfolder. Ensure that PDB_File argument is a simple 4 letter pdb code instead a woole file name", action="store_true")
parser.add_argument("-csv", "--loadFromCsvFile", help="use a csv file to perform the analysis on multiple entries. When this option is active, d, f, e, k, and PDB_File parameters are column names instead of values", action="store")
args = parser.parse_args()



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

def CheckParameters_andPreparePDB(f, e, k, PDB_File):

    if args.verbose:
        print("verbosity turned on")
        print("F from DFG is : "+ f)
        print("E from helix AlphaC is : "+ e)
        print("K from front clft Bsheet : "+ k)
        print("PDB file Path : -"+ PDB_File+"-")

        if path.isfile(PDB_File):
            print("File is found")
        if args.download and path.isdir("pdb") and path.isfile("pdb/pdb"+PDB_File+'.ent'):
            print("File "+PDB_File+".ent is found in the  pdb/ subfolder")
        elif args.download :
            print("File "+PDB_File+" will be downloaded in the pdb/ subfolder")

    if args.download :
        if not path.isdir("pdb"):
            if args.verbose:
                print("Creating pdb subfolder")
            mkdir("pdb")
        if(path.isfile("pdb/pdb"+PDB_File+'.ent')):
            PDB_File = "pdb/pdb"+PDB_File+'.ent';
        else:
            if args.verbose:
                print("Ready to download structure")
            pdbl = PDB.PDBList()
            pdbl.download_pdb_files([PDB_File.upper()], pdir='pdb', file_format='pdb')
            PDB_File = "pdb/pdb"+PDB_File+'.ent';
     
    # if pdb file is not found, exit
    if not path.isfile(PDB_File):
        if args.verbose:
            print('Error! file '+PDB_File+' not found!')
        return None
    # basic check on secondresidue definition
    f_split=f.split(':')
    if not StringRepresentsInt(f_split[0]):
        if args.verbose:
            print('Error with argument f format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
        return None
        
    if len(f_split) != 3:
        if args.verbose:
            print('Error with argument f format. format should be MODEL:CHAIN:RESNUM')
        return None
    # basic check on third residue definition
    e_split=e.split(':')
    if not StringRepresentsInt(e_split[0]):
        if args.verbose:
            print('Error with argument e format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
        return None
    if len(e_split) != 3:
        if args.verbose:
            print('Error with argument e format. format should be MODEL:CHAIN:RESNUM')
        return None
    # basic check on fourth residue definition
    k_split=k.split(':')
    if not StringRepresentsInt(k_split[0]):
        if args.verbose:
            print('Error with argument k format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
        return None
    if len(k_split) != 3:
        if args.verbose:
            print('Error with argument k format. format should be MODEL:CHAIN:RESNUM')
        return None
    if PDB_File.endswith('.gz') :
        if args.verbose:
            print("Ready to unzip the pdb file")
        with gzip.open(PDB_File, 'r') as f_in, open('_localPdbFile.pdb', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        PDB_File = '_localPdbFile.pdb'
    return PDB_File


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

def getResidueBasedOnRefResidueAndShift(stringDefOfResidue, shift):
    splitStr=stringDefOfResidue.split(':')
    if not StringRepresentsInt(splitStr[0]):
        sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM and here MODEL is not a number')
    if len(splitStr) != 3:
        sys.exit('Error with argument d format. format should be MODEL:CHAIN:RESNUM')

    model=splitStr[0]
    chain=splitStr[1]
    residue=int(splitStr[2]) + int(shift)
    return model+':'+chain+':'+str(residue)

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

def getAlphaCHelixType(strct, f, e):
    atom1 = getResidueAtom(strct, f, "CA")
    atom2 = getResidueAtom(strct, e, "CA")

    if atom1 is None or atom2 is None :
        return 'HaC:na'

    if args.verbose:
        print('PHE atom :')
        print(atom1)
        print('GLU atom :')
        print(atom2)

    atomVect1 = getAtomPosition(atom1)
    atomVect2 = getAtomPosition(atom2)

    if args.verbose:
        print('PHE of dfg position :')
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
    atom_K = getResidueAtom(strct, k, "CA")
    atom_string_id_e_plus_4 = getResidueBasedOnRefResidueAndShift(e,4) 
    atom_E_plus_4 = getResidueAtom(strct, atom_string_id_e_plus_4 , "CA")
    if atom_F is None or atom_K is None or atom_E_plus_4 is None :
        return 'DFG:na'

    atomVect_F = getAtomPosition(atom_F)
    atomVect_E = getAtomPosition(atom_E_plus_4)
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

def process(PDB_File, f, e, k):
    PDBFilePath = CheckParameters_andPreparePDB(f, e, k, str(PDB_File))
    if PDBFilePath is None:
        return "error" +  args.fieldSeparator+"error"
    if args.verbose:
        parser = PDB.PDBParser()    
    else:
        parser = PDB.PDBParser(QUIET=True)    

    structure = parser.get_structure('myPDB', PDBFilePath)

    # some checks over residue names are not required, as sequence is not constant
    #checkResidueName(structure, f, "PHE")
    checkResidueName(structure, e, "GLU")
    checkResidueName(structure, k, "LYS")

    helixType = getAlphaCHelixType(structure, f, e)
    dfgType =  getDFGType(structure, f, e, k)

    # ready for output (verbose, basic + distance, basic)
    if args.verbose:
        print('______________________')
        print('pdb: '+PDB_File)
        print('aC helix: ' + helixType)
        print('DFG: '+dfgType)
        print('______________________')
        return dfgType + args.fieldSeparator+helixType
    elif args.helixAlphaCDist:
        return dfgType + args.fieldSeparator +helixType +  args.fieldSeparator+ str( getAlphaCHelix_DFG_dist(structure, f, e) ) 
    else :
        return dfgType + args.fieldSeparator+helixType

def main():

    if args.loadFromCsvFile :
        dataframe = pandas.read_csv(args.loadFromCsvFile, sep=args.fieldSeparator)
        print(dataframe)
        print(len(dataframe))
        for index, row in dataframe.iterrows() :
            #print(row["FilePath"])
#Kinase_ID;Kinase_Name;smiles;structure_ID;pdb;alt;chain;missing_residues;ligand;allosteric_ligand;DFG;aC_helix;back;species;
            print(row[args.PDB_File],  row[args.f], row[args.e], row[args.k],process(row[args.PDB_File],  row[args.f], row[args.e], row[args.k] ))
            
    else :
        print(process(args.PDB_File,  args.f, args.e, args.k ))

main()  
