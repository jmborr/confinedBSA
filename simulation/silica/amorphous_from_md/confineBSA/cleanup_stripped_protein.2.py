import Bio.PDB as bPDB
from pdb import set_trace as tr

input_silica_structure_file = 'SiO2_107A_3d_NoConnect.pdb'
output_silica_structure_file = 'SiO2_107A_3d_clean.pdb'

parser = bPDB.PDBParser(PERMISSIVE=1)
structure = parser.get_structure('bsa', input_silica_structure_file)
tr()


io = PDBIO()
io.set_structure(structure)
io.save(output_silica_structure_file)
