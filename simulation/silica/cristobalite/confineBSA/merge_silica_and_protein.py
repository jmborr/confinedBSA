import MDAnalysis as mda
from pdb import set_trace as tr

input_protein_structure_file = 'stripped_system_h0.17.pdb'
output_protein_structure_file = 'stripped_clean_system_h0.17.pdb'
last_protein_atom_index = 9172

u = mda.Universe(input_protein_structure_file)

# Atomic numbers should start at 1 and with no gaps
for i, atom in enumerate(u.atoms):
    atom.id = 1 + i

# Residue numbers should start at 1 and with no gaps
for i, residue in enumerate(u.residues):
    residue.resid = 1 + i

# Change the residue name of water from "SOL" to "SPC"
for atom in u.select_atoms('resname SOL'):
    atom.residue.resname = 'SPC'

# Change the segment id to 'A'
u.segments[0].segid = 'A'

u.atoms.write(output_protein_structure_file)
last_system_atom_id = u.atoms[-1].id
last_system_residue_resid = u.residues[-1].resid

# MDAnalysis cannot write PDB file with shifted 'serial'
input_silica_structure_file = 'cristobalite_21_21_15.pdb'
output_silica_structure_file = 'cristobalite_21_21_15_clean.pdb'
v = mda.Universe(input_silica_structure_file)
atomfmt = ("ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
           "{chainID:1s}{resSeq:4d}{iCode:1s}"
           "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
           "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n")

v.residues[0].resid = last_system_residue_resid + 1
header = '''HEADER    
TITLE     MDANALYSIS FRAME 0: Created by PDBWriter
REMARK     THIS IS A SIMULATION BOX
CRYST1  104.406  104.406  103.835  90.00  90.00  90.00 P 1           1
'''
buffer = ''
for i, atom in enumerate(v.atoms):
    vals = dict()
    vals['serial'] = last_system_atom_id + 1 + i
    vals['name'] = ' SI ' if atom.name=='Si' else ' OI '
    vals['altLoc'] = atom.altLoc
    vals['resName'] = 'SIL'
    vals['chainID'] = 'C'
    vals['resSeq'] = last_system_residue_resid + 1
    vals['iCode'] = ' '
    vals['pos'] = atom.position
    vals['occupancy'] = atom.occupancy
    vals['tempFactor'] = atom.tempfactor
    vals['segID'] = ''
    vals['element'] = atom.type
    buffer +=  atomfmt.format(**vals)      
open(output_silica_structure_file, 'w').write(header + buffer + 'END\n')

# Concatenate the protein and silica
output_merged = 'SiO2_protein_0.17.pdb'
protein_lines = ''.join(open(output_protein_structure_file).readlines()[:-1])
open(output_merged, 'w').write(protein_lines + buffer)

