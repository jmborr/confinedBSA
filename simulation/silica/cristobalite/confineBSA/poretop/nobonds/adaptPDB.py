#!/usr/bin/env/python

"""
We treat each silica atom as its own residue. As a consequence,
we have more than 10000 residues in the system.
The silly PDB format only allows up to 9999 residues. We solve
this issue by manually creating a .gro file, which allows for
up to 99999 residues
"""

from __future__ import print_function
import re
import MDAnalysis as mda

pdbf = 'SiO2carved_ovl1.5_protein_0.17.pdb'

# retrieve atom info
u = mda.Universe(pdbf)
GRO_FMT = ('{resid:>5d}{resname:<5s}{name:>5s}{id:>5d}'
           '{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}'
           '\n')
gro = 'confined BSA, t= 0.0\n'
#natoml = '{:5d}\n'.format(len(u.atoms))
atoml = ''
iat=0
vals = dict()
last_resid = 0
for atom in u.atoms:
    iat += 1
    vals['id'] = atom.id
    vals['name'] = atom.name
    # residue name
    vals['resname'] = atom.resname
    if atom.name in ('SIO', 'OSI', 'OA'):
        vals['resname'] = atom.name
    elif atom.resname == 'SPC':
        vals['resname'] = 'SOL'
    # residue number
    vals['resid'] = atom.resid
    if vals['resname'] in ('SIO', 'OSI', 'OA'):
        last_resid += 1
        vals['resid'] = last_resid
    else:
        last_resid = atom.resid
    vals['pos'] = atom.position/10.0  # from Angstroms to nm
    atoml += GRO_FMT.format(**vals)
gro += '{:5d}\n'.format(iat)
gro += atoml
#retrieve the box size
pdb = open(pdbf).read()
RE_BOX = re.compile('CRYST1\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
xyz = [float(xi)/10.0 for xi in RE_BOX.search(pdb).groups()]
gro += ' {:9.4f} {:9.4f} {:9.4f}\n'.format(*xyz)

open('confinedBSA_0.gro', 'w').write(gro)