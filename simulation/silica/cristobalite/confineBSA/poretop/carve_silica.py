#!/usr/bin/env/python
from __future__ import print_function
import MDAnalysis as mda
import numpy as np
from silicanet import SilicaNetwork, contact_sample

if __name__ == '__main__':
    # Van der Waals radii for elements in the system
    vdw_radii = {'C': 1.52, 'O': 1.55, 'H': 1.2, 'Si': 2.1, 'S': 1.8, 'N': 1.7}
    u = mda.Universe('SiO2_protein_0.17.pdb')
    # Rename silica atoms in accordance with the force field:
    type2name = {'Si': 'SIO', 'O': 'OSI'}
    for atom in u.select_atoms('resname SIL'):
        atom.name = type2name[atom.type]
    SiOd = 1.7  # generous bond distance between Si and O
    maxcutoff = 2 * vdw_radii['Si']
    neighbor_silica = u.select_atoms('(around {} (not resname SIL))'.format(maxcutoff + 4*SiOd))
    ntw = SilicaNetwork(neighbor_silica)
    ntw.set_connections(SiOd)
    for atom in neighbor_silica:
        node = ntw[atom]
        if node.type == 'Si' and len(node.neighbors) == 4:
            node.assign_oxygens()

    sample = u.select_atoms('not resname SIL')
    # carve allowing for different overlaps
    for overlap in np.arange(1.0, 2.6, 0.5):
        for atom in neighbor_silica:
            node = ntw[atom]
            node.overlaps = True if contact_sample(atom, sample, vdw_radii, overlap=overlap) else False
            node.removed = False  # clear the flag for the next section
        # Remove the SiO2 molecules having some of its atoms overlapping protein+water
        for node in ntw.get_state('overlaps'):
            si_node = node if node.type == 'Si' else node.assigned[0]
            si_node.removed = True
            for assigned in si_node.assigned:
                assigned.removed = True
        # Rename O's in accordance with the force field:
        for node in ntw:
            if node.type == 'O':
                node.atom.name = 'OA' if any(n.removed for n in node.neighbors) else 'OSI'
        removed_atoms = [node.atom for node in ntw.get_state('removed')]
        indices = sorted([atom.index for atom in removed_atoms], reverse=True)
        ul = list(u.atoms)
        for index in indices:
            del ul[index]
        carved = mda.core.AtomGroup.AtomGroup(ul)
        outpdb = 'SiO2carved_ovl{:3.1f}_protein_0.17.pdb'.format(overlap)
        carved.write(outpdb)
        print('Wrote ', outpdb)
