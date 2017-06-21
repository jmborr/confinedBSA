#!/usr/bin/env/python
from __future__ import print_function
import MDAnalysis as mda
import numpy as np
from silicanet import SilicaNetwork

if __name__ == '__main__':
    u = mda.Universe('SiO2_protein_0.17.pdb')
    ntw = SilicaNetwork(u.select_atoms('resname SIL'))
    ntw.set_connections(1.7)  # generous bond distance between Si and O# generous bond distance between Si and O
    probe_count = 20000  # 10000 nodes is enough for statistics
    bond_Si_O = 0.0
    angle_O_Si_O = 0.0;  n_angle_O_Si_O = 0
    angle_Si_O_Si = 0.0; n_angle_Si_O_Si = 0
    for node in ntw[:probe_count]:
        if len(node.neighbors) < 2:
            continue
        xyz = node.atom.position
        end1 = node.neighbors[0].atom.position -xyz
        end2 = node.neighbors[1].atom.position -xyz
        # Si-O bond distance
        bond_Si_O += np.linalg.norm(end1)
        # Angle
        cosine_angle = np.dot(end1, end2) / (np.linalg.norm(end1) * np.linalg.norm(end2))
        angle = np.arccos(cosine_angle)
        if node.type == 'Si':
            angle_O_Si_O += angle; n_angle_O_Si_O += 1
        else:
            angle_Si_O_Si += angle; n_angle_Si_O_Si += 1
    print('bond_Si_O = ', bond_Si_O/probe_count)
    print('angle_O_Si_O = ', np.degrees(angle_O_Si_O/n_angle_O_Si_O))
    print('angle_Si_O_Si = ', np.degrees(angle_Si_O_Si/n_angle_Si_O_Si))
