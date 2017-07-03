#!/usr/bin/env/python
from __future__ import print_function
import MDAnalysis as mda
from silicanet import SilicaNetwork
from silicatop import GmxTopMixin

class Network(SilicaNetwork, GmxTopMixin):
    """
    Incorporate topology generator to SilicaNetwork
    """
    pass

if __name__ == '__main__':
    u = mda.Universe('SiO2carved_ovl2.5.pdb')
    ntw = Network(u.atoms)
    ntw.override_node_types({'SIO': 'Si'})  # MDAnalysis assigns the wrong type to SIO atoms
    SiOd = 1.7  # generous bond distance between Si and O
    ntw.set_connections(SiOd, box=u.dimensions[0:3])
    ntw.assign_oxygens(box=u.dimensions[0:3])
    #q_SIO = 1.113796741791171  # Renders neutral silica system
    #q_SIO = 1.1144720182324641  # Renders neutral silica+protein+water system
    #ntw.add_node_property('charge', {'SIO': q_SIO, 'OSI': -0.55, 'OA': -0.90})
    ntw.add_node_property('charge', {'SIO': 1.1, 'OSI': -0.55, 'OA': -0.55})
    ntw.add_node_property('mass', {'SIO': 28.0855, 'OSI': 15.9994, 'OA': 15.9994})
    bonds = {('SIO', 'OSI'): {'r': 0.157833041748, 'k': 251040.0},
             ('SIO',  'OA'): {'r': 0.157833041748, 'k': 251040.0},
            }
    angles = {('OSI', 'SIO', 'OSI'): {'theta': 109.149737561, 'k': 397.480},
              ( 'OA', 'SIO', 'OSI'): {'theta': 109.149737561, 'k': 397.480},
              ('SIO', 'OSI', 'SIO'): {'theta': 146.498743845, 'k': 397.480},
              ( 'OA', 'SIO',  'OA'): {'theta': 109.149737561, 'k': 397.480},
            }
    open('junk','w').write(ntw.header() + ntw.defaults() + ntw.moleculetype() + ntw.atoms() + ntw.bonds(bonds))
    open('SiO2carved_ovl2.5.itp', 'w').write(ntw.topology(bonds, angles))


