#!/usr/bin/env/python
from __future__ import print_function
import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix
import numbers
import operator
import numpy as np

class Node(object):
    """
    Represents a neighboring silica atom plus additional
    attributes for the carving process
    """
    def __init__(self, network, index):
        """
        :param network: NeighborNetwork containing this node
        :param index: index in the list of atoms making up the network
        """
        self.network = network
        self.index = index
        self.neighbors = list()
        self.assigned = list()
        self.overlaps = False
        self.removed = False
        self.is_assigned = False

    def __str__(self):
        return 'Node "{}" at index {}'.format(self.type, self.index)

    @property
    def atom(self):
        return self.network.atom_group[self.index]

    @property
    def type(self):
        return self.atom.type

    def insert_neighbor(self, neighbor):
        """
        Include neighbor if not already in the list, and if not
        reached maximum neighbor limit. Also include self in
        the neighbor's list of neighbors
        :param neighbor:
        :return True if inserted
        """
        if neighbor not in self.neighbors:
            self.neighbors.append(neighbor)
            neighbor.neighbors.append(self)

    def assign_oxygens(self):
        """
        In crystobalite, there two alternative Si positions, Si_1 and Si_2.
        Each Si is coordinated to four oxygens with the following coordinates,
        assuming Si is positioned at the origin:
        Oxygen relative coordinates for Si_1:
            [ 1.51200104  0.19600296  0.49399948]
            [-0.3029995   0.97300339 -1.23700047]
            [-0.19599915 -1.51199722 -0.49400043]
            [-0.97399902  0.3030014   1.23699951]
        Oxygen relative coordinates for Si_2:
            [ 0.97400284  0.30400085 -1.23700047]
            [ 0.30400085  0.97400284  1.23699951]
            [-1.51199722  0.19700241 -0.49300003]
            [ 0.19600296 -1.51199722  0.49399948]
        We will assign to each Si two oxygens, the one with the maximal X coordinate
        and the one with the maximal Y coordinate. This assignment produces a unique
        partitioning of the crystal into SiO_2 molecules.
        """
        if self.type != 'Si':
            raise TypeError('Node must be of tipe Si')
        origin = self.atom.position
        coords = [neighbor.atom.position - origin for neighbor in self.neighbors]
        # loop over the X(0) and Y(1) components of the coordinates
        for i in (0, 1):
            i_coords = [xyz[i] for xyz in coords]
            o_node = self.neighbors[i_coords.index(max(i_coords))]
            self.assigned.append(o_node)
            o_node.assigned.append(self)
            o_node.is_assigned = True
        self.is_assigned = True

class SilicaNetwork(object):
    """
    The neighbor_silica AtomGroup represented as a network of atoms in contact
    """
    def __init__(self, atom_group):
        self.atom_group = atom_group
        self.nodes = [Node(self, i) for i in range(len(atom_group))]
        self._ix2node={atom.index: self.nodes[i] for i, atom in enumerate(atom_group)}

    def __len__(self):
        return len(self.atom_group)

    def __getitem__(self, item):
        """
        Fetch a node
        :param item: node index or atom 
        :return: Node object
        """
        if isinstance(item, numbers.Integral):
            return self.nodes[item]
        elif isinstance(item, mda.core.groups.Atom):
            return self._ix2node[item.index]

    def set_connections(self, cutoff):
        """
        Find which nodes are connected with the help of a contact matrix
        :param cutoff: maximum distance for connecting two nodes
        :return: 
        """
        xyz = self.atom_group.positions  # cartesian coords of the atoms
        cm = contact_matrix(xyz, cutoff=cutoff, returntype='sparse')
        for i, indices in enumerate(cm.rows):
            node_a = self.nodes[i]
            for j in indices:
                if i == j:
                    continue
                node_b = self.nodes[j]
                node_a.insert_neighbor(node_b)

    def get_state_i(self, state_attribute, invert=False, atype=None):
        """
        :param state_attribute: attribute of node ('overlaps', or 'removed')
        :param invert: consider the negative of the value of the state_attribute
        :param atype: atom type, all types if None
        :return: list of atom_group indices for nodes of type overlapping the sample
        """
        indices = list()
        for i, node in enumerate(self.nodes):
            state = getattr(node, state_attribute)
            if invert:
                state = not state
            if state:
                if (not atype) or atype == node.type:
                    indices.append(i)
        return indices

    def get_state(self, state_attribute, invert=False, atype=None):
        """
        :param state_attribute: attribute of node we query (overlaps, or removed)
        :param invert: consider the negative of the value of the state_attribute
        :param atype: atom type, all types if None
        :return: list of nodes of type overlapping the sample
        """
        indices = self.get_state_i(state_attribute, invert=invert, atype=atype)
        return operator.itemgetter(*indices)(self.nodes)


def contact_sample(si, sample, vdw_radii, overlap=1.5):
    """
    Checks if atom of silica is in contact with the protein+water system
    :param si: silica Atom, can be Si or O
    :param sample: AtomGroup for the protein plus water
    :param vdw_radii: dictionary of Van der Waals atomic radii
    :param overlap: interpenetration between atoms, in Angstroms
    :returns True if contact found
    """
    ele1 = si.type
    for ele2 in ('H', 'C', 'O', 'N', 'S'):
        cutoff = vdw_radii[ele1] + vdw_radii[ele2] - overlap
        vals = {'type': ele2,
                'x': si.position[0], 'y': si.position[1], 'z': si.position[2],
                'co': cutoff}
        neighbors = sample.select_atoms('type {type} and point {x} {y} {z} {co} '.format(**vals))
        if len(neighbors):
            return True  # We found some sample atoms of type ele2 in contact with siatom
    return False


if __name__ == '__main__':
    # Van der Waals radii for elements in the system
    vdw_radii = {'C': 1.52, 'O': 1.55, 'H': 1.2, 'Si': 2.1, 'S': 1.8, 'N': 1.7}
    u = mda.Universe('SiO2_protein_0.17.pdb')
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
        removed_atoms = [node.atom for node in ntw.get_state('removed')]
        indices = sorted([atom.index for atom in removed_atoms], reverse=True)
        ul = list(u.atoms)
        for index in indices:
            del ul[index]
        carved = mda.core.AtomGroup.AtomGroup(ul)
        carved.write('SiO2carved_ovl{:3.1f}_protein_0.17.pdb'.format(overlap))

