#!/usr/bin/env/python
from __future__ import print_function
import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix
import numbers
import operator


def contact_sample(siatom, sample, vdw_radii, overlap=1.0):
    """
    Checks if atom of silica is in contact with the protein+water system
    :param siatom: silica Atom, can be Si or O
    :param sample: AtomGroup for the protein plus water
    :param vdw_radii: dictionary of Van der Waals atomic radii
    :param overlap: interpenetration between atoms, in Angstroms
    :returns True if contact found
    """
    ele1 = siatom.type
    for ele2 in ('H', 'C', 'O', 'N', 'Si'):
        cutoff = vdw_radii[ele1] + vdw_radii[ele2] - overlap
        vals = {'type': ele2, 'x': si.position[0], 'y':si.position[1], 'z':si.position[2], 'co': cutoff}
        neighbors = sample.select_atoms('type {type} and point {x} {y} {z} {co} '.format(**vals))
        if len(neighbors):
            return True  # We found some sample atoms of type ele2 in contact with siatom
    return False


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
        self.overlaps = False
        self.removed = False

    def __str__(self):
        return 'Node "{}" at index {} wrapping atom {}'.format(self.type, self.index, self.atom)

    @property
    def atom(self):
        return self.network.atom_group[self.index]

    @property
    def type(self):
        return self.atom.type

    @property
    def max_nn(self):
        if self.type == 'Si':
            return 4
        elif self.type == 'O':
            return 2

    def insert_neighbor(self, neighbor):
        """
        Include neighbor if not already in the list, and if not
        reached maximum neighbor limit. Also include self in
        the neighbor's list of neighbors
        :param neighbor:
        :return True if inserted
        """
        if len(self.neighbors) < self.max_nn:
            if len(neighbor.neighbors) < neighbor.max_nn:
                if neighbor not in self.neighbors:
                    self.neighbors.append(neighbor)
                    neighbor.neighbors.append(self)

    def exclude_neighbor(self, neighbor):
        """
        Severe the link between self and neighbor by mutually removing themselves
        from the list of neighbors
        :param neighbor: node to be excluded
        """
        self.neighbors.remove(neighbor)
        neighbor.neighbors.remove(self)

    def retain_potential(self, querying_node):
        """
        Potential for this node not to be removed
                      Scenario                                        Power
        len(neighbors)==1 overlaps                                     0
        len(neighbors)==1                                              1
        len(neighbors)==2 + overlaps + len(otherS.neighbors)==4        2
        the following two could be reversed
            len(neighbors)==2 + overlaps + len(otherS.neighbors)==3    3
            len(neighbors)==2 + len(otherS.neighbors)==4               4
        len(neighbors)==2 + len(otherS.neighbors)==3                   5
        
        type != 'O'                                                 TypeError
        :return: (int) potential, the higher the more potential to be retained
        :except: (TypeError) the node does not host an oxygen atom
        """
        ng = self.neighbors  # just a shortcut
        if self.type != 'O':
            raise TypeError('Node must be of type "O"')
        potential = -1
        if len(ng) == 1:
            if self.overlaps:
                potential = 0
            else:
                potential = 1
        elif len(ng) == 2:
            other_Si = ng[0] if ng[0] != querying_node else ng[1]
            l = len(other_Si.neighbors)
            if self.overlaps:
                if l == 4:
                    potential = 2
                elif l == 3:
                    potential = 3
            else:
                if l == 4:
                    potential = 4
                elif l == 3:
                    potential = 5
        if potential < 0:
            raise RuntimeError('Could not assign a retain potential for {}'.format(self))
        return potential

    def attempt_removal(self):
        """
        Check if this Si atom can be removed
        :except: TypeError if type of the node is not Si
        """
        if self.type != 'Si':
            raise TypeError('Atempting to remove a node of type different than "S"')
        # find the retain potential of the neighboring oxygen atoms
        ngp = [ng.retain_potential(self) for ng in self.neighbors]
        # sort neighbors using their retain potential, from lowest to highest
        sorted_neighbors = [neig for (pot, neig) in sorted(zip(ngp, self.neighbors))]
        # remove this Si node and associated oxygen(s) with smallest
        # retain potential to this Si node
        self.removed = True
        n_o_remove = 2  # number of oxygens to be removed
        n_o_remove = 1 if len(sorted_neighbors) == 3 else 2  # Some Si start with only three O neighbors
        for i in range(n_o_remove):
            o_node = sorted_neighbors[i]
            o_node.removed = True
            # We have to find the other Si atom bonded to this oxygen,
            # and severe the link to this o_node we are removing
            onn = o_node.neighbors
            if len(onn)>1:
                other_Si = onn[0] if onn[0] != self else onn[1]  # Si neighbor other than self
                o_node.exclude_neighbor(other_Si)
        # The other oxygen atoms are not neighbors of this Si node anymore
        for i in range(n_o_remove, len(sorted_neighbors)):
            o_node = sorted_neighbors[i]
            o_node.exclude_neighbor(self)


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

    def get_type(self, atype):
        """
        Get nodes of given atom type
        :param atype: valid atom type
        :return: list of nodes
        """
        return [node for node in self.nodes if node.type == atype]


if __name__ == '__main__':
    # silica block overlapping the protein plus water
    silica_protein_water = 'SiO2_protein_0.17.pdb'
    # carved silica block
    silica_holed = 'SiO2_holed.pdb'

    # Van der Waals radii for elements in the system
    vdw_radii = {'C': 1.52, 'O': 1.55, 'H': 1.2, 'Si': 1.8, 'N': 1.7}
    # minimum cutoff between silica and protein+water
    mincutoff = vdw_radii['O'] + vdw_radii['H']
    # maximum cutoff between silica and protein+water
    maxcutoff = 2 * vdw_radii['Si']

    silica_elements = ('Si', 'O')


    u = mda.Universe(silica_protein_water)
    print('Finding silica atoms close to the protein+water..')
    SiOd = 2.0  # generous bond distance between Si and O
    neighbor_silica = u.select_atoms('(around {} (not resname SIL))'.format(maxcutoff + 3*SiOd))
    print('....found {} atoms'.format(len(neighbor_silica)))

    ntw = SilicaNetwork(neighbor_silica)
    ntw.set_connections(SiOd)

    print('Finding silica atoms overlapping the protein+water..')
    sample = u.select_atoms('not resname SIL')
    for si in neighbor_silica:
        if contact_sample(si, sample, vdw_radii):
            si_node = ntw[si]
            if si_node.type == 'Si' and len(si_node.neighbors) < 3:
                raise RuntimeError('Si node {} does not have 4 neighbors'.format(si_node))
            si_node.overlaps = True  # Flag the network node
    print('....found {} overlapping atoms'.format(len(ntw.get_state('overlaps'))))
    overlapping_Si = ntw.get_state('overlaps', atype='Si')  # nodes of type S overlapping the sample

    # Remove the Si atoms overlapping protein+water, and the two associated O atoms
    for i, si_node in enumerate(overlapping_Si):
        print(len(overlapping_Si)-i, end='..')
        si_node.attempt_removal()

    removed_nodes = ntw.get_state('removed')
    unremoved_nodes = ntw.get_state('removed', invert=True)
    unremoved_atoms = [node.atom for node in unremoved_nodes]
    #unremoved_atoms = [node.atom for node in ntw.get_state('removed', invert=True)]
    unremoved_atomgroup = mda.core.AtomGroup.AtomGroup(unremoved_atoms)
    unremoved_atomgroup.write(silica_holed)

    print('Bye')
