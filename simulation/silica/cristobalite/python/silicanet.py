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
        self.type_override = None
        self.properties = {}

    def __str__(self):
        return 'Node "{}" at index {}'.format(self.type, self.index)

    @property
    def atom(self):
        return self.network.atom_group[self.index]

    @property
    def type(self):
        """
        Return the atom type as assigned by MDA, or the atom type manually overriden
        :return: 
        """
        if self.type_override:
            return self.type_override
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

    def assign_oxygens(self, box=None):
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
        # special case when only two neighbors are present
        if len(self.neighbors) == 2:
            for o_node in self.neighbors:
                self.assigned.append(o_node)
                o_node.assigned.append(self)
                o_node.is_assigned = True
            self.is_assigned = True
            return
        origin = self.atom.position
        coords = np.array([neighbor.atom.position - origin for neighbor in self.neighbors])
        # apply periodic boundary conditions
        if box is not None:
            tiled_box = np.tile(box, len(coords)).reshape(len(coords), 3)
            shifting = -np.sign(coords) * tiled_box
            borderline = abs(coords) - tiled_box / 2.0  # determines which coords should be shifted
            coords += shifting * (borderline > 0).astype(int)
        # loop over the X(0) and Y(1) components of the coordinates
        for i in (0, 1):
            i_coords = [xyz[i] for xyz in coords]
            o_node = self.neighbors[i_coords.index(max(i_coords))]
            self.assigned.append(o_node)
            o_node.assigned.append(self)
            o_node.is_assigned = True
        self.is_assigned = True
        # check validity of the assignements
        try:
            assert len(self.assigned) == 2
            assert self.assigned[0] != self.assigned[1]
        except AssertionError:
            raise RuntimeError('Assignment failed for {}'.format(self))


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
        if isinstance(item, numbers.Integral) or isinstance(item, slice):
            return self.nodes[item]
        elif isinstance(item, mda.core.groups.Atom):
            return self._ix2node[item.index]

    @property
    def molname(self):
        """
        :return: Name of the residue all silica atoms belong to 
        """
        return self.nodes[0].atom.resname

    def override_node_types(self, type_map):
        """
        Correct automatic assignment of types performed by MDAnalysis
        :param type_map: from atom name to atom type
        :return: 
        """
        for node in self:
            if node.atom.name in type_map:
                node.type_override = type_map[node.atom.name]

    def add_node_property(self, name, name2prop):
        """
        Insert a property value depending on the atom name
        :param name: name of the property
        :param name2prop: dictionary mapping atom name to property value
        """
        for node in self:
            node.properties[name] = name2prop[node.atom.name]

    def set_connections(self, cutoff, box=None):
        """
        Find which nodes are connected with the help of a contact matrix
        :param cutoff: maximum distance for connecting two nodes
        :return: 
        """
        xyz = self.atom_group.positions  # cartesian coords of the atoms
        cm = contact_matrix(xyz, cutoff=cutoff, returntype='sparse', box=box)
        for i, indices in enumerate(cm.rows):
            node_a = self.nodes[i]
            for j in indices:
                if i == j:
                    continue
                node_b = self.nodes[j]
                node_a.insert_neighbor(node_b)

    def assign_oxygens(self, box=None):
        """
        assign oxygen for every Si atom in the silica
        :param box: dimensions of the box, for periodic boundary conditions
        :return: 
        """
        for node in self:
            if node.type == 'Si':
                node.assign_oxygens(box=box)

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

