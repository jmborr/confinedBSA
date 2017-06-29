from __future__ import print_function
import bisect
import itertools

def find_bonding_info(types, bonding_info):
    """
    Find the parameters for the particular pair of atom types
    :param types: sequence containing either two or three atom types
    :param bonds_info: either bonds or angles
    :return: dictionary with bonding parameters
    """
    ttypes = tuple(types)
    try:
        return bonding_info[ttypes]
    except KeyError:
        try:
            return bonding_info[tuple(reversed(ttypes))]
        except KeyError:
            print(bonding_info)
            raise KeyError(str(ttypes)+'\n'+str(tuple(reversed(ttypes))))

class GmxTopMixin(object):
    """
    Mixin class to create a topology file suitable for Gromacs
    """
    def topology(self, bonds_info, angles_info):
        """
        Generate include topology from the conections
        :param bonds_info: dictionary with bond parameters for each pair of atom types
        :param angles_info: dictionary with angle parameters for each triad of atom types
        :return: include topology (itp) as a string
        """
        return self.header() + self.defaults() +\
               self.moleculetype() + self.atoms() +\
               self.bonds(bonds_info) + self.angles(angles_info)

    def header(self):
        """
        :return: comment stating this class 
        """
        return """;
; Created by GmxTopMin
;
\n"""

    def defaults(self):
        """
        :return: string for default settings 
        """
        return """[ defaults ]
; nbfunc    comb-rule    gen-pairs    fudgeLJ    fudgeQQ
     1           1           no           1.0        1.0
\n"""

    def moleculetype(self):
        """
        :return: Molecule name and exclusions
        """
        return """[ moleculetype ]
; name  nrexcl
{}   2
\n""".format(self.molname)

    def atoms(self):
        """
        Generate suitable list of atoms in the system
        :return:  
        """
        print('processing [ atoms ]')
        x = """[ atoms ]
;  nr type resnr residu atom cgnr charge      mass
"""
        qtot = 0.0
        fmt = '{nr:6d}  {type}  1  {residu}  {atom}  1  {charge:14.11f}  {mass:7.4f}  ; qtot {qtot}\n'
        for i, node in enumerate(self):
            vals = {}
            vals['nr'] = 1+i
            vals['type'] = node.atom.name
            vals['residu'] = self.molname
            vals['atom'] = node.atom.name
            vals['charge'] = node.properties['charge']
            qtot += node.properties['charge']
            vals['qtot'] = qtot
            vals['mass'] = node.properties['mass']
            x += fmt.format(**vals)
        return x + '\n'

    def bonds(self, bonds_info):
        """
        :param bonds_info: dictionary with bond parameters for each pair of atom types
        :return: [ bonds ] section for gromacs include topology file
        """
        print('processing [ bonds ]')
        x = """[ bonds ]
;  i       j  funct  length      force.c.      length      force.c.
"""
        bond_entries = list()
        for node in self:
            if node.type == 'O':
                continue
            i = 1 + node.atom.index
            for neighbor in node.neighbors:
                j = 1 + neighbor.atom.index
                parameters = find_bonding_info((node.atom.name, neighbor.atom.name), bonds_info)
                bisect.insort_left(bond_entries, (sorted([i, j]), parameters))
        for ((i,j), p) in bond_entries:
            x += '{0:6d}  {1:6d}  1  {r}  {k}  {r}  {k}\n'.format(i, j, **p)
        return x + '\n'

    def angles(self, angles_info):
        """
        :param angles_info: dictionary with angle parameters for each triad of atom types
        :return: [ angles ] section for gromacs include topology file
        """
        print('processing [ angles ]')
        x = """[ angles ]
;  i       j      k   funct   angle      force.c.
"""
        angle_entries = list()
        ijks = list()
        for m, node in enumerate(self):
            if m%1000 == 0:
                print('    Remaining ', len(self)-m, ' nodes to process')
            nn = node.neighbors
            if len(nn) < 2:
                continue
            for node_triad in ((ni, node, nk) for (ni, nk) in itertools.combinations(nn, 2)):
                i, j, k = (1 + n.atom.index for n in node_triad)
                if i==k:
                    print('ouch!')
                ijk = (i, j, k) if i < k else (k, j, i)
                if ijk in ijks:
                    continue
                bisect.insort_left(ijks, ijk)
                parameters = find_bonding_info((n.atom.name for n in node_triad), angles_info)
                bisect.insort_left(angle_entries, (ijk, parameters))
        for ((i, j, k), p) in angle_entries:
            x += '{0:6d}  {1:6d}  {2:6d}  1  {theta}  {k}  {theta}  {k}\n'.format(i, j, k, **p)
        return x + '\n'

