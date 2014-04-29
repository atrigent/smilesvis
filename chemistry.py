from collections import namedtuple, defaultdict
from linalg import Vector
import math

# these are the idealized geometries for the steric
# numbers from VSEPR theory - we will omit
# bonds (where there are lone pairs) as necessary
geometries = [
    # geometry for things with only one bond - steric number 1?
    [Vector(0, 1, 0)],
    # steric number 2 - linear
    [Vector(1, 0, 0),
     Vector(-1, 0, 0)],
    # steric number 3 - trigonal planar, bent
    [Vector(0, 1, 0),
     Vector(math.sqrt(3)/2, -1/2, 0),
     Vector(-math.sqrt(3)/2, -1/2, 0)],
    # steric number 4 - tetrahedral, trigonal pyramidal, bent
    [Vector(0, 1, 0),
     Vector((2 * math.sqrt(2)) / 3, -1/3, 0),
     Vector(-math.sqrt(2)/3, -1/3, math.sqrt(2/3)),
     Vector(-math.sqrt(2)/3, -1/3, -math.sqrt(2/3))],
    # steric number 5: trigonal bipyramidal, seesaw, t-shaped, linear
    [Vector(0, 1, 0),
     Vector(math.sqrt(3)/2, -1/2, 0),
     Vector(-math.sqrt(3)/2, -1/2, 0),
     Vector(0, 0, -1),
     Vector(0, 0, 1)],
    # steric number 6: octahedral, square pyramidal, square planar
    [Vector( 1,  0,  0),
     Vector(-1,  0,  0),
     Vector( 0,  1,  0),
     Vector( 0, -1,  0),
     Vector( 0,  0,  1),
     Vector( 0,  0, -1)]
    # steric number 7: pentagonal bipyramidal, pentagonal pyramidal,
    #                  pentagonal planar
    # steric number 8: square antiprismatic
    # steric number 9: tricapped trigonal prismatic,
    #                  capped square antiprismatic
]

def get_geometry(num, lone_pairs):
    return geometries[num-1][lone_pairs:]

Element = namedtuple('Element', 'symbol name smiles_organic color hybridizations')

elements = [
    Element('H', 'Hydrogen', False, (1, 1, 1), [
        {
            'steric_num': 1,
            'lone_pairs': 0
        }
    ]),
    Element('He', 'Helium', False, None, None),
    Element('Li', 'Lithium', False, None, None),
    Element('Be', 'Beryllium', False, None, None),
    Element('B', 'Boron', True, None, None),
    Element('C', 'Carbon', True, (0, 0, 0), [
        { # sp3
            'steric_num': 4,
            'lone_pairs': 0
        }, { # sp2
            'double': 1,
            'steric_num': 3,
            'lone_pairs': 0
        }, { # sp
            'double': 2,
            'steric_num': 2,
            'lone_pairs': 0,
        }, { # also sp
            'triple': 1,
            'steric_num': 2,
            'lone_pairs': 0
        }
    ]),
    Element('N', 'Nitrogen', True, (0, 0, 1), [
        {
            'double': 1,
            'steric_num': 3,
            'lone_pairs': 1
        }, {
            'steric_num': 4,
            'lone_pairs': 1
        }
    ]),
    Element('O', 'Oxygen', True, (1, 0, 0), [
        {
            'double': 1,
            'steric_num': 3,
            'lone_pairs': 2
        }, {
            'steric_num': 4,
            'lone_pairs': 2
        }
    ]),
    Element('F', 'Fluorine', True, (0, 1, 0), [
        {
            'steric_num': 4,
            'lone_pairs': 3
        }
    ]),
    Element('Ne', 'Neon', False, None, None),
    Element('Na', 'Sodium', False, None, None),
    Element('Mg', 'Magnesium', False, None, None),
    Element('Al', 'Aluminum', False, None, None),
    Element('Si', 'Silicon', False, None, None),
    Element('P', 'Phosphorus', True, None, None),
    Element('S', 'Sulfur', True, (0xdd/0xff, 0xdd/0xff, 0), [
        {
            'double': 1,
            'steric_num': 3,
            'lone_pairs': 1
        }, {
            'steric_num': 5,
            'lone_pairs': 1
        }, {
            'steric_num': 6,
            'lone_pairs': 0
        }, {
            'double': 3,
            'steric_num': 3,
            'lone_pairs': 0
        }
    ]),
    Element('Cl', 'Chlorine', True, (0, 1, 0), [
        {
            'steric_num': 4,
            'lone_pairs': 3
        }
    ]),
    Element('Ar', 'Argon', False, None, None),
    Element('K', 'Potassium', False, None, None),
    Element('Ca', 'Calcium', False, None, None),
    Element('Sc', 'Scandium', False, None, None),
    Element('Ti', 'Titanium', False, None, None),
    Element('V', 'Vanadium', False, None, None),
    Element('Cr', 'Chromium', False, None, None),
    Element('Mn', 'Manganese', False, None, None),
    Element('Fe', 'Iron', False, None, None),
    Element('Co', 'Cobalt', False, None, None),
    Element('Ni', 'Nickel', False, None, None),
    Element('Cu', 'Copper', False, None, None),
    Element('Zn', 'Zinc', False, None, None),
    Element('Ga', 'Gallium', False, None, None),
    Element('Ge', 'Germanium', False, None, None),
    Element('As', 'Arsenic', False, None, None),
    Element('Se', 'Selenium', False, None, None),
    Element('Br', 'Bromine', True, (0.5, 0, 0), [
        {
            'steric_num': 4,
            'lone_pairs': 3
        }
    ]),
    Element('Kr', 'Krypton', False, None, None),
    Element('Rb', 'Rubidium', False, None, None),
    Element('Sr', 'Strontium', False, None, None),
    Element('Y', 'Yttrium', False, None, None),
    Element('Zr', 'Zirconium', False, None, None),
    Element('Nb', 'Niobium', False, None, None),
    Element('Mo', 'Molybdenum', False, None, None),
    Element('Tc', 'Technetium', False, None, None),
    Element('Ru', 'Ruthenium', False, None, None),
    Element('Rh', 'Rhodium', False, None, None),
    Element('Pd', 'Palladium', False, None, None),
    Element('Ag', 'Silver', False, None, None),
    Element('Cd', 'Cadmium', False, None, None),
    Element('In', 'Indium', False, None, None),
    Element('Sn', 'Tin', False, None, None),
    Element('Sb', 'Antimony', False, None, None),
    Element('Te', 'Tellurium', False, None, None),
    Element('I', 'Iodine', True, None, None),
    Element('Xe', 'Xenon', False, None, None),
    Element('Cs', 'Caesium', False, None, None),
    Element('Ba', 'Barium', False, None, None),
    Element('La', 'Lanthanum', False, None, None),
    Element('Ce', 'Cerium', False, None, None),
    Element('Pr', 'Praseodymium', False, None, None),
    Element('Nd', 'Neodymium', False, None, None),
    Element('Pm', 'Promethium', False, None, None),
    Element('Sm', 'Samarium', False, None, None),
    Element('Eu', 'Europium', False, None, None),
    Element('Gd', 'Gadolinium', False, None, None),
    Element('Tb', 'Terbium', False, None, None),
    Element('Dy', 'Dysprosium', False, None, None),
    Element('Ho', 'Holmium', False, None, None),
    Element('Er', 'Erbium', False, None, None),
    Element('Tm', 'Thulium', False, None, None),
    Element('Yb', 'Ytterbium', False, None, None),
    Element('Lu', 'Lutetium', False, None, None),
    Element('Hf', 'Hafnium', False, None, None),
    Element('Ta', 'Tantalum', False, None, None),
    Element('W', 'Tungsten', False, None, None),
    Element('Re', 'Rhenium', False, None, None),
    Element('Os', 'Osmium', False, None, None),
    Element('Ir', 'Iridium', False, None, None),
    Element('Pt', 'Platinum', False, None, None),
    Element('Au', 'Gold', False, None, None),
    Element('Hg', 'Mercury', False, None, None),
    Element('Tl', 'Thallium', False, None, None),
    Element('Pb', 'Lead', False, None, None),
    Element('Bi', 'Bismuth', False, None, None),
    Element('Po', 'Polonium', False, None, None),
    Element('At', 'Astatine', False, None, None),
    Element('Rn', 'Radon', False, None, None),
    Element('Fr', 'Francium', False, None, None),
    Element('Ra', 'Radium', False, None, None),
    Element('Ac', 'Actinium', False, None, None),
    Element('Th', 'Thorium', False, None, None),
    Element('Pa', 'Protactinium', False, None, None),
    Element('U', 'Uranium', False, None, None),
    Element('Np', 'Neptunium', False, None, None),
    Element('Pu', 'Plutonium', False, None, None),
    Element('Am', 'Americium', False, None, None),
    Element('Cm', 'Curium', False, None, None),
    Element('Bk', 'Berkelium', False, None, None),
    Element('Cf', 'Californium', False, None, None),
    Element('Es', 'Einsteinium', False, None, None),
    Element('Fm', 'Fermium', False, None, None),
    Element('Md', 'Mendelevium', False, None, None),
    Element('No', 'Nobelium', False, None, None),
    Element('Lr', 'Lawrencium', False, None, None),
    Element('Rf', 'Rutherfordium', False, None, None),
    Element('Db', 'Dubnium', False, None, None),
    Element('Sg', 'Seaborgium', False, None, None),
    Element('Bh', 'Bohrium', False, None, None),
    Element('Hs', 'Hassium', False, None, None),
    Element('Mt', 'Meitnerium', False, None, None),
    Element('Ds', 'Darmstadtium', False, None, None),
    Element('Rg', 'Roentgenium', False, None, None),
    Element('Cn', 'Copernicium', False, None, None),
    Element('Uut', 'Ununtrium', False, None, None),
    Element('Fl', 'Flerovium', False, None, None),
    Element('Uup', 'Ununpentium', False, None, None),
    Element('Lv', 'Livermorium', False, None, None),
    Element('Uus', 'Ununseptium', False, None, None),
    Element('Uuo', 'Ununoctium', False, None, None)
]

def get_element(num):
    return elements[num-1]

class Bond:
    def __init__(self, bond_type, atom):
        self.bond_type = bond_type
        self.atom = atom

    def order(self):
        return self.bond_type.split('-')[0]

    def replace_atom(self, oldatom, newatom):
        if self.atom == oldatom:
            self.atom = newatom

class Atom:
    def __init__(self, element, isotope=None, tetrahedral=None, hydrogens=None, charge=None, substituents=None):
        self.element = element
        self.isotope = isotope
        self.tetrahedral = tetrahedral
        self.hydrogens = hydrogens
        self.charge = charge
        self.bonds = substituents

        self.bond_counts = defaultdict(lambda: 0)
        for bond in self.bonds:
            self.bond_counts[bond.bond_type] += 1

            if bond.bond_type != bond.order():
                self.bond_counts[bond.order()] += 1

        hybridization = None
        if self.element.hybridizations:
            for cur_hybrid in self.element.hybridizations:
                condition = len(self.bonds) <= (cur_hybrid['steric_num'] -
                                                cur_hybrid['lone_pairs'])

                non_singles = 0
                for bond_type in ['double', 'triple', 'quadruple']:
                    if bond_type in cur_hybrid:
                        condition = condition and self.bond_counts[bond_type] == cur_hybrid[bond_type]
                        non_singles += cur_hybrid[bond_type]

                condition = condition and self.bond_counts['single'] == len(self.bonds) - non_singles

                if condition:
                    hybridization = cur_hybrid
                    break

        if hybridization is None:
            print('Hybridizations unknown for {}, taking a wild guess...'
                  .format(self.element.name))

            self.steric_num = len(self.bonds) + (self.hydrogens
                                                 if self.hydrogens is not None
                                                 else 0)
            self.lone_pairs = 0
        else:
            self.steric_num = hybridization['steric_num']
            self.lone_pairs = hybridization['lone_pairs']

        num_hydrogens = (self.steric_num - self.lone_pairs) - len(self.bonds)
        if self.hydrogens in [None, num_hydrogens]:
            for _ in range(num_hydrogens):
                self._add_hydrogen()
        else:
            raise RuntimeError('Wrong number of hydrogens (expected {}, got {})'
                               .format(num_hydrogens, self.hydrogens))

        if self.steric_num is None or self.lone_pairs is None:
            raise RuntimeError('Unknown steric num or lone pairs count')

        if self.tetrahedral not in [None, 0]:
            if not (self.steric_num == 4 and self.lone_pairs == 0):
                raise RuntimeError('Tetrahedral stereoisomer information '
                                   'is only relevant for tetrahedral atoms!')

            # Ok. The last three bonds in the tetrahedral geometry are listed
            # clockwise, so things will be clockwise by default. We have to pick
            # a default, and this will do. I THINK this is correct when it comes
            # to SMILES, but I'm not positive.
            #
            # If we get @@, this is still clockwise, so no change needs to be made.
            # In the case of @, however, (tetrahedral=1), we need to flip the ordering
            # of the last three bonds.
            if self.tetrahedral == 1:
                self.bonds[1:] = reversed(self.bonds[1:])

        if (self.steric_num == 3 and self.lone_pairs == 0 and
            self.bond_counts['double'] == 1 and
            self.bond_counts['single'] == 2):
            # Deal with trigonal planar stereochemistry

            if self.bond_counts['single-left'] > 1 or self.bond_counts['single-right'] > 1:
                raise RuntimeError('Inconsistent trigonal planar with double '
                                   'bond stereochemistry specification!')

            if self.bond_counts['single-left'] > 0 or self.bond_counts['single-right'] > 0:
                types = {bond.bond_type: bond.atom for bond in self.bonds}
                if 'single' in types:
                    new_direction = None
                    if 'single-left' in types:
                        new_direction = 'right'
                    else:
                        new_direction = 'left'

                    types['single-' + new_direction] = types['single']
                    del types['single']

                self.bonds = [Bond(t, types[t])
                              for t in ['single-left',
                                        'single-right',
                                        'double']]

        self.geometry = get_geometry(self.steric_num, self.lone_pairs)

    def _add_hydrogen(self):
        self.bonds.append(Bond('single', Atom(get_element(1), substituents=[Bond('single', self)])))
        self.bond_counts['single'] += 1

    def __repr__(self):
        return 'Atom({}, tetrahedral={}, hydrogens={}, charge={})'.format(self.element.name,
                                             repr(self.tetrahedral),
                                             repr(self.hydrogens),
                                             repr(self.charge))

    def replace_bond(self, val, atom):
        for bond in self.bonds:
            bond.replace_atom(val, atom)

    def find_bond_to(self, bond_val):
        for i, bond in enumerate(self.bonds):
            if bond.atom == bond_val:
                return i

    def rotate(self, angle, rv):
        self.geometry = [v.rotate(angle, rv) for v in self.geometry]

    def plane_normal(self):
        if self.steric_num == 3 and self.lone_pairs < 2:
            return (self.geometry[0] ** self.geometry[1]).normalize()
        elif self.steric_num == 4 and self.lone_pairs == 0:
            return (self.geometry[0] ** self.geometry[1]).normalize()
