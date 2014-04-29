from functools import reduce
import pyparsing as pp
import operator

import chemistry

def make_element_literal(element):
    return (pp.Literal(element.symbol)
              .setParseAction(lambda: [element])
              .setResultsName('element'))

elements = [(make_element_literal(element), element)
            for element in chemistry.elements]

organic_element = reduce(operator.xor,
                         (literal
                          for literal, element in elements
                          if element.smiles_organic))
element = reduce(operator.xor,
                 (literal for literal, element in elements))

tetrahedral = (pp.Regex('@{0,2}')
                 .setParseAction(lambda t: len(t[0]))
                 .setResultsName('tetrahedral'))

hydrogens = (pp.Empty().setParseAction(lambda: 0) ^
             pp.Literal('H').suppress() +
             (pp.Empty().setParseAction(lambda: 1) ^
              pp.Word(pp.nums).setParseAction(lambda t: int(t[0])))
            ).setResultsName('hydrogens')

bond_type = (pp.Empty().setParseAction(lambda: 1) ^
             pp.Literal('-').setParseAction(lambda: 1) ^
             pp.Literal('\\').setParseAction(lambda: (1, 'left')) ^
             pp.Literal('/').setParseAction(lambda: (1, 'right')) ^
             pp.Literal('=').setParseAction(lambda: 2) ^
             pp.Literal('#').setParseAction(lambda: 3) ^
             pp.Literal('$').setParseAction(lambda: 4)
            )

def parse_charge(tokens):
    op, num = tokens

    return {
        '-': operator.neg,
        '+': operator.pos
    }[op](int(num))

charge = (pp.ZeroOrMore('+').setParseAction(lambda t: len(t)) ^
          pp.ZeroOrMore('-').setParseAction(lambda t: -len(t)) ^
          ((pp.Literal('+') ^ pp.Literal('-')) +
           pp.Word(pp.nums)).setParseAction(parse_charge)
         ).setResultsName('charge')

isotope = pp.Optional(pp.Word(pp.nums).setParseAction(lambda t: int(t[0]))
                                      .setResultsName('isotope'))

atom = pp.Or([organic_element,
              pp.Literal('[').suppress() +
                  isotope + element + tetrahedral + hydrogens + charge +
              pp.Literal(']').suppress()],
             True
             )#.setParseAction(lambda t: Atom(**t.atom.asDict())).setResultsName('atom')

substituent = pp.Forward()

full_atom = (atom +
             pp.Or([pp.OneOrMore(substituent),
                   pp.Empty().setParseAction(lambda: [])],
                   True
             ).setResultsName('substituents')
            ).setParseAction(lambda t: t.asDict())

atom_chain = full_atom + pp.ZeroOrMore(bond_type + full_atom)

substituent_num = (pp.Word(pp.nums, exact=1) ^
                   (pp.Literal('%').suppress() + pp.Word(pp.nums))
                  ).setParseAction(lambda t: int(t[0]))

substituent <<= (pp.Group(pp.Literal('(').suppress() +
                          bond_type + atom_chain +
                          pp.Literal(')').suppress()) ^
                 substituent_num)

smiles = (pp.StringStart() +
          pp.Group(atom_chain) + pp.ZeroOrMore(pp.Literal('.').suppress() + pp.Group(atom_chain)) +
          pp.StringEnd())

def construct_substituent(smiles_substituent, numbered_substituents, i=0):
    cur_bond = smiles_substituent[i]

    return chemistry.Bond(cur_bond, construct_atom_graph(smiles_substituent, numbered_substituents, i+1, cur_bond))

def construct_atom_graph(smiles_atoms, numbered_substituents, i=0, behind_bond_type=None):
    cur_atom = smiles_atoms[i]

    processed_substituents = []

    if behind_bond_type:
        processed_substituents.append(chemistry.Bond(behind_bond_type, 'behind'))

    for substituent in cur_atom['substituents']:
        if not isinstance(substituent, int):
            substituent = construct_substituent(substituent, numbered_substituents)
            processed_substituents.append(substituent)
        else:
            processed_substituents.append(chemistry.Bond(1, substituent))

    if i < len(smiles_atoms) - 1:
        processed_substituents.append(construct_substituent(smiles_atoms, numbered_substituents, i+1))

    cur_atom['substituents'] = processed_substituents

    atom = chemistry.Atom(**cur_atom)

    for bond in cur_atom['substituents']:
        if isinstance(bond.atom, int):
            numbered_substituents.setdefault(bond.atom, set()).add(atom)
        elif not isinstance(bond.atom, str):
            bond.atom.replace_bond('behind', atom)

    return atom

if __name__ == '__main__':
    tests = [
        'N#N',
        'CN=C=O',
        'OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1',
        'CC[C@H](O1)CC[C@@]12CCCO2',
        'CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2',
        '[2H]C(Cl)(Cl)Cl'
    ]

    from pprint import pprint
    for test in tests:
        print(test)
        print(repr(smiles.parseString(test)))
        print('---')
