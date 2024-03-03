load("constraint_network.sage")

from itertools import combinations as icomb

def all_different(C, *names):
    for (i, j) in icomb(names, 2):
        C[i, j] = lambda a, b: a != b

_zebra_names = ['Parliament', 'Kools', 'Chesterfield', 'OldGold', 'Lucky',
                'Zebra', 'Horse', 'Fox', 'Dog', 'Snails',
                'Japanese', 'Ukrainian', 'English', 'Spaniard', 'Norwegian',
                'Ivory', 'Yellow', 'Red', 'Green', 'Blue',
                'Water', 'Tea', 'Coffee', 'Milk', 'Orange']

_zebra_domain_dict = {name:[1..5] for name in _zebra_names}
_zebra_domain_dict['Milk'] = [3]
_zebra_domain_dict['Norwegian'] = [1]

_zebra_constraints = dict()
for row in [_zebra_names[5*i:5*i+5] for i in range(5)]:
    all_different(_zebra_constraints, *row)

_zebra_constraints['English', 'Red'] = lambda x, y: x == y
_zebra_constraints['Spaniard', 'Dog'] = lambda x, y: x == y
_zebra_constraints['Coffee', 'Green'] = lambda x, y: x == y
_zebra_constraints['Ukrainian', 'Tea'] = lambda x, y: x == y
_zebra_constraints['Green', 'Ivory'] = lambda x, y: x - 1 == y
_zebra_constraints['OldGold', 'Snails'] = lambda x, y: x == y
_zebra_constraints['Kools', 'Yellow'] = lambda x, y: x == y
_zebra_constraints['Chesterfield', 'Fox'] = lambda x, y: abs(x - y) == 1
_zebra_constraints['Yellow', 'Horse'] = lambda x, y: abs(x - y) == 1
_zebra_constraints['Norwegian', 'Blue'] = lambda x, y: abs(x - y) == 1
_zebra_constraints['Lucky', 'Orange'] = lambda x, y: x == y
_zebra_constraints['Japanese', 'Parliament'] = lambda x, y: x == y

def _C():
    return ConstraintNetwork(_zebra_domain_dict.copy(), _zebra_constraints.copy())

C = _C()

def print_zebra_domains(C):
    rows = [_zebra_names[5*i:5*i+5] for i in range(5)]
    max_domain_size = max(len(v) for v in C.domains.values())
    field_width = 3 * max_domain_size + 1
    for row in rows:
        replacement_field = '{!s:^' + str(field_width) + '}'
        format_string = replacement_field * 5
        print(format_string.format(*[C.domains[name] for name in row]))
            
