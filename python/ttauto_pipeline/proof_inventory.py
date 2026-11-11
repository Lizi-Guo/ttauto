# SPDX-License-Identifier: GPL-3.0-or-later
"""Complete formal strata and the necessary homology representation in each."""
from ttauto_cert.lefschetz import partitions
from ttauto_cert.cover_orbits import disc_cover
from .config import keys, positive


def inventory(config):
    surface, selected = config['surface'], config['strata']
    cases = []
    if surface['kind'] == 'disc':
        n = surface['punctures']
        if selected == 'all':
            selected = []
            # Filling non-one-pronged internal punctures gives a containing
            # family on m<=n punctures, with the distinguished outer mark kept.
            for m in ([n] if surface['one_pronged'] else range(3, n+1)):
                for outer in range(1, m-1):
                    for orders in partitions(m-2-outer):
                        selected.append(dict(n=m, interior_prongs=[d+2 for d in orders], outer_prongs=outer))
        for row in selected:
            keys(row, {'n', 'interior_prongs', 'outer_prongs'}, 'disc stratum')
            m = positive(row.get('n'), 'stratum.n', 3)
            if m > n or (surface['one_pronged'] and m != n):
                raise ValueError('stratum outside requested puncture scope')
            outer = positive(row.get('outer_prongs'), 'outer_prongs')
            inner = row.get('interior_prongs')
            if not isinstance(inner, list) or any(type(p) is not int or p < 3 for p in inner):
                raise ValueError('interior singular prongs must be integers >=3')
            inner = sorted(inner)
            if -m+sum(p-2 for p in inner)+outer-2 != -4:
                raise ValueError('disc stratum violates Euler formula')
            cover = disc_cover(m, inner, outer)
            cases.append(dict(kind='disc', id=f'D{m}-I{inner}-O{outer}', n=m,
                              interior_prongs=inner, outer_prongs=outer, cover=cover,
                              rank=2*cover['genus'], homology='sphere_orientation_cover'))
        theorem = 'fill non-one-pronged internal punctures; enumerate all m=3..n and Euler partitions'
        if surface['one_pronged']:
            theorem = 'all Euler partitions with n one-pronged interior punctures'
    else:
        g = surface['genus']
        if selected == 'all':
            selected = [dict(type=kind, orders=list(p)) for kind, degree in (('H', 2*g-2), ('Q', 4*g-4))
                        for p in partitions(degree)]
        for row in selected:
            keys(row, {'type', 'orders'}, 'closed stratum')
            kind, orders = row.get('type'), row.get('orders')
            if kind not in ('H', 'Q') or not isinstance(orders, list) or not orders or any(type(d) is not int or d < 1 for d in orders):
                raise ValueError('closed stratum requires H/Q and positive zero orders')
            orders = sorted(orders)
            if sum(orders) != (2*g-2 if kind == 'H' else 4*g-4):
                raise ValueError('closed orders violate Euler formula')
            rank = 2*g if kind == 'H' else 2*g-2+sum(d % 2 for d in orders)
            cases.append(dict(kind='closed', id=f'g{g}-{kind}{orders}', type=kind, orders=orders,
                              genus=g, rank=rank, homology='base' if kind == 'H' else 'anti_invariant'))
        theorem = ('orientable-foliation homology; nonorientable case uses the integral anti-invariant '
                   'lattice of rank 2g-2+b, with invariant base homology of rank 2g')
    if len({r['id'] for r in cases}) != len(cases):
        raise ValueError('duplicate proof strata')
    return {'schema': 'ttauto-candidate-space-v1', 'surface': surface, 'cases': cases,
            'coverage': {'scope': 'all_requested_strata' if config['strata'] == 'all' else 'selected_strata',
                         'theorem_dependency': theorem,
                         'formal_types_may_be_empty': True, 'automaton_required': False}}
