"""Exact necessary conditions and finite lower-bound proofs on orientating covers.

This is a lower-bound research program, not a realizability decision procedure.
All elimination steps use integer/Fraction arithmetic and necessary conditions.
Coefficients are descending. No complex-root numerical filter is used.
"""
from fractions import Fraction
from itertools import product
from math import gcd
from .exact_math import count_real_roots, compare_largest_real_roots, compare_largest_real_root

SOURCE = 'https://arxiv.org/pdf/1004.5344'

def _polynomial(p, reciprocal=False):
    if len(p)<2 or p[0]!=1 or any(type(c) is not int for c in p):
        raise ValueError('a monic polynomial with integer coefficients is required')
    if reciprocal and (len(p)%2!=1 or list(p)!=list(reversed(p))):
        raise ValueError('an even-degree monic reciprocal polynomial is required')

def _horizon(stop):
    if type(stop) is not int or stop<1:raise ValueError('trace horizon must be a positive integer')

def _groups(groups):
    if any(len(block)!=2 or any(type(x) is not int for x in block) or
           block[0]<4 or block[0]%2 or block[1]<1 for block in groups):
        raise ValueError('singular blocks must be (even prongs >= 4, positive multiplicity)')

def traces(p, stop):
    _polynomial(p);_horizon(stop)
    d=len(p)-1
    s=[d]
    for k in range(1,stop+1):
        s.append(-sum(p[j]*s[k-j] for j in range(1,min(k,d+1)))
                 -(k*p[k] if k<=d else 0))
    return s

def reciprocal_superset(g,bound,stop=24):
    """Include every monic reciprocal P of degree 2g whose house is <= bound
    and which has a positive real root in (1,bound]. In particular this includes
    the characteristic polynomial of every selected positive pA lift below
    the bound; polynomials with no positive expanding real root are omitted.

    Reciprocal pairs imply |Tr(A^k)| <= g*(bound^k+bound^-k).
    Newton identities give finite coefficient ranges recursively.
    Extra trace/real-root conditions remove only impossible polynomials.
    """
    if isinstance(bound,float):raise TypeError('use an exact rational house bound')
    bound=Fraction(bound)
    _horizon(stop)
    if type(g) is not int or g<1 or bound<=1 or stop<g:
        raise ValueError('invalid genus, bound or trace horizon')
    B=[None]+[g*(bound**k+bound**(-k)) for k in range(1,stop+1)]
    def extend(coeff, s):
        k=len(coeff)
        if k>g:
            p=coeff+coeff[-2::-1]
            tr=traces(p,stop)
            if any(abs(tr[i])>B[i] for i in range(1,stop+1)): return
            if not count_real_roots(p,1,bound): return
            if count_real_roots(p,bound,None) or count_real_roots(p,None,-bound,include_right=False): return
            yield p
            return
        t=-sum(coeff[j]*s[k-j] for j in range(1,k))
        lo=(Fraction(t)-B[k])/k
        hi=(Fraction(t)+B[k])/k
        for a in range(-((-lo.numerator)//lo.denominator),hi.numerator//hi.denominator+1):
            yield from extend(coeff+[a],s+[t-k*a])
    yield from extend([1],[2*g])

def partitions(n,lo=1):
    if type(n) is not int or type(lo) is not int or n<0 or lo<1:
        raise ValueError('invalid integer partition arguments')
    if not n:
        yield ()
    for a in range(lo,n+1):
        for tail in partitions(n-a,a):yield (a,)+tail

def singular_patterns(groups,negative,stop,budget=None):
    """groups are (prongs,number) independently permutable blocks.

    Keep all cyclic separatrix shifts with the right orientation parity.
    A block may be split to enforce marked singularities being fixed.
    """
    _groups(groups);_horizon(stop)
    if type(negative) is not bool:raise ValueError('negative must be Boolean')
    def options(prongs,num):
        for part in partitions(num):
            for shifts in product(*[range(l%2 if negative else 0,prongs,2) for l in part]):
                if budget is not None:budget.tick()
                contribution=[0]*(stop+1)
                fixed=[0]*(stop+1)
                for l,shift in zip(part,shifts):
                    for m in range(l,stop+1,l):
                        fixed[m]+=l
                        contribution[m]+=l*((1-prongs) if shift*(m//l)%prongs==0 else 1)
                yield contribution,fixed,tuple(zip(part,shifts))
    def visit(i,selected):
        if i==len(groups):
            if budget is not None:budget.tick()
            yield ([sum(x[0][m] for x in selected) for m in range(stop+1)],
                   [sum(x[1][m] for x in selected) for m in range(stop+1)],
                   [x[2] for x in selected])
            return
        for option in options(*groups[i]):yield from visit(i+1,selected+[option])
    yield from visit(0,[])

def compatible(p,groups,negative=False,regular_fixed_min=0,stop=24,*,budget=None):
    _polynomial(p,reciprocal=True);_groups(groups);_horizon(stop)
    # Historical sphere-cover API: p is the positive-lift polynomial and the
    # other lift acts as -A. General covers must use actual homological traces.
    trace=traces(p,stop)[1:]
    if negative:
        trace=[(-1 if k%2 else 1)*t for k,t in enumerate(trace,1)]
    return compatible_traces(trace,len(p)-1,groups,negative,regular_fixed_min,budget=budget)


def compatible_traces(actual_traces,rank,groups,negative=False,regular_fixed_min=0,*,budget=None):
    """Necessary orbit conditions with trace and foliation sign independent.

    actual_traces[k-1] is Tr(B**k) of the ACTUAL lift B. `negative` only
    controls the local foliation orientation and singularity rotations; it
    never changes B or its traces. On a general cover B may be J*A, not -A.
    """
    _groups(groups);_horizon(len(actual_traces))
    if type(rank) is not int or rank<2 or rank%2 or any(type(t) is not int for t in actual_traces):
        raise ValueError('even homology rank and exact integer traces required')
    if type(negative) is not bool:
        raise ValueError('negative must be Boolean')
    if sum((prongs-2)*number for prongs,number in groups)!=2*rank-4:
        raise ValueError('singularities do not satisfy the Euler formula for the cover genus')
    if not callable(regular_fixed_min) and (type(regular_fixed_min) is not int or regular_fixed_min<0):
        raise ValueError('regular fixed point lower bound must be nonnegative')
    stop=len(actual_traces)
    trace=[rank]+list(actual_traces)
    feasible=[]
    for singular, fixed, description in singular_patterns(groups,negative,stop,budget):
        orbits=[0]*(stop+1);regular=[0]*(stop+1)
        for m in range(1,stop+1):
            neg=negative and m%2
            lef=2-trace[m]
            regular[m]=(lef-singular[m]) if neg else (singular[m]-lef)
            required=regular_fixed_min(m) if callable(regular_fixed_min) else regular_fixed_min
            if type(required) is not int or required<0:raise ValueError('invalid fixed point lower bound')
            if regular[m]<required:break
            points=regular[m]-sum(d*orbits[d] for d in range(1,m) if m%d==0)
            if points<0 or points%m:break
            orbits[m]=points//m
        else:feasible.append({'singular_orbits':description,'regular_orbits':orbits[1:], 'singular_fixed':fixed[1:]})
    return feasible

def rotation_symmetry(p):
    """If P=Q(x^d), d>1, a nonzero root has equal-modulus companions.

    This excludes the STRICT dominant homological eigenvalue of an orientable
    pseudo-Anosov. It must not be used to classify arbitrary PF graph matrices.
    """
    _polynomial(p)
    d=0
    for exponent,coefficient in enumerate(reversed(p)):
        if coefficient:d=gcd(d,exponent)
    return d

def equivariant_regular_compatible(positive,negative,branch_regular_count,branch_fixed_min=0,marked_free_pair=False):
    """Necessary commuting-deck-involution constraints, on regular orbits.

    Remove the same branch-point permutation from both lifts. Off the branch
    locus tau acts freely. For odd d, paired phi d-orbits become one psi 2d-orbit;
    a tau-invariant phi 2d-orbit becomes two psi d-orbits. At lengths divisible by
    four the two orbit counts agree. At 2d (d odd), residual counts are equal,
    even and nonnegative. This uses psi=tau*phi, tau^2=1, [tau,phi]=1.
    It is necessary only: singular-point/deck compatibility is still relaxed.
    """
    stop=len(positive['regular_orbits'])
    if stop<2 or len(negative['regular_orbits'])!=stop or any(type(x) is not int or x<0
            for x in positive['regular_orbits']+negative['regular_orbits']):
        raise ValueError('nonnegative orbit counts with equal horizons >= 2 are required')
    if type(branch_regular_count) is not int or branch_regular_count<0 or type(branch_fixed_min) is not int or not 0<=branch_fixed_min<=branch_regular_count:
        raise ValueError('invalid branch regular-point counts')
    for part in partitions(branch_regular_count):
        if part.count(1)<branch_fixed_min:continue
        branch=[0]*(stop+1)
        for d in part:
            if d<=stop:branch[d]+=1
        a=[0]+[positive['regular_orbits'][d-1]-branch[d] for d in range(1,stop+1)]
        b=[0]+[negative['regular_orbits'][d-1]-branch[d] for d in range(1,stop+1)]
        if any(x<0 for x in a+b):continue
        # A distinguished fixed regular point downstairs lifts to a tau-pair.
        # It is either a paired phi fixed orbit, or a paired psi fixed orbit.
        if marked_free_pair and a[1]+b[1]<2:continue
        good=True
        for d in range(1,stop+1):
            if d%2:
                if a[d]%2 or b[d]%2:good=False;break
                if 2*d<=stop:
                    x=a[2*d]-b[d]//2;y=b[2*d]-a[d]//2
                    if x<0 or x!=y or x%2:good=False;break
            elif d%4==0 and a[d]!=b[d]:good=False;break
        if good:return True
    return False

def branch_singular_compatible(positive,negative,groups,indices):
    # At a branch singularity the deck involution rotates by p/2 prongs.
    # For a length-l orbit: (tau*phi)^l=tau^l*phi^l.
    _groups(groups)
    for i in indices:
        if type(i) is not int or not 0<=i<len(groups) or groups[i][0]%4!=2:
            raise ValueError('a branch singularity must have prongs congruent to 2 modulo 4')
        prongs=groups[i][0]
        expected=sorted((l,(shift+(l%2)*(prongs//2))%prongs)
                        for l,shift in positive['singular_orbits'][i])
        actual=sorted(tuple(x) for x in negative['singular_orbits'][i])
        if expected!=actual:return False
    return True

def self_test():
    def check(condition):
        if not condition:raise AssertionError('Lefschetz arithmetic self-check failed')
    check(traces([1,-3,1],6)==[2,3,7,18,47,123,322])
    check([1,-2,1,-2,1] in list(reciprocal_superset(2,Fraction(1884,1000))))
    check([1,0,-1,-4,-1,0,1] in list(reciprocal_superset(3,Fraction(1840,1000))))
    check(rotation_symmetry([1,0,0,-4,0,0,1])==3)
    check(rotation_symmetry([1,-2,1,-2,1])==1)
    # A 12-prong point may rotate by four prongs, with order three (not six).
    patterns=list(singular_patterns([(12,1)],False,6))
    check(any(x[0][1:]==[1,1,-11,1,1,-11] for x in patterns))
    # Positive/negative torus Anosov Lefschetz counts must survive.
    check(compatible([1,-3,1],[],False) and compatible([1,-3,1],[],True))
    check(any(equivariant_regular_compatible(a,b,4) for a in compatible([1,-3,1],[],False)
               for b in compatible([1,-3,1],[],True)))
    # Independently test the deck-orbit identities on random permutations
    # commuting with the free involution (0 1)(2 3)... .
    import random
    rng=random.Random(602025)
    def cycles(perm):
        counts=[0]*24;seen=set()
        for v in range(len(perm)):
            if v in seen:continue
            length=0
            while v not in seen:seen.add(v);length+=1;v=perm[v]
            counts[length-1]+=1
        return {'regular_orbits':counts}
    for pairs in range(1,9):
        for repeat in range(20):
            order=list(range(pairs));rng.shuffle(order)
            flip=[rng.randrange(2) for _ in order]
            phi=[2*order[i//2]+((i%2)^flip[i//2]) for i in range(2*pairs)]
            psi=[j^1 for j in phi]
            check(equivariant_regular_compatible(cycles(phi),cycles(psi),0))

def run(include_star=False,include_g4=False,progress=None):
    """Prove selected n=6 lower bounds or report the remaining obstruction set.

    An empty obstruction set establishes a universal stratum lower bound under
    the cover assumptions recorded in each row. A nonempty set says nothing
    about realization. Upper witnesses must be certified separately.
    """
    # Upper rational enclosures are conservative, exact candidate comparisons follow.
    specifications=[
        ('s1',2,Fraction(1884,1000),[1,-2,1,-2,1],[(4,2)],0,6,0,False,[]),
        ('s2',3,Fraction(1840,1000),[1,0,-1,-4,-1,0,1],[(6,1),(6,1)],0,6,0,False,[0,1]),
        ('s3',2,Fraction(1884,1000),[1,-2,1,-2,1],[(4,2)],0,6,0,True,[]),
        ('s6',3,Fraction(1884,1000),[1,-2,1,-2,1],[(6,1),(4,2)],1,7,1,False,[0]),
    ]
    if include_star:
        specifications += [
            ('s4',3,Fraction(2082,1000),[1,-1,-2,-1,1],[(10,1)],1,7,1,False,[0]),
            ('s5',3,Fraction(2082,1000),[1,-1,-2,-1,1],[(6,2)],lambda m:2 if m%2==0 else 0,6,0,True,[0]),
        ]
    if include_g4:
        specifications += [
            ('s7',4,Fraction(2172,1000),[1,-2,1,-4,4,-4,1,-2,1],[(6,3)],1,7,1,False,[0]),
        ]
    results={}
    for name,g,bound,target,groups,marked,branch_regular,branch_fixed,marked_pair,branch_singular in specifications:
        if compare_largest_real_root(target,bound)>=0:
            raise ValueError('house enclosure must strictly exceed the target root')
        candidates=[p for p in reciprocal_superset(g,bound) if compare_largest_real_roots(p,target)<0]
        survivors=[];rejected=[]
        for p in candidates:
            d=rotation_symmetry(p)
            if d>1:
                rejected.append({'polynomial':p,'reason':'root_rotation_symmetry','order':d})
                continue
            opposite=[a*((-1)**i) for i,a in enumerate(p)]
            if compare_largest_real_roots(p,opposite)<=0:
                rejected.append({'polynomial':p,'reason':'negative_root_at_least_as_large'})
                continue
            pos=compatible(p,groups,False,marked)
            neg=compatible(p,groups,True,marked)
            equivariant=any(branch_singular_compatible(a,b,groups,branch_singular) and
                            equivariant_regular_compatible(a,b,branch_regular,branch_fixed,marked_pair)
                            for a in pos for b in neg)
            if pos and neg and equivariant:survivors.append({'polynomial':p,'positive':pos,'negative':neg})
            else:rejected.append({'polynomial':p,'positive_feasible':bool(pos),'negative_feasible':bool(neg),'deck_equivariant_feasible':equivariant})
        results[name]={'schema':'ttauto-lefschetz-lower-bound-v1','n':6,'stratum_index':int(name[1:])-1,
            'source':SOURCE,'source_locator':'Sections 2.4-2.6; Appendix A.2 Table 6',
            'assumptions':['orientation-preserving pseudo-Anosov on a disc',
                'six interior punctures each have one prong','outer marked point fixed',
                'canonical orientating cover and positive-leading-eigenvalue lift',
                'closed-cover homology is integral symplectic with a unique dominant eigenvalue'],
            'bound':str(bound),'target':target,'genus':g,'singular_groups':groups,'trace_horizon':24,
            'branch_regular_count':branch_regular,'branch_fixed_min':branch_fixed,
            'marked_free_pair':marked_pair,'branch_singular_groups':branch_singular,
            'regular_fixed_min_first_24':[marked(m) if callable(marked) else marked for m in range(1,25)],
            'candidate_superset':candidates,'rejected':rejected,'survivors':survivors,
            'lower_bound_established':not survivors,'realizability_claimed':False}
        if progress is not None:progress(name,len(candidates),len(survivors),[x['polynomial'] for x in survivors])
    return results
