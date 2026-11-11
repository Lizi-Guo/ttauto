"""Enumerate ALL embedded multigon-tree types, independently of fold reachability.

Completeness is within the specified tree/multigon model. It does not assert
that every pA has an invariant representative in that model. No Rauzy method
is used here. Integer tuples match the C++ ttcode-v1 fields exactly.
"""
from itertools import permutations, product
from functools import lru_cache
import heapq
from math import comb, factorial, prod


def partitions(n, least=1):
    """Nondecreasing integer partitions in the historic order, without recursion."""
    if n < 0 or least < 1:
        return
    prefix = []
    stack = [(n, iter(range(least, n+1)))]
    while stack:
        remaining, choices = stack[-1]
        if remaining == 0:
            yield tuple(prefix)
            stack.pop()
            if prefix: prefix.pop()
            continue
        try:
            a = next(choices)
        except StopIteration:
            stack.pop()
            if prefix: prefix.pop()
            continue
        prefix.append(a)
        stack.append((remaining-a, iter(range(a, remaining-a+1))))


def strata(n):
    if n < 3:
        raise ValueError("at least three interior punctures are required")
    for defect in range(n-2):
        for part in partitions(defect):
            yield {"n":n, "interior_prongs":tuple(a+2 for a in part),
                   "outer_prongs":n-2-defect}


def positive_compositions(n, k):
    if k < 1 or n < k:
        return
    prefix = []
    stack = [(n, k, iter(range(1, n-k+2)))]
    while stack:
        remaining, slots, choices = stack[-1]
        if slots == 1:
            yield tuple(prefix)+(remaining,)
            stack.pop()
            if prefix: prefix.pop()
            continue
        try:
            a = next(choices)
        except StopIteration:
            stack.pop()
            if prefix: prefix.pop()
            continue
        prefix.append(a)
        stack.append((remaining-a, slots-1, iter(range(1, remaining-a-(slots-1)+2))))


def prufer_tree(word, size):
    degree = [1]*size
    for a in word:
        degree[a] += 1
    leaves = [a for a in range(size) if degree[a] == 1]
    heapq.heapify(leaves)
    adjacent = [[] for _ in range(size)]
    for a in word:
        leaf=heapq.heappop(leaves)
        adjacent[a].append(leaf)
        adjacent[leaf].append(a)
        degree[a] -= 1
        if degree[a] == 1:
            heapq.heappush(leaves,a)
    a,b=leaves
    adjacent[a].append(b)
    adjacent[b].append(a)
    return adjacent


@lru_cache(None)
def local_orders(neighbors, prongs):
    """Every linear order at each prong, modulo rotation of whole prongs."""
    found=set()
    for perm in permutations(neighbors):
        for lengths in positive_compositions(len(neighbors),prongs):
            groups=[]
            offset=0
            for length in lengths:
                groups.append(perm[offset:offset+length])
                offset+=length
            groups=tuple(groups)
            found.add(min(groups[j:]+groups[:j] for j in range(prongs)))
    return tuple(sorted(found))


def canonical_code(orders):
    """Minimum DFS code among univalent monogon roots; preserve orientation."""
    codes=[]
    for root,row in enumerate(orders):
        if len(row) != 1 or len(row[0]) != 1:
            continue
        code=[0,1,0,0,1]

        def incidences_from_parent(vertex, parent):
            groups=orders[vertex]
            incoming=next((p,e) for p,g in enumerate(groups)
                          for e,v in enumerate(g) if v == parent)
            incidences=[(p,e,v) for p,g in enumerate(groups) for e,v in enumerate(g)]
            start=next(j for j,(p,e,_) in enumerate(incidences) if (p,e)==incoming)
            for p,e,target in incidences[start:]+incidences[:start]:
                yield ((p-incoming[0])%len(groups),len(groups),0,e,len(groups[p])), target

        vertex = row[0][0]
        stack = [(vertex, root, iter(incidences_from_parent(vertex, root)))]
        while stack:
            vertex, parent, choices = stack[-1]
            try:
                block, target = next(choices)
            except StopIteration:
                stack.pop()
                continue
            code.extend(block)
            if target != parent:
                stack.append((target, vertex, iter(incidences_from_parent(target, vertex))))
        codes.append(tuple(code))
    if not codes:
        raise ValueError("tree has no uncusped monogon root")
    return min(codes)


def _validate_stratum(n, interior_prongs):
    k = tuple(interior_prongs)
    if n < 3 or any(a < 3 for a in k) or n-2-sum(a-2 for a in k) < 1:
        raise ValueError("invalid paper stratum")
    return k


def degree_constrained_prufer_words(prongs):
    """Every Prüfer word whose vertex degrees are at least ``prongs``.

    Vertex i occurs deg(i)-1 times in its Prüfer word. A prefix is discarded
    exactly when its unmet occurrence quotas exceed the remaining positions.
    No unlabeled-tree or folding-graph information enters this generator.
    """
    prongs = tuple(prongs)
    size = len(prongs)
    if size < 2 or any(p < 1 for p in prongs):
        raise ValueError("need at least two vertices of positive minimum degree")
    deficits = [p-1 for p in prongs]
    missing = sum(deficits)
    length = size-2
    if missing > length:
        return
    word = []

    next_letters = [0]
    used_counts = []
    while next_letters:
        remaining = length-len(word)
        if remaining == 0:
            yield tuple(word)
            next_letters.pop()
            if word:
                letter, used = word.pop(), used_counts.pop()
                deficits[letter] += used
                missing += used
            continue
        if next_letters[-1] == size:
            next_letters.pop()
            if word:
                letter, used = word.pop(), used_counts.pop()
                deficits[letter] += used
                missing += used
            continue
        letter = next_letters[-1]
        next_letters[-1] += 1
        used = int(deficits[letter] > 0)
        if missing-used > remaining-1:
            continue
        deficits[letter] -= used
        missing -= used
        word.append(letter)
        used_counts.append(used)
        next_letters.append(0)


def enumerate_tree_types_prufer(n, interior_prongs, *, degree_constrained=False):
    """Exhaustive labeled trees + local orders, then exact isomorphism quotient."""
    k=_validate_stratum(n,interior_prongs)
    size=n+len(k)
    if size>7 and not degree_constrained:
        raise ValueError("this audited Prufer implementation is limited to seven multigons")
    prongs=(1,)*n+k
    codes=set()
    counts={"prufer_words":0,"degree_admissible_trees":0,"decorations":0}
    words = (degree_constrained_prufer_words(prongs) if degree_constrained
             else product(range(size),repeat=size-2))
    for word in words:
        counts["prufer_words"]+=1
        degrees=[1]*size
        for a in word:
            degrees[a]+=1
        if any(a<b for a,b in zip(degrees,prongs)):
            continue
        counts["degree_admissible_trees"]+=1
        neighbors=prufer_tree(word,size)
        choices=[local_orders(tuple(sorted(row)),p) for row,p in zip(neighbors,prongs)]
        for orders in product(*choices):
            counts["decorations"]+=1
            codes.add(canonical_code(orders))
    counts["types"]=len(codes)
    if degree_constrained:
        counts["algorithm"]="degree_constrained_prufer"
        counts["unconstrained_word_count"]=size**(size-2)
    return sorted(codes),counts


def enumerate_tree_types_planted(n, interior_prongs):
    """Enumerate rooted decorated trees bottom-up, without vertex labels.

    A planted subtree has a distinguished incoming edge. Rotate the prongs so
    that this edge lies in prong zero. That prong consists of an ordered list
    before the edge, the incoming edge, and an ordered list after it. Each of
    the other p-1 prongs is a nonempty ordered list of child subtrees. This is
    a unique recursive description. Finally attach a univalent monogon root
    and quotient its choices using the original exact DFS canonical code.

    Tables are indexed by the multiset of vertex prongs, in increasing total
    size. Tree nodes contain integer child IDs, not recursively nested tuples.
    Neither generation, hashing nor conversion uses a tree-depth call stack.
    The algorithm never consults transition matrices or graph reachability.
    """
    k=_validate_stratum(n,interior_prongs)
    colors=(1,)+tuple(sorted(set(k)))
    inventory=(n-1,)+tuple(k.count(p) for p in colors[1:])
    generated={"planted_subtrees":0,"ordered_forests":0}

    def subtract(a,b):
        return tuple(x-y for x,y in zip(a,b))

    def forest_defect(counts):
        # For c planted components: sum(2-prongs) >= c. In particular a
        # nonempty forest with nonpositive defect is impossible.
        return sum(count*(2-p) for count,p in zip(counts,colors))

    def subinventories(counts):
        return (sub for sub in product(*(range(c+1) for c in counts))
                if any(sub) and forest_defect(sub) >= 1)

    @lru_cache(None)
    def cut_lengths(length, prongs):
        if length < prongs-1:
            return ()
        if prongs == 1:
            return tuple((before,length-before) for before in range(length+1))
        return tuple((before,after)+rest
                     for before in range(length-prongs+2)
                     for after in range(length-before-prongs+2)
                     for rest in positive_compositions(length-before-after,prongs-1))

    empty=(0,)*len(inventory)
    tree_table={empty:()}
    forest_table={empty:((),)}
    nodes=[]
    node_ids={}
    inventories=sorted(product(*(range(c+1) for c in inventory)),key=lambda c:(sum(c),c))
    for counts in inventories:
        if counts==empty:
            continue
        if forest_defect(counts)<1:
            tree_table[counts]=()
            forest_table[counts]=()
            continue
        found=[]
        for color,p in enumerate(colors):
            if counts[color] == 0:
                continue
            remainder=list(counts)
            remainder[color]-=1
            if sum(remainder) < p-1:
                continue
            for children in forest_table[tuple(remainder)]:
                for lengths in cut_lengths(len(children),p):
                    groups=[]
                    offset=0
                    for length in lengths:
                        groups.append(children[offset:offset+length])
                        offset+=length
                    node=(p,groups[0],groups[1],tuple(groups[2:]))
                    if node in node_ids:
                        raise AssertionError('duplicate planted-tree construction')
                    node_ids[node]=len(nodes)
                    found.append(len(nodes))
                    nodes.append(node)
        generated["planted_subtrees"]+=len(found)
        tree_table[counts]=tuple(found)
        # A nonempty forest is uniquely its first tree followed by a forest.
        # trees(counts) is already available; every nonempty tail is smaller.
        found=[]
        for first_inventory in subinventories(counts):
            remainder=subtract(counts,first_inventory)
            tails=forest_table[remainder]
            for tree in tree_table[first_inventory]:
                for tail in tails:
                    found.append((tree,)+tail)
        generated['ordered_forests']+=len(found)
        forest_table[counts]=tuple(found)

    def to_orders(tree):
        orders=[((1,),),None]
        stack=[(tree,1,0)]
        while stack:
            node,here,parent=stack.pop()
            p,before,after,other=nodes[node]
            groups=[]
            for group in (before+(None,)+after,)+other:
                neighbors=[]
                for child in group:
                    if child is None:
                        neighbors.append(parent)
                    else:
                        child_vertex=len(orders)
                        orders.append(None)
                        neighbors.append(child_vertex)
                        stack.append((child,child_vertex,here))
                groups.append(tuple(neighbors))
            if len(groups)!=p:
                raise AssertionError("planted-tree prong count mismatch")
            orders[here]=tuple(groups)
        return tuple(orders)

    rooted=tree_table[inventory]
    # Independent coefficient check from formal Lagrange inversion:
    # T = sum_p x_p T^(p-1)/(1-T)^(p+1).
    total=sum(inventory)
    deficit=forest_defect(inventory)
    numerator=factorial(total-1)*comb(3*total-2,deficit-1)
    denominator=prod(factorial(c) for c in inventory)
    expected,remainder=divmod(numerator,denominator)
    if remainder or len(rooted) != expected or len(set(rooted)) != len(rooted):
        raise AssertionError("planted enumeration disagrees with the exact coefficient formula")
    codes=sorted({canonical_code(to_orders(tree)) for tree in rooted})
    counts=dict(generated,algorithm="planted_tree_bottom_up",
                rooted_types=len(rooted),types=len(codes),
                rooted_count_formula=expected,
                subtree_inventory_states=len(tree_table),
                forest_inventory_states=len(forest_table),
                tree_depth_uses_python_recursion=False)
    return codes,counts


def enumerate_tree_types(n, interior_prongs):
    """Preserve the n<=5 certificate API; use bottom-up enumeration above it."""
    if n <= 5:
        return enumerate_tree_types_prufer(n,interior_prongs)
    return enumerate_tree_types_planted(n,interior_prongs)
