from itertools import accumulate
from enum import Enum
from pprint import pformat
import networkx as nx
import matplotlib


class Color(Enum):
    BLACK = 0
    WHITE = 1
    BOTH = 2  # todo: decide if I really need a BOTH Color

    def __add__(self, other):
        if self == other:
            return self
        if self == BOTH or other == BOTH:
            return BOTH
        if self == WHITE and other == BLACK:
            return BOTH
        if self == BLACK and other == WHITE:
            return BOTH


BLACK, WHITE, BOTH = Color.BLACK, Color.WHITE, Color.BOTH

constants = [
    ((None, None), (BLACK, None)),
    ((BLACK, None), (None, None)),
    ((None, BLACK), (None, None)),
    ((None, None), (None, BLACK)),
    ((None, None), (WHITE, None)),
    ((WHITE, None), (None, None)),
    ((None, WHITE), (None, None)),
    ((None, None), (None, WHITE)),
]


def merge(a, b):
    # binary operation that is commutative, associative and idempotent
    # a and b are elements of the algebra. They can be constants, terms or atoms.
    if a == b:
        return a
    if a == 0:
        return b
    if b == 0:
        return a
    result = []
    for i in range(len(a)):
        row = []
        for j in range(len(a[0])):
            ea, eb = a[i][j], b[i][j]
            if ea is None:
                row.append(eb)
            elif eb is None:
                row.append(ea)
            else:
                row.append(ea + eb)
        result.append(tuple(row))
    return tuple(result)


def merge_all(element_set):
    return tuple(accumulate(element_set, merge))[-1]


def less(graph, atoms_list, a, b):
    # a < b iff (for all phi in atoms_list, (phi not -> a) or (phi -> b))
    for phi in atoms_list:
        neighbors = graph[phi]
        if not (
                (a not in neighbors) or (b in neighbors)
        ):
            return False
    return True


def dual(node):
    # todo
    return None, node


def gla(graph, atoms_list, x):
    # Graph Lower Segment Atoms
    # the following identity holds: gla(merge(x,y)) === gla(x) | gla(y)
    in_neighs = list(zip(*graph.in_edges(x)))[0]
    return set(in_neighs) & set(atoms_list)


def crossing():
    # todo
    pass


# each term is a merge of constants
# each constant is a merge of atoms


T1p = {constants[2 - 1], constants[7 - 1], constants[1 - 1], constants[8 - 1]}
T2p = {constants[6 - 1], constants[3 - 1], constants[5 - 1], constants[4 - 1]}
T1n = {constants[2 - 1], constants[7 - 1], constants[5 - 1], constants[4 - 1]}
T2n = {constants[6 - 1], constants[3 - 1], constants[5 - 1], constants[8 - 1]}
T3n = {constants[6 - 1], constants[7 - 1], constants[5 - 1], constants[4 - 1]}

positive_class_def, negative_class_def = [T1p, T2p], [T1n, T2n, T3n]
# positive_class, negative_class = [merge_all(x) for x in positive_class_def], [merge_all(x) for x in negative_class_def]
training_examples_def = positive_class_def + negative_class_def

# v is a constant (the vertical bar in this example) which we want to describe to our algebra

# training_set R definition:
#   for T in positive_class: v < T
#   for T in negative_class: not (v < T)

# the graph G contains as nodes all T in both classes as well as v.
# (a -> b) \implies (a < b), i. e., a subset of the partial order < is represented by G. todo
# if all c_i of a term T also are part of a term S we add the edge (T -> S)
# we always use edges if any of the elements involved are atoms: (\phi -> b) \iff (phi < b) todo


if __name__ == '__main__':

    # the graph G contains as nodes all T in both classes as well as v.
    G_start = nx.DiGraph()
    master_atoms = [0]
    #   for every term T which is DEFINED as the sum os constants c_i we add the edges (c_i -> T)
    G_start.add_edge(0, 'v')
    for c in constants:
        G_start.add_edge(0, c)
    for T_def in training_examples_def:
        for c in T_def:
            G_start.add_edge(c, merge_all(T_def))

    # if all c_i of a term T also are part of a term S we add the edge (T -> S)
    for T_def in training_examples_def:
        for S_def in training_examples_def:
            if T_def is S_def:
                continue
            bond = True
            for c in T_def:
                if c not in S_def:
                    bond = False
                    break
            if bond:
                T, S = merge_all(T_def), merge_all(S_def)
                G_start.add_edge(T, S)

    # the graph G of our algebra M which will evolve during the learning process
    G_master = nx.algorithms.transitive_closure(G_start, reflexive=None)

    # # # # # DUALS
    # Dual of constant or term --> constant
    # Dual of atoms --> dual-of-atom

    # M* has constants, dual-of-atoms and atoms, but does not have terms.
    # Atoms in M* are not the dual of anything in M.

    # axiom: (a -> b) \implies ([b] -> [a])
    G_dual = G_master.reverse()
    mapping = dict(zip(G_dual.nodes, map(dual, G_dual.nodes)))
    G_dual = nx.relabel_nodes(G_dual, mapping)

    dual_atoms = ['0*']
    # adds the '0*' atom
    for T_def in training_examples_def:
        T = merge_all(T_def)
        G_dual.add_edge('O*', dual(T))

    # In M* we add additional edges for the positive order relations of R such as v < T1p
    # add [T1p] -> [v]
    for T_def in positive_class_def:
        T = merge_all(T_def)
        G_dual.add_edge(dual(T), dual('v'))
    G_dual = nx.algorithms.transitive_closure(G_dual, False)

    # # # Trace and Trace Constraints
    # The trace maps an element x of M to a set of atoms in M*
    # Tr(x) = intersection_of_list([gla(dual(atom)) for atom in gla(x)])
    # therefore Tr(phi) = gla(dual(phi))
    # linearity of trace: Tr(merge(a,b)) = intersection(Tr(a), Tr(b))
    # (a < b) \implies (Tr(b) \subset Tr(a))
    # strategy:
    #   for v < T_i+ enforce Tr(T_i+) \subset Tr(v)
    #   for not (v < T_i-) enforce (Tr(T_i+) \not\subset Tr(v))

    # optional step: enforce \not ([T_i-] < [v]) by adding atoms \epsilon_i -> [T_i-]
    i = 1
    for T_def in negative_class_def:
        T = merge_all(T_def)
        G_dual.add_edge("epsilon_" + str(i), dual(T))
        i += 1
    del i
    G_dual = nx.algorithms.transitive_closure(G_dual, False)

    consistent = True
    for T_def in training_examples_def:
        for S_def in training_examples_def:
            if T_def is S_def:
                continue
            bond = True
            for c in T_def:
                if c not in S_def:
                    bond = False
                    break
            if bond:  # all c_i of T are also of S
                dualT, dualS = dual(merge_all(T_def)), dual(merge_all(S_def))
                neighs = G_dual(dualT)
                if dualS not in neighs:
                    consistent = False
                    break
        if not consistent:
            break

    # we can now reduce mutual elements to be the same: (a->b) and (b->a) \implies (a === b) todo

    # # # # # DEBUG
    print(pformat(list(G_dual.edges)))
    print(consistent)
    print(dual('v') in G_dual[dual(merge_all(T1n))])
    # nx.draw(G_master, pos=nx.spring_layout(G_master))
    nx.draw(G_dual, pos=nx.spring_layout(G_dual))
    matplotlib.pyplot.show()
    # print(pformat(len(G_master.edges)))
    # print(pformat(positive_class + negative_class))
    # print(less(G, atoms, merge_all(negative_class_def[2]), merge_all(positive_class_def[1])))
