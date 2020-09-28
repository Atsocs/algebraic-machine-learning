import random
import numpy
import networkx as nx
from pprint import pformat
from matplotlib import pyplot
from itertools import accumulate
from enum import Enum
from ordered_set import OrderedSet

seed = 123456  # or any fixed integer
random.seed(seed)
numpy.random.seed(seed)

target, zero, zero_star = 'v', '0', '0*'
closure = True


def draw_colored_graph(g, show=False):
    nx.draw(g, node_color=list((nx.get_node_attributes(g, "color")).values()))
    if show:
        pyplot.show()


class Color(Enum):
    BLACK = 0
    WHITE = 1
    BOTH = 2

    def __add__(self, other):
        if self == other:
            return self
        if self == Color.BOTH or other == Color.BOTH:
            return Color.BOTH
        if self == Color.WHITE and other == Color.BLACK:
            return Color.BOTH
        if self == Color.BLACK and other == Color.WHITE:
            return Color.BOTH

    def __lt__(self, other):
        if other is None:
            print("wow")
            return False

    def __gt__(self, other):
        if other is None:
            print("wow2")
            return True


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


def term(term_def, return_only_name=True):
    name, definition = term_def
    if return_only_name:
        return name
    else:
        return name, tuple(accumulate(definition, merge))[-1]


def dual_of(node):
    return "dual", node


def un_dual_of(dual_of_node):
    return dual_of_node[1]


def subSet_def(t_def, s_def):
    # is T subSet of S?
    is_subSet = True
    for c in t_def[1]:
        if c not in s_def[1]:
            is_subSet = False
            break
    return is_subSet


def create_master():
    training_def = pos_def + neg_def

    master_atoms_list = [zero]
    master_graph = nx.DiGraph()
    master_graph.add_edge(zero, target)
    master_graph.add_edges_from([(zero, c) for c in consts])
    for t_def in training_def:
        name, t = term(t_def, False)
        master_graph.add_edges_from([(c, name) for c in t_def[1]])
    for t_def in training_def:
        for s_def in training_def:
            if s_def is t_def:
                continue
            if subSet_def(t_def, s_def):
                master_graph.add_edge(term(t_def), term(s_def))

    if closure:
        master_graph = nx.algorithms.transitive_closure(master_graph)

    # -------- COLORING --------
    nx.set_node_attributes(master_graph, "k", "color")
    nx.set_node_attributes(master_graph, {target: 'gray', zero: 'b'}, "color")

    dict_pos = {term(t_pos): 'lime' for t_pos in pos_def}
    nx.set_node_attributes(master_graph, dict_pos, "color")

    dict_neg = {term(t_neg): 'red' for t_neg in neg_def}
    nx.set_node_attributes(master_graph, dict_neg, "color")
    # --------------------------

    return master_graph, master_atoms_list


def create_dual():
    training_def = pos_def + neg_def
    dual_atoms_list = [zero_star]

    dual_graph = master.reverse()
    dual_graph = nx.relabel_nodes(dual_graph, dict(zip(dual_graph.nodes, map(dual_of, dual_graph.nodes))))
    dual_graph.add_edges_from([(zero_star, dual_of(term(t_def))) for t_def in training_def])
    dual_graph.add_edges_from([(dual_of(term(t_def)), dual_of(target)) for t_def in pos_def])

    if closure:
        dual_graph = nx.algorithms.transitive_closure(dual_graph)

    # -------- COLORING --------
    nx.set_node_attributes(dual_graph, "k", "color")
    nx.set_node_attributes(dual_graph, {dual_of(target): 'gray', zero_star: 'aqua', dual_of(zero): 'b'}, "color")

    dict_pos = {dual_of(term(t_pos)): 'lime' for t_pos in pos_def}
    nx.set_node_attributes(dual_graph, dict_pos, "color")

    dict_neg = {dual_of(term(t_neg)): 'red' for t_neg in neg_def}
    nx.set_node_attributes(dual_graph, dict_neg, "color")
    # --------------------------

    return dual_graph, dual_atoms_list


def include_zeta_atoms():
    global dual, dual_atoms
    n = len(dual_atoms)
    new_atoms = ["zeta_" + str(i + n) for i in range(len(neg_def))]
    duals_of_t_neg = [dual_of(term(t_def)) for t_def in neg_def]
    dual.add_edges_from((zip(new_atoms, duals_of_t_neg)))

    if closure:
        dual = nx.algorithms.transitive_closure(dual)

    # -------- COLORING --------
    nx.set_node_attributes(dual, {atom: 'orange' for atom in new_atoms}, "color")
    # --------------------------

    dual_atoms += new_atoms


def there_is_edge(g, node1, node2):
    neighs = list(g[node1])
    return node2 in neighs


def Tr(S, x):
    if S == "master":
        g, S2 = master, "dual"
    elif S == "dual":
        g, S2 = dual, "master"
    else:
        raise Exception('Tr is callable only with "master" or "dual" as argument.')

    phi_Set = GLa(S, x)
    Setlist = []
    for phi in phi_Set:
        gla = GLa(S2, dual_of(phi))
        Setlist.append(gla)

    return OrderedSet.intersection(*Setlist)


def GU(g, node):
    return OrderedSet(g[node])


def GL(g, x):
    in_neighs = list(zip(*g.in_edges(x)))[0]
    return OrderedSet(in_neighs) | {x}


def GLa(S, x):
    if S == "master":
        g = master
    elif S == "dual":
        g = dual
    else:
        raise Exception('GLa is callable only with "master" or "dual" as argument.')
    return GL(g, x) & A(S)


def A(S):
    if S == "master":
        return OrderedSet(master_atoms)
    if S == "dual":
        return OrderedSet(dual_atoms)
    raise Exception('A is callable only with "master" or "dual" as argument.')


def C(S):
    if S == "master":
        ret = consts.copy()
        ret += [target]
        return OrderedSet(ret)
    if S == "dual":
        ret = [dual_of(x) for x in consts]
        ret.append(dual_of(target))
        ret += [dual_of(term(t_def)) for t_def in pos_def + neg_def]
        return OrderedSet(ret)
    raise Exception('C is callable only with "master" or "dual" as argument.')


def random_pop(Set_s):
    x = random.choice(tuple(Set_s))
    Set_s.remove(x)
    return x, Set_s


def find_strongly_discriminant_constant(a, b):
    global master_atoms, dual_atoms
    omega_a = OrderedSet(dual_of(c) for c in GL(master, a) & C("master"))
    u = Tr("master", b)
    while u:
        zeta, u = random_pop(u)
        if omega_a - GU(dual, zeta):
            dual_of_c = random.choice(tuple(omega_a - GU(dual, zeta)))
            return un_dual_of(dual_of_c)
    return None


def enforce_negative_trace_constraints():
    global master, dual, master_atoms, dual_atoms
    a = target
    for t_def in neg_def:
        b = term(t_def)
        tr_b = Tr("master", b)
        tr_a = Tr("master", a)
        if tr_b.issubset(tr_a):
            while True:
                c = find_strongly_discriminant_constant(a, b)
                if c is None:
                    s1 = GL(dual, dual_of(b))
                    s2 = C("dual")
                    GLc_dual_of_b = s1 & s2
                    choose_h_from_Set = GLc_dual_of_b - GL(dual, dual_of(a))
                    h = random.choice(tuple(choose_h_from_Set))

                    zeta = "zeta_" + str(len(dual_atoms) + 1)
                    dual_atoms.append(zeta)
                    dual.add_edge(zeta, h)
                    if closure:
                        dual = nx.algorithms.transitive_closure(dual)
                    # ---- COLORING ------
                    nx.set_node_attributes(dual, {zeta: 'sienna'}, "color")
                    # --------------------
                else:
                    break
            phi = "phi_" + str(len(master_atoms))
            master_atoms.append(phi)
            master.add_edge(phi, c)
            dual.add_edge(dual_of(c), dual_of(phi))
            if closure:
                master = nx.algorithms.transitive_closure(master)
                dual = nx.algorithms.transitive_closure(dual)
            # ---- COLORING ------
            nx.set_node_attributes(master, {phi: 'orchid'}, "color")
            nx.set_node_attributes(dual, {dual_of(phi): 'orchid'}, "color")
            # --------------------


def enforce_positive_trace_constraints():
    global master, dual, master_atoms, dual_atoms
    d = target
    for t_def in pos_def:
        e = term(t_def)
        tr_e = Tr("master", e)
        tr_d = Tr("master", d)
        is_subSet = tr_e.issubset(tr_d)
        Set_diff = tr_e - tr_d
        while not is_subSet:
            zeta = random.choice(tuple(Set_diff))
            c_is_in = GL(master, e) & C("master")
            Gamma = OrderedSet()
            for c in c_is_in:
                gl = GL(dual, dual_of(c))
                if zeta not in gl:
                    Gamma.add(c)
            if not Gamma:
                dual.add_edge(zeta, dual_of(d))
                if closure:
                    dual = nx.algorithms.transitive_closure(dual)
            else:
                c = random.choice(tuple(Gamma))

                phi = "phi_" + str(len(master_atoms))
                master_atoms.append(phi)
                master.add_edge(phi, c)
                dual.add_edge(dual_of(c), dual_of(phi))
                if closure:
                    master = nx.algorithms.transitive_closure(master)
                    dual = nx.algorithms.transitive_closure(dual)
                # ---- COLORING ------
                nx.set_node_attributes(master, {phi: 'purple'}, "color")
                nx.set_node_attributes(dual, {dual_of(phi): 'purple'}, "color")
                # --------------------

            tr_e = Tr("master", e)
            tr_d = Tr("master", d)
            is_subSet = tr_e.issubset(tr_d)
            Set_diff = tr_e - tr_d


def is_consistent():
    global dual
    # function is not working properly. If we test against a non-consistent dataSet we get True anyway
    consistent = True
    training_def = pos_def + neg_def

    for t1def in training_def:
        for t2def in training_def:
            if t2def is t1def:
                continue
            if subSet_def(t1def, t2def):
                dual_of_t1, dual_of_t2 = dual_of(term(t1def)), dual_of(term(t2def))
                if not there_is_edge(dual, dual_of_t2, dual_of_t1):
                    consistent = False
                    break
        if not consistent:
            break
    return consistent


if __name__ == '__main__':
    # ------- DATA ---------
    consts = [
        ((None, None), (Color.BLACK, None)),
        ((Color.BLACK, None), (None, None)),
        ((None, Color.BLACK), (None, None)),
        ((None, None), (None, Color.BLACK)),
        ((None, None), (Color.WHITE, None)),
        ((Color.WHITE, None), (None, None)),
        ((None, Color.WHITE), (None, None)),
        ((None, None), (None, Color.WHITE)),
    ]

    T1_pos = ('T1_pos', {consts[2 - 1], consts[7 - 1], consts[1 - 1], consts[8 - 1]})
    T2_pos = ('T2_pos', {consts[6 - 1], consts[3 - 1], consts[5 - 1], consts[4 - 1]})
    T1_neg = ('T1_neg', {consts[2 - 1], consts[7 - 1], consts[5 - 1], consts[4 - 1]})
    T2_neg = ('T2_neg', {consts[6 - 1], consts[3 - 1], consts[5 - 1], consts[8 - 1]})
    T3_neg = ('T3_neg', {consts[6 - 1], consts[7 - 1], consts[5 - 1], consts[4 - 1]})

    pos_def, neg_def = [T1_pos, T2_pos], [T1_neg, T2_neg, T3_neg]
    # ----------------------

    master, master_atoms = create_master()

    dual, dual_atoms = create_dual()
    include_zeta_atoms()
    # todo: fix and print(is_consistent(dual, pos_class, neg_class))
    #  merge pairs of dual nodes which are neighbor of each other. Two elements in M can share the same dual

    enforce_negative_trace_constraints()
    enforce_positive_trace_constraints()
    draw_colored_graph(master, show=True)
    draw_colored_graph(dual, show=True)

    print("#master_atoms =", len(master_atoms))
    print("#dual_atoms =", len(dual_atoms))
    print(pformat(list(master.edges)))
