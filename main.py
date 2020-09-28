import random
from itertools import accumulate
from enum import Enum
from pprint import pformat
import networkx as nx
from matplotlib import pyplot


def draw_colored_graph(g, show=False):
    nx.draw(g, node_color=list((nx.get_node_attributes(g, "color")).values()))
    if show:
        pyplot.show()


target, zero, zero_star = 'v', '0', '0*'
closure = True


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


def subset_def(t_def, s_def):
    # is T subset of S?
    is_subset = True
    for c in t_def[1]:
        if c not in s_def[1]:
            is_subset = False
            break
    return is_subset


def create_master(consts, pos_def, neg_def):
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
            if subset_def(t_def, s_def):
                master_graph.add_edge(term(t_def), term(s_def))

    if closure:
        master_graph = nx.algorithms.transitive_closure(master_graph)

    # -------- COLORING --------
    nx.set_node_attributes(master_graph, "k", "color")
    nx.set_node_attributes(master_graph, {target: 'fuchsia', zero: 'b'}, "color")

    dict_pos = {term(t_pos): 'lime' for t_pos in pos_def}
    nx.set_node_attributes(master_graph, dict_pos, "color")

    dict_neg = {term(t_neg): 'red' for t_neg in neg_def}
    nx.set_node_attributes(master_graph, dict_neg, "color")
    # --------------------------

    return master_graph, master_atoms_list


def create_dual(master_graph, pos_def, neg_def):
    training_def = pos_def + neg_def
    dual_atoms_list = [zero_star]

    dual_graph = master_graph.reverse()
    dual_graph = nx.relabel_nodes(dual_graph, dict(zip(dual_graph.nodes, map(dual_of, dual_graph.nodes))))
    dual_graph.add_edges_from([(zero_star, dual_of(term(t_def))) for t_def in training_def])
    dual_graph.add_edges_from([(dual_of(term(t_def)), dual_of(target)) for t_def in pos_def])

    if closure:
        dual_graph = nx.algorithms.transitive_closure(dual_graph)

    # -------- COLORING --------
    nx.set_node_attributes(dual_graph, "k", "color")
    nx.set_node_attributes(dual_graph, {dual_of(target): 'fuchsia', zero_star: 'aqua', dual_of(zero): 'b'}, "color")

    dict_pos = {dual_of(term(t_pos)): 'lime' for t_pos in pos_def}
    nx.set_node_attributes(dual_graph, dict_pos, "color")

    dict_neg = {dual_of(term(t_neg)): 'red' for t_neg in neg_def}
    nx.set_node_attributes(dual_graph, dict_neg, "color")
    # --------------------------

    return dual_graph, dual_atoms_list


def include_zeta_atoms(dual_graph, dual_atoms_list, neg_def):
    n = len(dual_atoms_list)
    new_atoms = ["zeta_" + str(i + n) for i in range(len(neg_def))]
    duals_of_t_neg = [dual_of(term(t_def)) for t_def in neg_def]
    dual_graph.add_edges_from((zip(new_atoms, duals_of_t_neg)))

    if closure:
        dual_graph = nx.algorithms.transitive_closure(dual_graph)

    # -------- COLORING --------
    nx.set_node_attributes(dual_graph, {atom: 'silver' for atom in new_atoms}, "color")
    # --------------------------

    dual_atoms_list += new_atoms
    return dual_graph, dual_atoms_list


def there_is_edge(g, node1, node2):
    neighs = list(g[node1])
    return node2 in neighs


def Tr(S, x):
    if S == "master":
        g, g2, S2 = master, dual, "dual"
    elif S == "dual":
        g, g2, S2 = dual, master, "master"
    else:
        raise Exception('A is callable only with "master" or "dual" as argument.')

    phi_set = GLa(g, S, x)
    setlist = [GLa(g2, S2, dual_of(phi)) for phi in phi_set]

    return set.intersection(*setlist)


def GU(g, node):
    return set(g[node])


def GL(g, x):
    in_neighs = list(zip(*g.in_edges(x)))[0]
    return set(in_neighs)


def GLa(g, S, x):
    return GL(g, x) & A(S)


def A(S):
    if S == "master":
        return set(master_atoms)
    if S == "dual":
        return set(dual_atoms)
    raise Exception('A is callable only with "master" or "dual" as argument.')


def C(S):
    if S == "master":
        ret = constants.copy()
        ret += [target]
        return set(ret)
    if S == "dual":
        ret = [dual_of(x) for x in constants]
        ret.append(dual_of(target))
        ret += [dual_of(term(t_def)) for t_def in pos_class + neg_class]
        return set(ret)
    raise Exception('C is callable only with "master" or "dual" as argument.')


def random_pop(set_s):
    x = random.choice(tuple(set_s))
    set_s = set_s - {x}
    return x, set_s


def find_strongly_discriminant_constant(a, b):
    omega_a = {dual_of(c) for c in GL(master, a) & C("master")}
    u = Tr("master", b)
    while u:
        zeta, u = random_pop(u)
        if omega_a - GU(dual, zeta):
            dual_of_c = random.choice(tuple(omega_a - GU(dual, zeta)))
            return un_dual_of(dual_of_c)
    return None


def enforce_negative_trace_constraints(neg_def, master_atoms_list, dual_atoms_list):
    # todo: make this work!
    a = target
    for t_def in neg_def:
        b = term(t_def)
        if Tr("master", b).issubset(Tr("master", a)):
            while True:
                c = find_strongly_discriminant_constant(a, b)
                if c is None:
                    s1 = GL(dual, dual_of(b))
                    s2 = C("dual")
                    GLc_dual_of_b = s1 & s2
                    choose_h_from_set = GLc_dual_of_b - GL(dual, dual_of(a))
                    h = random.choice(tuple(choose_h_from_set))

                    zeta = "zeta_" + str(len(dual_atoms_list) + 1)
                    dual_atoms_list.append(zeta)
                    dual.create_edge(zeta, h)
                else:
                    break
            phi = "phi_" + str(len(master_atoms_list) + 1)
            master_atoms_list.append(phi)
            master.create_edge(phi, c)

    return master_atoms_list, dual_atoms_list


def is_consistent(dual_graph, pos_def, neg_def):
    # fixme: function not working. If we test against a non-consistent dataset we get True anyway
    consistent = True
    training_def = pos_def + neg_def

    for t1def in training_def:
        for t2def in training_def:
            if t2def is t1def:
                continue
            if subset_def(t1def, t2def):
                dual_of_t1, dual_of_t2 = dual_of(term(t1def)), dual_of(term(t2def))
                if not there_is_edge(dual_graph, dual_of_t2, dual_of_t1):
                    consistent = False
                    break
        if not consistent:
            break
    return consistent


if __name__ == '__main__':
    # ------- DATA ---------
    constants = [
        ((None, None), (Color.BLACK, None)),
        ((Color.BLACK, None), (None, None)),
        ((None, Color.BLACK), (None, None)),
        ((None, None), (None, Color.BLACK)),
        ((None, None), (Color.WHITE, None)),
        ((Color.WHITE, None), (None, None)),
        ((None, Color.WHITE), (None, None)),
        ((None, None), (None, Color.WHITE)),
    ]

    T1_pos = ('T1_pos', {constants[2 - 1], constants[7 - 1], constants[1 - 1], constants[8 - 1]})
    T2_pos = ('T2_pos', {constants[6 - 1], constants[3 - 1], constants[5 - 1], constants[4 - 1]})
    T1_neg = ('T1_neg', {constants[2 - 1], constants[7 - 1], constants[5 - 1], constants[4 - 1]})
    T2_neg = ('T2_neg', {constants[6 - 1], constants[3 - 1], constants[5 - 1], constants[8 - 1]})
    T3_neg = ('T3_neg', {constants[6 - 1], constants[7 - 1], constants[5 - 1], constants[4 - 1]})

    pos_class, neg_class = [T1_pos, T2_pos], [T1_neg, T2_neg, T3_neg]

    # ----------------------

    # c_m_star = [dual_of(target)] + [dual_of(x) for x in c_m] + [dual_of(term(x)) for x in pos_class + neg_class]

    master, master_atoms = create_master(constants, pos_class, neg_class)
    # draw_colored_graph(master, show=True)

    dual, dual_atoms = create_dual(master, pos_class, neg_class)
    dual, dual_atoms = include_zeta_atoms(dual, dual_atoms, neg_class)
    # todo: fix and print(is_consistent(dual, pos_class, neg_class))
    # todo: merge pairs of dual nodes which are neighbor of each other. Two elements in M can share the same dual
    # draw_colored_graph(dual, show=True)

    # todo: master_atoms, dual_atoms = enforce_negative_trace_constraints(neg_class, master_atoms, dual_atoms)

    # todo: enforce_positive_trace_constraints
