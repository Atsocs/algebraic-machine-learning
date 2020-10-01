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


def draw_colored_graph(g, show=False, spring=False):
    if spring:
        nx.draw(g, pos=nx.spring_layout(g), node_color=list((nx.get_node_attributes(g, "color")).values()))
    else:
        nx.draw(g, pos=nx.random_layout(g), node_color=list((nx.get_node_attributes(g, "color")).values()))

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
    m = len(dual_atoms)
    new_atoms = ["zeta_" + str(i + m) for i in range(len(neg_def))]
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
        raise Exception('Tr is callable only with "master" or "dual" as S argument.')

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
        raise Exception('GLa is callable only with "master" or "dual" as S argument.')
    return GL(g, x) & A(S)


def A(S):
    if S == "master":
        return OrderedSet(master_atoms)
    if S == "dual":
        return OrderedSet(dual_atoms)
    raise Exception('A is callable only with "master" or "dual" as S argument.')


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
    raise Exception('C is callable only with "master" or "dual" as S argument.')


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
    nodes_created = 0
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
                    nodes_created += 1
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
            nodes_created += 1
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
    return nodes_created


def enforce_positive_trace_constraints(just_one_pos_term=False, t=None):
    global master, dual, master_atoms, dual_atoms
    nodes_created = 0
    d = target
    if not just_one_pos_term:
        pos_terms_to_enforce = map(term, pos_def)
    else:
        pos_terms_to_enforce = [t]
    for e in pos_terms_to_enforce:
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
                nodes_created += 1
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

            tr_e = Tr("master", e)
            tr_d = Tr("master", d)
            is_subSet = tr_e.issubset(tr_d)
            Set_diff = tr_e - tr_d
    return nodes_created


def same(ordered_set_1, ordered_set_2):
    if len(ordered_set_2) != len(ordered_set_1):
        return False
    for a in ordered_set_1:
        if a not in ordered_set_2:
            return False
    return True


def sparse_crossing(a, b):
    global master, dual, master_atoms
    nodes_created = 0
    AA = GLa("master", a) - GL(master, b)
    U = OrderedSet()
    for phi in AA:
        U, B, Del = OrderedSet(), GLa("master", b), A("dual") - GL(dual, dual_of(phi))
        while True:
            eps = random.choice(tuple(B))
            Del_prime = Del & GL(dual, dual_of(eps))
            if (not Del) or (not same(Del_prime, Del)):
                psi = "psi_" + str(len(master_atoms) + 1)
                nodes_created += 1
                master.add_edge(psi, phi)
                master.add_edge(psi, eps)
                dual.add_edge(dual_of(phi), dual_of(psi))
                dual.add_edge(dual_of(eps), dual_of(psi))
                master_atoms.append(psi)
                if closure:
                    master = nx.algorithms.transitive_closure(master)
                    dual = nx.algorithms.transitive_closure(dual)
                # ---- COLORING ------
                nx.set_node_attributes(master, {psi: 'orchid'}, "color")
                nx.set_node_attributes(dual, {dual_of(psi): 'orchid'}, "color")
                # --------------------

                Del = Del_prime
                U.append(eps)
            B.remove(eps)
            if not Del:
                break
    for eps in U:
        eps_prime = "eps_prime_" + str(len(master_atoms) + 1)
        nodes_created += 1
        master.add_edge(eps_prime, eps)
        dual.add_edge(dual_of(eps_prime), dual_of(eps))
        master_atoms.append(eps_prime)
        if closure:
            master = nx.algorithms.transitive_closure(master)
            dual = nx.algorithms.transitive_closure(dual)
        # ---- COLORING ------
        nx.set_node_attributes(master, {eps_prime: 'orchid'}, "color")
        nx.set_node_attributes(dual, {dual_of(eps_prime): 'orchid'}, "color")
        # --------------------
    nodes_to_remove = U | AA
    nodes_removed = len(nodes_to_remove)
    duals = map(dual_of, nodes_to_remove)

    master.remove_nodes_from(U | AA)
    dual.remove_nodes_from(duals)
    master_atoms = [elem for elem in master_atoms if (elem not in nodes_to_remove)]

    if closure:
        master = nx.algorithms.transitive_closure(master)
        dual = nx.algorithms.transitive_closure(dual)

    ret = nodes_created - nodes_removed
    ret += enforce_negative_trace_constraints()
    ret += enforce_positive_trace_constraints(True, b)
    return ret


def sparse_crossing_all_positive_relations():
    a = target
    ret = 0
    for t_def in pos_def:
        b = term(t_def)
        ret += sparse_crossing(a, b)
    return ret


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


def less(S, a, b):
    global master, dual, master_atoms, dual_atoms
    if S == "master":
        graph, atoms_list = master, master_atoms
    elif S == "dual":
        graph, atoms_list = dual, dual_atoms
    else:
        raise Exception('less is callable only with "master" or "dual" as S argument.')

    # a < b iff (for all phi in atoms_list, (phi not -> a) or (phi -> b))
    for phi in atoms_list:
        neighbors = graph[phi]
        if not (
                (a not in neighbors) or (b in neighbors)
        ):
            return False
    return True


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

    while True:
        n = 0
        n += enforce_negative_trace_constraints()
        n += enforce_positive_trace_constraints()
        if n > 0:
            print("enforcement creates", n, "atoms")
        else:
            break

    n = sparse_crossing_all_positive_relations()
    print("sparse_crossing creates minus removes", n, "atoms")

    drawing = True
    if drawing:
        fig = pyplot.figure(figsize=(10, 5))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        ax1.title.set_text('Master')
        ax2.title.set_text('Dual')
        pyplot.subplot(ax1)
        draw_colored_graph(master)
        pyplot.subplot(ax2)
        draw_colored_graph(dual)
        pyplot.show()

    print("#master_atoms =", len(master_atoms))
    print("#dual_atoms =", len(dual_atoms))
    print("#master_edges =", len(master.edges))
    print("#dual_edges =", len(dual.edges))
    print("(target < T_i_neg) =", [less("master", target, term(x)) for x in neg_def])  # Should be a list of False's
    print("(target < T_i_pos) =", [less("master", target, term(x)) for x in pos_def])

    print_edges = False
    if print_edges:
        print(pformat(list(master.edges)))
