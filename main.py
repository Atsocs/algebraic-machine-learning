import random
from pprint import pformat

from ordered_set import OrderedSet
import networkx as nx
import data
import numpy
from matplotlib import pyplot

# -------------------------------------
# NOTATION
# -------------------------------------


def d(x):
    # returns the dual of an node
    return '[' + x + ']'


def und(x):
    return x[1:-1]


def ggl(g, x):
    in_edges = g.in_edges(x)
    if not in_edges:
        return OrderedSet({x})
    in_neighs = list(zip(*in_edges))[0]
    return OrderedSet(in_neighs) | {x}


def master_gl(x):
    return ggl(master, x)


def dual_gl(x):
    return ggl(dual, x)


def master_atoms():
    types = nx.get_node_attributes(master, 'type')
    return OrderedSet(node for node in types if types[node] == 'atom')


def dual_atoms():
    types = nx.get_node_attributes(dual, 'type')
    return OrderedSet(node for node in types if types[node] == d('atom'))


def master_constants():
    return OrderedSet(data.get_constants()) | {data.target}


def dual_constants():
    terms = data.get_terms()
    terms = OrderedSet(terms[0] + terms[1])
    return OrderedSet(map(d, master_constants() | terms))


def master_gla(x):
    return master_gl(x) & master_atoms()


def dual_gla(x):
    return dual_gl(x) & dual_atoms()


def master_glc(x):
    return master_gl(x) & master_constants()


def dual_glc(x):
    return dual_gl(x) & dual_constants()


def ggu(g, x):
    return OrderedSet(g[x])


def master_gu(x):
    return ggu(master, x)


def dual_gu(x):
    return ggu(dual, x)


def master_u(x):
    return OrderedSet(y for y in master.nodes if less(x, y))


def dual_u(x):
    return OrderedSet(y for y in dual.nodes if less(x, y))


def trace(x):
    # x is always a member of master
    phi_set = master_gla(x)
    if phi_set:
        set_list = [dual_gla(d(phi)) for phi in phi_set]
        return OrderedSet.intersection(*set_list)
    else:
        return None
        # raise Exception("trace(x): x contains no atoms")


def master_dis(a, b):
    return master_gla(a) - master_gla(b)


def dual_dis(a, b):
    return dual_gla(a) - dual_gla(b)


def less(a, b):
    for phi in master_atoms():
        neighbors = master[phi]
        if not (
                (a not in neighbors) or (b in neighbors)
        ):
            return False
    return True


def same(ordered_set_1, ordered_set_2):
    if len(ordered_set_1) != len(ordered_set_2):
        return False
    for a in ordered_set_1:
        if a not in ordered_set_2:
            return False
    return True


def new_atom(name):
    return name


def new_phi():
    master.graph["phi"] += 1
    return r"$\phi_{" + str(master.graph["phi"]) + r"}$"


def new_Tphi():
    master.graph["Tphi"] += 1
    return r"$T\phi_{" + str(master.graph["Tphi"]) + r"}$"


def new_epsilon():
    master.graph["epsilon"] += 1
    return r"$\epsilon_{" + str(master.graph["epsilon"]) + r"}$"


def new_zeta():
    dual.graph["zeta"] += 1
    return r"$\zeta_{" + str(dual.graph["zeta"]) + r"}$"

# -------------------------------------
# CREATION
# -------------------------------------


def populate_master():
    global master
    consts = data.get_constants()
    terms = data.get_terms(expand=True)
    terms = terms[0] + terms[1]

    master.add_node(data.zero, type="atom")
    master.add_node(data.target, type="constant")
    master.add_edge(data.zero, data.target)
    master.add_nodes_from([(c, {"type": "constant"}) for c in consts])
    master.add_edges_from([(data.zero, c) for c in consts])
    master.add_nodes_from([(t[0], {"type": "term"}) for t in terms])
    for t in terms:
        master.add_edges_from([(c, t[0]) for c in t[1]])
    # add between-terms edges:
    for t in terms:
        for s in terms:
            if (t is not s) and (set(t[1])).issubset(set(s[1])):
                master.add_edge(t[0], s[0])

    master = nx.transitive_closure(master)
    master.graph["phi"] = 0
    master.graph["epsilon"] = 0
    master.graph["Tphi"] = 0


def populate_dual():
    global master, dual
    pos, neg = data.get_terms()
    terms = pos + neg

    dual = master.reverse()
    dual = nx.relabel_nodes(dual, dict(zip(dual.nodes, map(d, dual.nodes))))
    dual.add_node(data.zero_star, type=d('atom'))  # todo: fix all types
    dual.add_edges_from([(data.zero_star, d(t)) for t in terms])
    dual.add_edges_from([(d(p), d(data.target)) for p in pos])

    dual = nx.transitive_closure(dual)
    dual.graph["zeta"] = 0


def include_zeta_atoms():
    global dual
    pos, neg = data.get_terms()
    new_atoms = []
    for i in range(len(neg)):
        new_atoms.append(new_zeta())
    dual.add_nodes_from([(zeta, {"type": d("atom")}) for zeta in new_atoms])
    dual.add_edges_from(zip(new_atoms, map(d, neg)))
    dual = nx.transitive_closure(dual)
    return len(new_atoms)


# -------------------------------------
# ALGORITHMS
# -------------------------------------


def find_strongly_discriminant_constant(a, b):
    omega = OrderedSet(d(c) for c in master_glc(a))
    u = trace(b)
    while u:
        zeta = random.choice(tuple(u))
        u.remove(zeta)
        choose_from = omega - dual_gu(zeta)
        if choose_from:
            dc = random.choice(tuple(choose_from))
            return und(dc)
    return None


def algorithm1(r_neg):
    global master, dual
    zetas_added, phis_added = 0, 0
    for (a, b) in r_neg:
        if trace(b).issubset(trace(a)):
            while True:
                c = find_strongly_discriminant_constant(a, b)
                if c is None:
                    choose_from = dual_glc(d(b)) - dual_gl(d(a))
                    h = random.choice(tuple(choose_from))  # todo: maybe not random here is okay
                    zeta = new_zeta()
                    dual.add_node(zeta, type=d("atom"))
                    dual.add_edge(zeta, h)
                    dual = nx.transitive_closure(dual)
                    zetas_added += 1
                else:
                    break
            phi = new_atom(r"$\phi$")  # alternatively: phi = new_phi()
            master.add_node(phi, type="atom")
            dual.add_node(d(phi), type=d("dual-of-atom"))
            master.add_edge(phi, c)
            dual.add_edge(d(c), d(phi))
            master, dual = nx.transitive_closure(master), nx.transitive_closure(dual)
            phis_added += 1
    return phis_added, zetas_added


def algorithm2(r_pos):
    global master, dual
    phis_added = 0
    for (D, E) in r_pos:
        while not trace(E).issubset(trace(D)):
            choose_from = trace(E) - trace(D)
            zeta = random.choice(tuple(choose_from))
            gamma = OrderedSet(c for c in master_glc(E) if zeta not in dual_gl(d(c)))
            if not gamma:
                dual.add_edge(zeta, d(D))
                dual = nx.transitive_closure(dual)
            else:
                c = random.choice(tuple(gamma))
                epsilon = new_epsilon()
                master.add_node(epsilon, type="atom")
                dual.add_node(d(epsilon), type=d("dual-of-atom"))
                master.add_edge(epsilon, c)
                dual.add_edge(d(c), d(epsilon))
                master, dual = nx.transitive_closure(master), nx.transitive_closure(dual)
                phis_added += 1
    return phis_added


def delete_atoms(atoms_ordered_set):
    master.remove_nodes_from(atoms_ordered_set)
    dual.remove_nodes_from(map(d, atoms_ordered_set))
    return len(atoms_ordered_set)


def algorithm3(a, b):
    global master, dual
    psi_added, eps_prime_added = 0, 0
    A = master_gla(a) - master_gl(b)
    U = OrderedSet()
    for phi in A:
        U, B, delta = OrderedSet(), master_gla(b), dual_atoms() - dual_gl(d(phi))
        while True:
            eps = random.choice(tuple(B))
            delta_prime = delta & dual_gl(d(eps))
            if (not delta) or (not same(delta, delta_prime)):
                psi = new_phi()
                master.add_node(psi, type="atom")
                dual.add_node(d(psi), type=d("dual-of-atom"))
                edges = [(psi, phi), (psi, eps)]
                master.add_edges_from(edges)
                dual.add_edges_from([(d(y), d(x)) for (x, y) in edges])
                master, dual = nx.transitive_closure(master), nx.transitive_closure(dual)
                psi_added += 1

                delta = delta_prime
                U.append(eps)
            B.remove(eps)
            if not delta:
                break

    for eps in U:
        eps_prime = new_phi()
        master.add_node(eps_prime, type="atom")
        dual.add_node(d(eps_prime), type=d("dual-of-atom"))
        master.add_edge(eps_prime, eps)
        dual.add_edge(d(eps), d(eps_prime))
        master, dual = nx.transitive_closure(master), nx.transitive_closure(dual)
        eps_prime_added += 1

    removed = delete_atoms(U | A)
    return psi_added, eps_prime_added, removed


def algorithm4():
    # this is a stochastic algorithm to reduce the number of atoms.
    # this method can be called anytime, and should reduce the number of atoms in a few calls
    Q, A = OrderedSet(), [c for c in master_constants() if trace(c) is not None]
    while True:
        c = random.choice(tuple(A))
        A.remove(c)
        S = master_gl(c) & Q
        if not S:
            W = dual_atoms()
        else:
            W = OrderedSet.intersection(*[dual_gla(d(phi)) for phi in S])
        PHI = OrderedSet(d(phi) for phi in master_gla(c))
        tr = trace(c)
        while (tr is not None) and (not same(W, tr)):
            choose_from = W - tr
            xi = random.choice(tuple(choose_from))
            choose_from = PHI - dual_gu(xi)
            phi = und(random.choice(tuple(choose_from)))
            Q.append(phi)
            W &= dual_gla(d(phi))
        if not A:
            break
    removed = delete_atoms(master_atoms() - Q)

    return removed


def algorithm5(r_neg):
    Q, S = OrderedSet(), r_neg
    while S:
        r = random.choice(tuple(S))
        S.remove(r)
        a, b = r
        dis = dual_dis(d(b), d(a))
        if not (dis & Q):
            xi = random.choice(tuple(dis))
            Q.append(xi)
    removed = delete_atoms(dual_constants() - Q)

    return removed


def algorithm6(r_pin=None):
    if r_pin is None:
        r_pin = []
    for phi in master.nodes:
        H = master_constants() - master_u(phi)
        if not H:  # todo: review this
            continue
        T_phi = new_Tphi()
        composition = OrderedSet.union(*[master_gla(x) for x in H])
        master.add_node(T_phi, type="term", composition=composition)
        dual.add_node(d(T_phi), type="constant", composition=d(composition))
        for c in master_constants() & master_u(phi):
            r_pin.append((c, T_phi))
        return r_pin

# -------------------------------------
# DRAWING
# -------------------------------------


def draw_graph(g, node_size=900):
    types = list((nx.get_node_attributes(g, "type")).values())
    mapping = {"atom": "b", "constant": "g", "term": "r",
               d("atom"): "c", d("constant"): "m", d("dual-of-atom"): "y"}
    colors = [mapping[x] for x in types]
    # pos = nx.spring_layout(g)
    pos = nx.drawing.nx_agraph.graphviz_layout(g, prog="dot")
    g = nx.transitive_reduction(g)
    nx.draw(g, pos, node_size=node_size, edgecolors=colors, node_color="white", linewidths=2.0)
    nx.draw_networkx_labels(g, pos, font_size=10, font_family='sans-serif')


def draw(oneFigure=False):
    if oneFigure:
        fig = pyplot.figure(figsize=(8, 10))
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
        ax1.title.set_text('Master')
        ax2.title.set_text('Dual')
        pyplot.subplot(ax1)
        draw_graph(master, node_size=600)
        pyplot.subplot(ax2)
        draw_graph(dual, node_size=600)
        pyplot.show()
    else:
        pyplot.suptitle("Master")
        draw_graph(master)
        pyplot.show()
        pyplot.suptitle("Dual")
        draw_graph(dual)
        pyplot.show()


def main():
    enforce_negative_constraints = algorithm1
    enforce_positive_constraints = algorithm2
    sparse_crossing = algorithm3
    atom_set_reduction = algorithm4
    atom_set_reduction_for_the_dual_algebra = algorithm5
    generation_of_pinning_terms_and_relations = algorithm6

    pos, neg = data.get_terms()
    r_neg = [(data.target, n) for n in neg]
    r_pos = [(data.target, p) for p in pos]

    populate_master()
    populate_dual()
    # draw()
    include_zeta_atoms()
    # draw()
    r_pin = generation_of_pinning_terms_and_relations()  # fixme: this produces no effect somehow
    print(r_pin)
    # draw()
    enforce_negative_constraints(r_neg)
    # draw()
    enforce_positive_constraints(r_pos)
    # draw()
    sparse_crossing(data.target, pos[0])
    sparse_crossing(data.target, pos[1])
    # draw()
    for i in range(3):
        atom_set_reduction()
        # draw()
    for i in range(3):
        atom_set_reduction_for_the_dual_algebra(r_neg)
        # draw()
    draw(True)

    print("master.edges:", pformat(list(master.edges)))
    print("dual.edges:", pformat(list(dual.edges)))
    print("Pos Tr.Constr.:", pformat([trace(p).issubset(trace(data.target)) for p in pos]))
    print("Neg Tr.Constr.:", pformat([not trace(n).issubset(trace(data.target)) for n in neg]))
    print("R+:", pformat([less(data.target, p) for p in pos]))
    print("R-:", pformat([not less(data.target, n) for n in neg]))

    return len(master_atoms()), len(dual_atoms())


if __name__ == '__main__':
    seed = 123456
    master = nx.DiGraph()
    dual = nx.DiGraph()

    numpy.random.seed(seed)
    random.seed(seed)

    print(main())
