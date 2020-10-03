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
    in_neighs = list(zip(*g.in_edges(x)))[0]
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


def trace(x):
    # x is always member of master
    phi_set = master_gla(x)
    set_list = [dual_gla(d(phi)) for phi in phi_set]
    return OrderedSet.intersection(*set_list)


def less(a, b):
    for phi in master_atoms():
        neighbors = master[phi]
        if not (
                (a not in neighbors) or (b in neighbors)
        ):
            return False
    return True


def new_atom(name):
    return name


def new_phi():
    master.graph["phi"] += 1
    return r"$\phi_{" + str(master.graph["phi"]) + r"}$"


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


def populate_dual():
    global master, dual
    pos, neg = data.get_terms()
    terms = pos + neg

    dual = master.reverse()
    dual = nx.relabel_nodes(dual, dict(zip(dual.nodes, map(d, dual.nodes))))
    dual.add_node(data.zero_star, type=d('atom'))
    dual.add_edges_from([(data.zero_star, d(t)) for t in terms])
    dual.add_edges_from([(d(p), d(data.target)) for p in pos])

    dual = nx.transitive_closure(dual)
    dual.graph["zeta"] = 0


def include_zeta_atoms():
    global dual
    pos, neg = data.get_terms()
    new_atoms = [new_zeta() for i in range(len(neg))]
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


def draw_graph(g):
    types = list((nx.get_node_attributes(g, "type")).values())
    mapping = {"atom": "b", "constant": "g", "term": "r",
               d("atom"): "c", d("constant"): "m", d("dual-of-atom"): "y"}
    colors = [mapping[x] for x in types]
    pos = nx.circular_layout(g)
    g = nx.transitive_reduction(g)
    nx.draw(g, pos, node_size=1000, edgecolors=colors, node_color="white", linewidths=2.0)
    nx.draw_networkx_labels(g, pos, font_size=10, font_family='sans-serif')


def draw():
    fig = pyplot.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    ax1.title.set_text('Master')
    ax2.title.set_text('Dual')
    pyplot.subplot(ax1)
    draw_graph(master)
    pyplot.subplot(ax2)
    draw_graph(dual)
    pyplot.show()


def main():
    seed = 123456
    numpy.random.seed(seed)
    random.seed(seed)

    enforce_negative_constraints = algorithm1
    enforce_positive_constraints = algorithm2

    pos, neg = data.get_terms()

    populate_master()
    populate_dual()
    # draw()
    include_zeta_atoms()
    # draw()
    r_neg = [(data.target, n) for n in neg]
    enforce_negative_constraints(r_neg)
    # draw()
    r_pos = [(data.target, p) for p in pos]
    enforce_positive_constraints(r_pos)
    draw()
    print("master.edges:", pformat(list(master.edges)))
    print("dual.edges:", pformat(list(dual.edges)))
    print("R+ satisfied:", pformat([less(data.target, p) for p in pos]))
    print("R- satisfied:", pformat([not less(data.target, n) for n in neg]))


if __name__ == '__main__':
    master = nx.DiGraph()
    dual = nx.DiGraph()
    main()
