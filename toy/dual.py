import networkx as nx

from toy.master import Master
from toy import data
from toy.LaTeX import node_code


def d(x):
    return '[' + str(x) + ']'


def dtype(x):
    if x == "atom":
        return "dual-of-atom"
    elif x == "term" or x == "constant":
        return "constant"
    else:
        assert False


class Dual:
    """Dual Algebra and Graph"""

    def add_node(self, kind, node, latex):
        assert (node not in self.graph.nodes)
        kinds = {"atom": self.atoms, "constant": self.constants, "dual-of-atom": self.dual_of_atoms}
        assert (kind in kinds)
        node_list = kinds[kind]

        if latex is None:
            latex = node
        self.graph.add_node(node, type=kind, latex=latex)
        node_list.append(node)
        return node

    add_atom = Master.__dict__['add_atom']
    add_edge = Master.__dict__['add_edge']
    close_graph = Master.__dict__['close_graph']

    def draw(self, avoid_clutter=True, node_size=900, pos=None):
        types = list((nx.get_node_attributes(self.graph, "type")).values())
        mapping = {"atom": "c", "constant": "m", "dual-of-atom": "y"}
        colors = [mapping[x] for x in types]
        if pos is None:
            pos = self.get_pos()
        labels_items = list((nx.get_node_attributes(self.graph, "latex")).values())
        assert len(pos) == len(labels_items)
        labels = dict(zip(pos.keys(), labels_items))
        if avoid_clutter:
            g = nx.transitive_reduction(self.graph)
        else:
            g = self.graph
        nx.draw(g, pos, node_size=node_size, edgecolors=colors, node_color="white", linewidths=1.0)
        nx.draw_networkx_labels(g, pos, labels=labels, font_size=10, font_family='sans-serif')

    def get_pos(self):
        pos = nx.drawing.nx_agraph.graphviz_layout(self.graph, prog="dot")
        y_values = sorted(list(set(y for (x, y) in pos.values())))
        new_pos = []
        for level in range(len(y_values)):
            pairs = [(n, (x, y)) for (n, (x, y)) in pos.items() if y == y_values[level]]
            n, p = zip(*pairs)
            n, p = sorted(n), sorted(p)
            new_pos += list(zip(n, p))
        pos = dict(new_pos)
        flipped_pos = {node: (x, -y) for (node, (x, y)) in pos.items()}
        return flipped_pos

    def reverted_negative_relations(self):
        for (i, ne) in enumerate(data.negative_examples):
            zeta = self.add_atom(f"zeta_{i}", r"$\zeta_{" + f'{i}' + "}$")
            self.add_edge(zeta, d(ne))
        self.close_graph()

    def __init__(self, master):
        self.constants = []  # duals of master's constants or dual_of_atoms
        self.dual_of_atoms = []  # duals of master's atoms
        self.atoms = []  # not duals of any element in the master

        self.constants += [d(c) for c in master.constants]
        self.constants += [d(t) for t in master.terms]
        self.dual_of_atoms += [d(a) for a in master.atoms]

        g = master.graph.reverse()
        g = nx.relabel_nodes(g, dict(zip(g.nodes, map(d, g.nodes))))
        for n in g.nodes:
            g.nodes[n]['type'] = dtype(g.nodes[n]['type'])
            g.nodes[n]['latex'] = '$' + d('{' + g.nodes[n]['latex'][1: -1] + '}') + '$'
        self.graph = g

        self.add_atom(data.zero_star, node_code.dual_atom(data.zero_star))
        for t in master.terms:
            self.add_edge(data.zero_star, d(t))
        for pe in data.positive_examples:
            self.add_edge(d(pe), d(data.target))

        self.close_graph()
