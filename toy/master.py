import os

import networkx as nx
from matplotlib import pyplot
from ordered_set import OrderedSet

from toy.LaTeX import node_code
from toy import data, drawing


class Master:
    """Master Algebra and Graph"""

    def add_node(self, kind, node, latex):
        assert (node not in self.graph.nodes)
        assert (kind in self.kinds)
        node_list = self.kinds[kind]

        if latex is None:
            latex = node
        self.graph.add_node(node, type=kind, latex=latex)
        node_list.append(node)
        if self.draw_flag:
            self.draw_and_save()
        return node

    def add_constant(self, constant, latex):
        return self.add_node("constant", constant, latex)

    def add_term(self, term, latex):
        return self.add_node("term", term, latex)

    def add_atom(self, atom, latex):
        return self.add_node("atom", atom, latex)

    def add_edge(self, a, b):
        assert (a in self.graph.nodes)
        assert (b in self.graph.nodes)
        ret = self.graph.add_edge(a, b)
        if self.draw_flag:
            self.draw_and_save()
        return ret

    def remove_atom(self, atom):
        self.atoms.remove(atom)
        self.graph.remove_node(atom)
        if self.draw_flag:
            self.draw_and_save()
        return True

    def remove_atoms_from(self, atom_list):
        removed = 0
        for atom in atom_list:
            removed += self.remove_atom(atom)
        return removed

    def close_graph(self):
        self.graph = nx.transitive_closure(self.graph)
        return self.graph

    def gl(self, x):
        in_edges = self.graph.in_edges(x)
        if not in_edges:
            return OrderedSet({x})
        in_neighs = list(zip(*in_edges))[0]
        return OrderedSet(in_neighs) | {x}

    def gla(self, x):
        return self.gl(x) & self.atoms

    def glc(self, x):
        return self.gl(x) & self.constants

    def gu(self, x):
        return OrderedSet(self.graph[x])

    def dis(self, a, b):
        return self.gla(a) - self.gla(b)

    def less(self, a, b):
        for phi in self.atoms:
            neighbors = self.graph[phi]
            if not (
                    (a not in neighbors) or (b in neighbors)
            ):
                return False
        return True

    def u(self, x):
        return OrderedSet(y for y in self.atoms + self.constants + self.terms if self.less(x, y))

    def __init__(self):
        self.img_dir = "img/master/"
        self.drawing_mapping = {"atom": "r", "constant": "g", "term": "b"}
        self.relations = data.R['+'] + data.R['-']
        self.pinning_relations = []
        self.terms = []
        self.constants = []
        self.atoms = []
        self.kinds = {"atom": self.atoms, "constant": self.constants, "term": self.terms}
        self.phi_counter = 0
        self.epsilon_counter = 0
        self.psi_counter = 0
        self.epsilon_prime_counter = 0
        self.pinning_term_counter = 0
        self.draw_flag = True

        self.fig_counter = 0

        self.graph = nx.DiGraph()
        self.add_atom(data.zero, node_code.master_atom(data.zero))
        for constant in data.constants:
            self.add_constant(constant, node_code.master_constant(constant))

        for (target, relation, term) in self.relations:
            self.add_term(term, node_code.master_term(term))

        self.add_edge(data.zero, data.target)
        for constant in self.constants:
            self.add_edge(data.zero, constant)
        for t in self.terms:
            for (n, color) in enumerate(t):
                constant = color + str(n // 2 + 1) + str(n % 2 + 1)
                self.add_edge(constant, t)

        self.close_graph()

    def draw_and_save(self):
        pyplot.clf()
        self.draw()
        pyplot.savefig(self.img_dir + str(drawing.draw_and_save_counter).zfill(4) + ".png")
        drawing.draw_and_save_counter += 1

    def draw(self, avoid_clutter=True, node_size=900, pos=None):
        types = list((nx.get_node_attributes(self.graph, "type")).values())
        colors = [self.drawing_mapping[x] for x in types]
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
        flipped_pos = {node: (x, -y) for (node, (x, y)) in pos.items()}
        return flipped_pos
