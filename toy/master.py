import networkx as nx
from toy import data
from toy.LaTeX import node_code


class Master:
    """Master Algebra and Graph"""

    def add_node(self, kind, node, latex):
        kinds = {"atom": self.atoms, "constant": self.constants, "term": self.terms}
        assert (kind in kinds)
        node_list = kinds[kind]

        if latex is None:
            latex = node
        self.graph.add_node(node, type=kind, latex=latex)
        node_list.append(node)

    def add_atom(self, atom, latex):
        self.add_node("atom", atom, latex)

    def add_constant(self, constant, latex):
        self.add_node("constant", constant, latex)

    def add_term(self, term, latex):
        self.add_node("term", term, latex)

    def add_edge(self, a, b):
        assert (a in self.graph.nodes)
        assert (b in self.graph.nodes)
        self.graph.add_edge(a, b)

    def close_graph(self):
        self.graph = nx.transitive_closure(self.graph)

    def __init__(self):
        self.relations = data.relations
        self.terms = []
        self.constants = []
        self.atoms = []

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

    def draw(self, avoid_clutter=True, node_size=900, pos=None):
        types = list((nx.get_node_attributes(self.graph, "type")).values())
        mapping = {"atom": "r", "constant": "g", "term": "b"}
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
        # pos = nx.spring_layout(g)
        pos = nx.drawing.nx_agraph.graphviz_layout(self.graph, prog="dot")
        atoms_pos = sorted([p for (n, p) in pos.items() if n in self.atoms])
        constants_pos = sorted([p for (n, p) in pos.items() if n in self.constants])
        terms_pos = sorted([p for (n, p) in pos.items() if n in self.terms])
        pos = dict(zip(self.atoms + self.constants + self.terms,
                       atoms_pos + constants_pos + terms_pos))
        flipped_pos = {node: (x, -y) for (node, (x, y)) in pos.items()}
        return flipped_pos
