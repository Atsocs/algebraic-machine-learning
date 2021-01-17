import networkx as nx

from toy.master import Master
from toy import data, drawing
from toy.LaTeX import node_code


def d(x):
    return '[' + str(x) + ']'


def und(x):
    return x[1:-1]


def dlatex(node_attr):
    return '$' + d('{' + node_attr['latex'][1: -1] + '}') + '$'


def dtype(node_attr):
    node_type = node_attr['type']
    if node_type == "atom":
        return "dual-of-atom"
    elif node_type == "term" or node_type == "constant":
        return "constant"
    else:
        assert False


class Dual:
    """Dual Algebra and Graph"""

    add_node = Master.__dict__['add_node']
    rename_node = Master.__dict__['rename_node']

    def add_dual_of_atom(self, dual_of_atom, latex):
        return self.add_node("dual-of-atom", dual_of_atom, latex)

    def add_constant(self, constant, latex):
        return self.add_node("constant", constant, latex)

    add_atom = Master.__dict__['add_atom']
    add_edge = Master.__dict__['add_edge']
    close_graph = Master.__dict__['close_graph']
    draw = Master.__dict__['draw']
    get_pos = Master.__dict__['get_pos']

    def reverted_negative_relations(self, negative_relations):
        negative_examples = [ne for (v, minus, ne) in negative_relations]
        for ne in negative_examples:
            self.zeta_counter += 1
            zeta = self.add_atom(f"zeta_{self.zeta_counter}", r"$\zeta_{" + f'{self.zeta_counter}' + "}$")
            self.add_edge(zeta, d(ne))
        self.close_graph()

    gl = Master.__dict__['gl']
    gla = Master.__dict__['gla']
    glc = Master.__dict__['glc']
    gu = Master.__dict__['gu']
    dis = Master.__dict__['dis']

    def remove_dual_of_atom(self, dual_of_atom):
        self.dual_of_atoms.remove(dual_of_atom)
        self.graph.remove_node(dual_of_atom)
        if self.draw_flag:
            self.draw_and_save()
        return True

    def remove_dual_of_atoms_from(self, dual_of_atom_list):
        removed = 0
        for dual_of_atom in dual_of_atom_list:
            removed += self.remove_dual_of_atom(dual_of_atom)
        return removed

    remove_atom = Master.__dict__['remove_atom']
    remove_atoms_from = Master.__dict__['remove_atoms_from']

    def __init__(self, master):
        self.img_dir = drawing.path + "/img/dual/"
        self.draw_flag = drawing.draw_flag
        self.fig_counter = 0
        self.drawing_mapping = {"atom": "c", "constant": "m", "dual-of-atom": "y"}
        self.constants = []  # duals of master's constants or terms
        self.dual_of_atoms = []  # duals of master's atoms
        self.atoms = []  # not duals of any element in the master
        self.kinds = {"atom": self.atoms, "constant": self.constants, "dual-of-atom": self.dual_of_atoms}

        self.constants += [d(c) for c in master.constants]
        self.constants += [d(t) for t in master.terms]
        self.dual_of_atoms += [d(a) for a in master.atoms]
        self.zeta_counter = 0

        g = master.graph.reverse()
        g = nx.relabel_nodes(g, dict(zip(g.nodes, map(d, g.nodes))))
        for n in g.nodes:
            g.nodes[n]['type'] = dtype(g.nodes[n])
            g.nodes[n]['latex'] = dlatex(g.nodes[n])
        self.graph = g

        self.add_atom(data.zero_star, node_code.dual_atom(data.zero_star))
        for t in master.terms:
            self.add_edge(data.zero_star, d(t))
        positive_examples = [pe for (v, plus, pe) in data.R['+']]
        for pe in positive_examples:
            self.add_edge(d(pe), d(data.target))

        self.close_graph()

    draw_and_save = Master.__dict__['draw_and_save']
