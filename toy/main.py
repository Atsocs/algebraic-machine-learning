from matplotlib import pyplot

import toy.master
import toy.dual
from toy import algorithm

master = toy.master.Master()
dual = toy.dual.Dual(master)

dual.reverted_negative_relations()
# TODO: corrigir visualização das reverted negative relations
# TODO: check for consistency after dual.reverted_negative_relations()

negative_relations = [(toy.data.target, ne) for ne in toy.data.negative_examples]

algorithm.enforce_negative_trace_constraints(master, dual, negative_relations)


def draw(save_name=None):
    fig = pyplot.figure(figsize=(8, 10))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    ax1.title.set_text('Master')
    ax2.title.set_text('Dual')
    pyplot.subplot(ax1)
    master.draw(node_size=600)
    pyplot.subplot(ax2)
    dual.draw(node_size=600)
    if save_name is None:
        pyplot.show()
    else:
        pyplot.savefig(save_name)
