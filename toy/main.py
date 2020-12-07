from matplotlib import pyplot

import toy.master
import toy.dual
from toy import algorithm
from toy.data import target, positive_examples, negative_examples

master = toy.master.Master()
dual = toy.dual.Dual(master)

dual.reverted_negative_relations()
# TODO: corrigir visualização das reverted negative relations
# TODO: check for consistency after dual.reverted_negative_relations()

negative_relations = [(target, ne) for ne in negative_examples]
positive_relations = [(target, pe) for pe in positive_examples]

algorithm.enforce_negative_trace_constraints(master, dual, negative_relations)
algorithm.enforce_positive_trace_constraints(master, dual, positive_relations)

if algorithm.test_trace_constraints(master, dual):
    print("Trace constraints fufilled")
else:
    print("ERROR: Trace constraints NOT fufilled")

for pe in positive_examples:
    algorithm.sparse_crossing(master, dual, target, pe)


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
