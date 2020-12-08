import random
from pprint import pformat

import numpy
from matplotlib import pyplot

import toy.master
import toy.dual
from toy import algorithm
from toy.data import target, positive_examples, negative_examples

seed = 123456
numpy.random.seed(seed)
random.seed(seed)

master = toy.master.Master()
dual = toy.dual.Dual(master)

dual.reverted_negative_relations()
# TODO: corrigir visualização
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


master.fig_counter += 1000
dual.fig_counter += 1000

ctraces_before = [algorithm.trace(master, dual, c) for c in master.constants]
algorithm.atom_set_reduction(master, dual)
ctraces_after = [algorithm.trace(master, dual, c) for c in master.constants]
print(ctraces_after == ctraces_before)

# algorithm.atom_set_reduction_for_the_dual_algebra(dual)


