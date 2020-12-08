import random

import numpy

import toy.dual
import toy.master
from toy import algorithm, drawing
from toy.drawing import draw, erase_img
from toy.data import R

erase_img()

# def main():
# seed = 123456
seed = 234256
numpy.random.seed(seed)
random.seed(seed)

master = toy.master.Master()
dual = toy.dual.Dual(master)
dual.reverted_negative_relations()
# TODO: corrigir visualização
# TODO: check for consistency after dual.reverted_negative_relations()
drawing.draw_and_save_counter += 1000
algorithm.enforce_negative_trace_constraints(master, dual, R['-'])
drawing.draw_and_save_counter += 1000
algorithm.enforce_positive_trace_constraints(master, dual, R['+'])

if algorithm.test_trace_constraints(master, dual):
    print("Trace constraints fufilled 1")
else:
    print("ERROR: Trace constraints NOT fufilled 1")

drawing.draw_and_save_counter += 1000
for (a, relation, b) in R['+']:
    algorithm.sparse_crossing(master, dual, a, b)

if algorithm.test_trace_constraints(master, dual):
    print("Trace constraints fufilled 2")
else:
    print("ERROR: Trace constraints NOT fufilled 2")

drawing.draw_and_save_counter += 1000
algorithm.atom_set_reduction(master, dual)
drawing.draw_and_save_counter += 1000
algorithm.atom_set_reduction_for_the_dual_algebra(dual)
drawing.draw_and_save_counter += 1000
algorithm.generation_of_pinning_terms_and_relations(master, dual)

draw(master, dual)


# if __name__ == '__main__':
# main()
