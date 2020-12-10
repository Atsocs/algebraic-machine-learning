import random
import numpy

import toy.dual
import toy.master
from toy import algorithm, drawing
from toy.drawing import draw, erase_img


def main():
    if drawing.draw_flag:
        erase_img()

    seed = 0
    numpy.random.seed(seed)
    random.seed(seed)

    master = toy.master.Master()
    dual = toy.dual.Dual(master)
    master, dual = batch(master, dual)

    print("Goal reached = " + str(algorithm.test_goal(master, master.relations)))

    draw(master, dual)

    return master, dual


def batch(master, dual):
    print("epoch " + str(master.epoch))

    dual.reverted_negative_relations(master.negative_relations)
    print("Consistency = " + str(algorithm.is_consistent(master, dual, master.relations)))
    drawing.draw_and_save_counter += 1000
    algorithm.enforce_negative_trace_constraints(master, dual, master.negative_relations)
    drawing.draw_and_save_counter += 1000
    algorithm.enforce_positive_trace_constraints(master, dual, master.positive_relations)

    if algorithm.test_trace_constraints(master, dual, master.relations):
        print("Trace constraints fufilled 1")
    else:
        print("ERROR: Trace constraints NOT fufilled 1")

    drawing.draw_and_save_counter += 1000
    for (a, relation, b) in master.positive_relations:
        algorithm.sparse_crossing(master, dual, a, b)

    if algorithm.test_trace_constraints(master, dual, master.relations):
        print("Trace constraints fufilled 2")
    else:
        print("ERROR: Trace constraints NOT fufilled 2")

    drawing.draw_and_save_counter += 1000
    algorithm.atom_set_reduction(master, dual)
    drawing.draw_and_save_counter += 1000
    algorithm.atom_set_reduction_for_the_dual_algebra(dual, master.negative_relations)
    drawing.draw_and_save_counter += 1000
    algorithm.generation_of_pinning_terms_and_relations(master, dual)

    master.negative_relations += master.pinning_relations
    master.relations += master.pinning_relations
    master.epoch += 1
    return master, dual


if __name__ == '__main__':
    algebras = main()
