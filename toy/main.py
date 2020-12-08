import random

import numpy

import toy.dual
import toy.master
from toy import algorithm
from toy.drawing import draw
from toy.data import R


def main():
    seed = 123456
    numpy.random.seed(seed)
    random.seed(seed)

    master = toy.master.Master()
    dual = toy.dual.Dual(master)
    dual.reverted_negative_relations()
    # TODO: corrigir visualização
    # TODO: check for consistency after dual.reverted_negative_relations()
    algorithm.enforce_negative_trace_constraints(master, dual, R['-'])
    algorithm.enforce_positive_trace_constraints(master, dual, R['+'])
    for pr in R['+']:
        algorithm.sparse_crossing(master, dual, *pr)
    algorithm.atom_set_reduction(master, dual)

    draw(master, dual)


if __name__ == '__main__':
    main()
