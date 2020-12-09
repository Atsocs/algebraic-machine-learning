import toy.data as data


def master_constant(c):
    if c is not data.target:
        latex = '$' + c[0] + '_{' + c[1:] + '}$'
    else:
        latex = "$v$"
    return latex


def master_term(t):
    # TODO: print the term's image in latex
    positive_examples = [pe for (v, plus, pe) in data.R['+']]
    negative_examples = [ne for (v, minus, ne) in data.R['-']]
    if t in positive_examples:
        i = 1 + positive_examples.index(t)
        return "$T_{" + str(i) + "}^{+}$"
    elif t in negative_examples:
        i = 1 + negative_examples.index(t)
        return "$T_{" + str(i) + "}^{-}$"
    else:
        assert False


def master_atom(a):
    # TODO
    return "$" + str(a) + "$"


def dual_atom(a):
    # TODO
    return "$" + str(a) + "$"
