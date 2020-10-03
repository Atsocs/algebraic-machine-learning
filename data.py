from pprint import pformat

zero = "0"
zero_star = "0*"
target = "v"


def get_constants(expand=False):
    if not expand:
        return [let + str(i) for let in ('B', 'W') for i in range(4)]

    no_colors = ['N', 'N', 'N', 'N']
    constants = []
    for color in ['B', 'W']:
        for i in range(4):
            c = no_colors.copy()
            c[i] = color
            constants.append((color + str(i), c))
    return constants


def get_terms(expand=False):
    pos = [('B', 'B', 'W', 'W'),
           ('W', 'W', 'B', 'B')]
    neg = [('W', 'B', 'W', 'B'),
           ('W', 'W', 'B', 'W'),
           ('W', 'W', 'W', 'B')]

    pos_names = [r"$T_{" + str(i+1) + r"}^{+}$" for i in range(len(pos))]
    neg_names = [r"$T_{" + str(i+1) + r"}^{-}$" for i in range(len(neg))]
    if not expand:
        return pos_names, neg_names

    pos_def, neg_def = [], []
    for (lst, lst_def) in ((pos, pos_def), (neg, neg_def)):
        for term in lst:
            term_def = []
            for i in range(4):
                color = term[i]
                term_def.append(color + str(i))
            lst_def.append(term_def)
    pos_def = list(zip(pos_names, pos_def))
    neg_def = list(zip(neg_names, neg_def))
    return pos_def, neg_def


def main():
    constants = get_constants(expand=True)
    pos_def, neg_def = get_terms(expand=True)
    print("constants = " + pformat(constants))
    print("pos_def = " + pformat(pos_def))
    print("neg_def = " + pformat(neg_def))


if __name__ == '__main__':
    main()
