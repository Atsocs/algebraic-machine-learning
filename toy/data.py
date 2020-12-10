zero = "0"
zero_star = "0*"
target = "v"

black = "B"
white = "W"
mixed = "M"
empty = "E"

constants = [color + str(i) + str(j) for color in [black, white] for i in [1, 2] for j in [1, 2]]
constants = [target] + constants
constants = tuple(constants)

positive_examples = (
    (black, white,
     black, white),
    (white, black,
     white, black),
)

negative_examples = (
    (black, white,
     white, black),
    (white, black,
     white, white),
    (white, white,
     white, black),
)

negative_relations = [(target, '-', ne) for ne in negative_examples]
positive_relations = [(target, '+', pe) for pe in positive_examples]

R = {'-': negative_relations, '+': positive_relations}

del(positive_examples, negative_examples, negative_relations, positive_relations)
