zero = "0"
zero_star = "0*"
target = "v"

black = "B"
white = "W"

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

negative_relations = [(target, ne) for ne in negative_examples]
positive_relations = [(target, pe) for pe in positive_examples]

R = {'-': negative_relations, '+': positive_relations}

# print("Welcome to AML Toy!".center(50, "="))
# print("constants: {}".format(constants),
#       "positive_examples: {}".format(positive_examples),
#       "negative examples: {}".format(negative_examples),
#       "R: {}".format(R), sep="\n")
# print("".center(50, "="))
