import random

myRandom = random.Random(0)


def random_pop(set_s):
    x = random_choice(set_s)
    set_s.remove(x)
    return x, set_s


def random_choice(set_s):
    # todo: use random module and random.seed(0) to guarantee the same result always
    seq = list(set_s)
    seq.sort()
    return seq[0]


def f():
    global s
    x, s = random_pop(s)
    print(x, s)


s = set(range(10))

for i in range(10):
    f()
