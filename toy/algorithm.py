import random

from ordered_set import OrderedSet

from toy import drawing
from toy.dual import d, und, dlatex
from toy.data import R, target, zero, zero_star, mixed, empty

aaaa = True
if aaaa:
    choice_feeder = [1, 0] + [0] * 6 + [1] * 3 + [0] * (13 + 6)
    choices_made = 0


    def random_choice(tup):
        global choices_made
        if choices_made < len(choice_feeder):
            ret = tup[choice_feeder[choices_made]]
            choices_made += 1
            return ret
        n = None
        while n is None:
            print(*[f'{i} - ' + str(option) for (i, option) in enumerate(tup)], sep='\n')
            n = int(input())
            if 0 <= n < len(tup):
                return tup[n]
            else:
                print('try again')
                n = None
else:
    random_choice = random.choice

choice = random_choice  # todo: write function choice


def same_set(ordered_set_1, ordered_set_2):
    if len(ordered_set_1) != len(ordered_set_2):
        return False
    for a in ordered_set_1:
        if a not in ordered_set_2:
            return False
    return True


def trace(master, dual, x):
    assert(x in master.graph.nodes)
    phi_set = master.gla(x)
    if phi_set:
        set_list = [dual.gla(d(phi)) for phi in phi_set]
        return OrderedSet.intersection(*set_list)
    else:
        return dual.atoms


def test_trace_constraints(master, dual):
    positive_trace_constraints = [trace(master, dual, pe).issubset(trace(master, dual, target)) for (v, plus, pe) in
                                  R['+']]
    negative_trace_constraints = [not trace(master, dual, ne).issubset(trace(master, dual, target)) for (v, minus, ne)
                                  in
                                  R['-']]
    return all(positive_trace_constraints + negative_trace_constraints)


def find_strongly_discriminant_constant(master, dual, a, b):
    omega = OrderedSet(d(c) for c in master.glc(a))
    u = trace(master, dual, b)
    while u:
        zeta = choice(tuple(u))
        print(drawing.draw_and_save_counter, zeta)
        u.remove(zeta)
        choose_from = omega - dual.gu(zeta)
        if choose_from:
            dc = choice(tuple(choose_from))
            print(drawing.draw_and_save_counter, dc)

            return und(dc)
    return None


def enforce_negative_trace_constraints(master, dual, negative_relations):
    phi_counter_before = master.phi_counter
    zeta_counter_before = dual.zeta_counter
    for (a, relation, b) in negative_relations:
        if trace(master, dual, b).issubset(trace(master, dual, a)):
            while True:
                c = find_strongly_discriminant_constant(master, dual, a, b)
                if c is None:
                    choose_from = dual.glc(d(b)) - dual.gl(d(a))
                    h = choice(tuple(choose_from))
                    print(drawing.draw_and_save_counter, h)

                    # add new atom zeta to dual and edge zeta -> h
                    dual.zeta_counter += 1
                    zeta = dual.add_atom(f"zeta_{dual.zeta_counter}", r"$\zeta_{" + f'{dual.zeta_counter}' + "}$")
                    dual.add_edge(zeta, h)
                    dual.close_graph()
                else:
                    break
            # add new atom phi to M and edge phi -> c
            master.phi_counter += 1
            phi = master.add_atom(f"phi_{master.phi_counter}", r"$\phi_{" + f'{master.phi_counter}' + "}$")
            dual.add_dual_of_atom(d(phi), dlatex(master.graph.nodes[phi]))
            master.add_edge(phi, c)
            dual.add_edge(d(c), d(phi))

            master.close_graph()
            dual.close_graph()
    phis_created = master.phi_counter - phi_counter_before
    zetas_created = dual.zeta_counter - zeta_counter_before
    return phis_created, zetas_created


def enforce_positive_trace_constraints(master, dual, positive_relations):
    epsilons_before = master.epsilon_counter
    for (D, relation, E) in positive_relations:
        while not trace(master, dual, E).issubset(trace(master, dual, D)):
            choose_from = trace(master, dual, E) - trace(master, dual, D)
            zeta = random_choice(tuple(choose_from))
            print(drawing.draw_and_save_counter, zeta)
            gamma = OrderedSet(c for c in master.glc(E) if zeta not in dual.gl(d(c)))
            if not gamma:
                dual.add_edge(zeta, d(D))
                dual.close_graph()
            else:
                c = random_choice(tuple(gamma))
                print(drawing.draw_and_save_counter, c)

                # add new atom epsilon to M and edge epsilon -> c
                master.epsilon_counter += 1
                epsilon = master.add_atom(f"eps_{master.epsilon_counter}",
                                          r"$\varepsilon_{" + f'{master.epsilon_counter}' + "}$")
                dual.add_dual_of_atom(d(epsilon), dlatex(master.graph.nodes[epsilon]))
                master.add_edge(epsilon, c)
                dual.add_edge(d(c), d(epsilon))

                master.close_graph()
                dual.close_graph()
    epsilons_created = master.epsilon_counter - epsilons_before
    return epsilons_created


def sparse_crossing(master, dual, a, b):
    psi_counter_before = master.psi_counter
    epsilon_prime_counter_before = master.epsilon_prime_counter

    A = master.gla(a) - master.gl(b)
    U = OrderedSet()
    for phi in A:
        U, B, delta = OrderedSet(), master.gla(b), dual.atoms - dual.gl(d(phi))
        while True:
            eps = random_choice(tuple(B))
            print(drawing.draw_and_save_counter, eps)
            delta_prime = delta & dual.gl(d(eps))
            if (not delta) or (not same_set(delta, delta_prime)):
                # add new atom psi to M and edges psi -> phi and psi -> epsilon
                master.psi_counter += 1
                psi = master.add_atom(f"psi_{master.psi_counter}",
                                      r"$\psi_{" + f'{master.psi_counter}' + "}$")
                dual.add_dual_of_atom(d(psi), dlatex(master.graph.nodes[psi]))
                master.add_edge(psi, phi)
                master.add_edge(psi, eps)
                dual.add_edge(d(phi), d(psi))
                dual.add_edge(d(eps), d(psi))

                master.close_graph()
                dual.close_graph()

                delta = delta_prime
                U.append(eps)
            B.remove(eps)
            if not delta:
                break

    for eps in U:
        # create new atom eps_prime and edge eps_prime -> eps
        master.epsilon_prime_counter += 1
        eps_prime = master.add_atom(f"eps_prime_{master.epsilon_prime_counter}",
                                    r"$\varepsilon_{" + f'{master.epsilon_prime_counter}' + "}' $")
        dual.add_dual_of_atom(d(eps_prime), dlatex(master.graph.nodes[eps_prime]))
        master.add_edge(eps_prime, eps)
        dual.add_edge(d(eps), d(eps_prime))

        master.close_graph()
        dual.close_graph()

    to_remove = list(U | A)
    master.remove_atoms_from(to_remove)
    dual.remove_dual_of_atoms_from(map(d, to_remove))
    master.close_graph()
    dual.close_graph()

    psi_added = master.psi_counter - psi_counter_before
    eps_prime_added = master.epsilon_prime_counter - epsilon_prime_counter_before
    atoms_removed = len(to_remove)

    return psi_added, eps_prime_added, atoms_removed


def atom_set_reduction(master, dual, keep_zero=False):
    # this is a stochastic algorithm to reduce the number of atoms.
    # this method can be called anytime, and should reduce the number of atoms in a few calls
    # Q, A = OrderedSet(), master.constants  # todo: verify which line is correct, this one or the line below
    Q, A = OrderedSet(), [c for c in master.constants if trace(master, dual, c) is not None]
    if keep_zero:
        Q.append(zero)
    while True:
        c = random_choice(tuple(A))
        print(drawing.draw_and_save_counter, c)
        A.remove(c)
        S = master.gl(c) & Q
        if not S:
            W = dual.atoms
        else:
            W = OrderedSet.intersection(*[dual.gla(d(phi)) for phi in S])
        PHI = OrderedSet(d(phi) for phi in master.gla(c))
        tr = trace(master, dual, c)
        while (tr is not None) and (not same_set(W, tr)):
            choose_from = W - tr
            xi = random_choice(tuple(choose_from))
            print(drawing.draw_and_save_counter, xi)
            choose_from = PHI - dual.gu(xi)
            phi = und(random_choice(tuple(choose_from)))
            print(drawing.draw_and_save_counter, phi)
            Q.append(phi)
            W &= dual.gla(d(phi))
        if not A:
            break

    to_remove = master.atoms - Q
    master.remove_atoms_from(to_remove)
    dual.remove_dual_of_atoms_from(map(d, to_remove))

    master.close_graph()
    dual.close_graph()
    return len(to_remove)


def atom_set_reduction_for_the_dual_algebra(dual, keep_zero=False):
    Q = OrderedSet()
    if keep_zero:
        Q.append(zero_star)
    S = R['-'].copy()
    while S:
        r = random_choice(S)
        print(drawing.draw_and_save_counter, r)
        S.remove(r)
        a, relation, b = r
        dis = dual.dis(d(b), d(a))
        if not (dis & Q):
            xi = choice(tuple(dis))
            print(drawing.draw_and_save_counter, xi)
            Q.append(xi)
    to_remove = dual.atoms - Q
    dual.remove_atoms_from(to_remove)
    dual.close_graph()


def combine_from(master, dual, combine_list):
    target_present = target in combine_list
    CL = combine_list.copy()
    if target_present:
        CL.remove(target)
    assert (len(CL) > 0)
    while len(CL) > 1:
        a, b = CL[0], CL[1]
        c = combine(master, a, b)
        CL = CL[1:]
        CL[0] = c
    ret = CL[0]
    if ret not in master.graph.nodes:
        master.combined_terms_counter += 1
        master.add_term(ret, "$T_{" + f"{master.combined_terms_counter}" + "}$")
        dual.add_constant(d(ret), dlatex(master.graph.nodes[ret]))
        for x in combine_list:
            master.add_edge(x, ret)
            dual.add_edge(d(ret), d(x))
        master.close_graph()
        dual.close_graph()
    return ret


def combine(master, a, b):
    assert all((x in master.constants or len(x) == 4) for x in (a, b))
    assert all(x != target for x in (a, b))
    # first, convert a and b to terms
    r = [a, b]
    for i in range(2):
        if r[i] in master.constants:
            T = [empty] * 4
            index = (int(r[i][1]) - 1) * 2 + (int(r[i][2]) - 1)
            T[index] = r[i][0]
            r[i] = T
    a, b = r
    return combine_terms(a, b)


def combine_terms(a, b):
    T = [empty] * 4
    for i in range(4):
        if a[i] == empty:
            T[i] = b[i]
        elif b[i] == empty:
            T[i] = a[i]
        elif a[i] == mixed:
            T[i] = mixed
        elif b[i] == mixed:
            T[i] = mixed
        elif a[i] == b[i]:
            T[i] = a[i]
        else:
            T[i] = mixed
    print(*a[0:2], ' ', *b[0:2], ' ', *T[0:2])
    print(*a[2:4], '+', *b[2:4], '=', *T[2:4], end='\n\n')
    return tuple(T)


def generation_of_pinning_terms_and_relations(master, dual):
    for phi in master.atoms:
        H = master.constants - master.u(phi)

        if H:
            T_phi = combine_from(master, dual, list(H))
            master.pinning_term_counter += 1
            old_name = T_phi
            new_name = 'T_' + phi
            master.rename_node(old_name, new_name,
                               '$T_{' + f'{(master.graph.nodes[phi]["latex"])[1:-1]}' + '}$')
            dual.rename_node(d(old_name), d(new_name), dlatex(master.graph.nodes[new_name]))
            T_phi = new_name
            for c in master.constants & master.u(phi):
                r = (c, '-', T_phi)
                master.pinning_relations.append(r)


def is_consistent(master, dual, input_order_relations):
    mentioned_terms = [b for (a, rel, b) in input_order_relations if b in master.terms]
    consistent = True
    for T1 in mentioned_terms:
        for T2 in mentioned_terms:
            if T1 != T2:
                constants1 = master.glc(T1)
                constants2 = master.glc(T2)
                if constants1.issubset(constants2):
                    consistent = (d(T2), d(T1)) in dual.graph.edges
                    if not consistent:
                        return False
    return consistent
