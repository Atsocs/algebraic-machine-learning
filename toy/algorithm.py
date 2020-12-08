import random

from ordered_set import OrderedSet

from toy import drawing
from toy.dual import d, und, dlatex
from toy.data import positive_examples, negative_examples, target


def same_set(ordered_set_1, ordered_set_2):
    if len(ordered_set_1) != len(ordered_set_2):
        return False
    for a in ordered_set_1:
        if a not in ordered_set_2:
            return False
    return True


def trace(master, dual, x):
    # x is always a member of master
    phi_set = master.gla(x)
    if phi_set:
        set_list = [dual.gla(d(phi)) for phi in phi_set]
        return OrderedSet.intersection(*set_list)
    else:
        return dual.atoms
        # raise Exception("trace(x): x contains no atoms")


def test_trace_constraints(master, dual):
    positive_trace_constraints = [trace(master, dual, pe).issubset(trace(master, dual, target)) for pe in
                                  positive_examples]
    negative_trace_constraints = [not trace(master, dual, ne).issubset(trace(master, dual, target)) for ne in
                                  negative_examples]
    return all(positive_trace_constraints + negative_trace_constraints)


def find_strongly_discriminant_constant(master, dual, a, b):
    omega = OrderedSet(d(c) for c in master.glc(a))
    u = trace(master, dual, b)
    while u:
        zeta = random.choice(tuple(u))  # todo: no need of random here
        print(drawing.draw_and_save_counter, zeta)
        u.remove(zeta)
        choose_from = omega - dual.gu(zeta)
        if choose_from:
            dc = random.choice(tuple(choose_from))  # todo: no need of random here
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
                    h = random.choice(tuple(choose_from))  # todo: maybe not random here is okay
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
            zeta = random.choice(tuple(choose_from))
            print(drawing.draw_and_save_counter, zeta)
            gamma = OrderedSet(c for c in master.glc(E) if zeta not in dual.gl(d(c)))
            if not gamma:
                dual.add_edge(zeta, d(D))
                dual.close_graph()
            else:
                c = random.choice(tuple(gamma))
                print(drawing.draw_and_save_counter, c)

                # add new atom epsilon to M and edge epsilon -> c #todo: maybe change epislon to phi too?
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
            eps = random.choice(tuple(B))
            print(drawing.draw_and_save_counter, eps)
            delta_prime = delta & dual.gl(d(eps))
            if (not delta) or (not same_set(delta, delta_prime)):
                # add new atom psi to M and edges psi -> phi and psi -> epsilon #todo: maybe change epislon to phi too?
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


def atom_set_reduction(master, dual):
    # this is a stochastic algorithm to reduce the number of atoms.
    # this method can be called anytime, and should reduce the number of atoms in a few calls
    # Q, A = OrderedSet(), master.constants
    Q, A = OrderedSet(), [c for c in master.constants if trace(master, dual, c) is not None]
    while True:
        c = random.choice(tuple(A))
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
            xi = random.choice(tuple(choose_from))
            print(drawing.draw_and_save_counter, xi)
            choose_from = PHI - dual.gu(xi)
            phi = und(random.choice(tuple(choose_from)))
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


def atom_set_reduction_for_the_dual_algebra(dual):
    Q = OrderedSet()
    S = [(ne, target) for ne in negative_examples]
    while S:
        r = random.choice(S)
        print(drawing.draw_and_save_counter, r)
        S.remove(r)
        a, b = r
        dis = dual.dis(d(b), d(a))
        if dis and not (dis & Q):
            xi = random.choice(tuple(dis))  # todo no need of random
            print(drawing.draw_and_save_counter, xi)
            Q.append(xi)
    to_remove = dual.atoms - Q
    dual.remove_atoms_from(to_remove)
    dual.close_graph()


def merge_into(master, dual, merge_set, x):
    to_union = [master.gla(x) for x in merge_set]
    gla = OrderedSet().union(*to_union)
    for atom in gla:
        master.add_edge(atom, x)
        dual.add_edge(d(x), d(atom))
    master.close_graph()
    return x


def generation_of_pinning_terms_and_relations(master, dual):
    for phi in master.atoms:
        H = master.constants - master.u(phi)
        master.pinning_term_counter += 1
        T = master.add_term('PT_' + str(master.pinning_term_counter),
                            '$PT_{' + f'{master.pinning_term_counter}' + '}$')
        dual.add_constant(d(T), dlatex(master.graph.nodes[T]))
        merge_into(master, dual, H, T)
        for c in master.constants & master.u(phi):
            r = (c, '-', T)
            master.pinning_relations.append(r)
