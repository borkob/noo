# encoding=utf8


import time
from sys import argv
import numpy as np


def random_sol(rng, size):
    return rng.choice([-1, 1], size, p=[0.5, 0.5])


def eval_e_x(x):
    e = 0
    for k in range(1, len(x)):
        c = 0
        for i in range(len(x) - k): c += x[i] * x[i + k]
        e += c * c
    return e


def eval_mf(x):
    return len(x) ** 2 / (2 * eval_e_x(x))


def eval_e(x, c):
    e = 0
    for k in range(1, len(x)):
        c[k] = 0
        for i in range(len(x) - k): c[k] += x[i] * x[i + k]
        e += c[k] * c[k]
    return e


def neighbor_e(x, c, i):
    e, k, lmt = 0, 1, np.max([len(x) - i, i + 1])
    while k < lmt:
        ck = c[k]
        if i + k < len(x): ck -= 2 * x[i] * x[k + i]
        if k <= i: ck -= 2 * x[i - k] * x[i]
        e += ck * ck
        k += 1
    while k < len(x):
        e += c[k] * c[k]
        k += 1
    return e


def update_labs(x, c, i):
    k, lmt = 1, np.max([len(x) - i, i + 1])
    while k < lmt:
        ck = c[k]
        if i + k < len(x): ck -= 2 * x[i] * x[k + i]
        if k <= i: ck -= 2 * x[i - k] * x[i]
        c[k] = ck
        k += 1
    x[i] = - x[i]


def labs_search_e_neighborhood(L, nfes=1e4, seed=1):
    """Search for bit sequences for the sequence with the lowers energy.

    Args: 
        L (int): Length of the search sequence
        nfes (int): Number of allowed fucntion evaluations.
        seed (int): Seed for random generator

    Returns:
        Tuple[np.ndarray, float]:
            1. Best sulution found
            2. Best solution energy
    """
    rng = np.random.RandomState(seed)
    best_x = random_sol(rng, L)
    c = np.zeros(L)
    best_e = eval_e(best_x, c)
    n, curr_x, curr_e = 1, np.copy(best_x), best_e
    while n < nfes:
        best_neighbor_e, best_neighbor_i = np.Inf, -1
        for i in range(L):
            if n + 1 >= nfes: return best_x, best_e
            else: n += 1
            e = neighbor_e(curr_x, c, i)
            if e < best_neighbor_e: best_neighbor_i, best_neighbor_e = i, e
        if best_neighbor_e >= curr_e:
            if n + 1 >= nfes: return best_x, best_e
            else: n += 1
            curr_x = random_sol(rng, L)
            curr_e = eval_e(curr_x, c)
        else:
            update_labs(curr_x, c, best_neighbor_i)
            curr_e = best_neighbor_e
            if curr_e < best_e: best_x, best_e = np.copy(curr_x), curr_e
    return best_x, best_e


def eval_psl_x(x):
    psl = 0
    for k in range(1, len(x)):
        c = 0
        for i in range(len(x) - k): c += x[i] * x[i + k]
        if np.abs(c) > psl: psl = np.abs(c)
    return psl


def eval_psl(x, c):
    psl = 0
    for k in range(1, len(x)):
        c[k] = 0
        for i in range(len(x) - k): c[k] += x[i] * x[i + k]
        if np.abs(c[k]) > psl: psl = np.abs(c[k])
    return psl


def neighbor_psl(x, c, i):
    psl, k, lmt = 0, 1, np.max([len(x) - i, i + 1])
    while k < lmt:
        ck = c[k]
        if i + k < len(x): ck -= 2 * x[i] * x[k + i]
        if k <= i: ck -= 2 * x[i - k] * x[i]
        if np.abs(ck) > psl: psl = np.abs(ck)
        k += 1
    while k < len(x):
        if np.abs(c[k]) > psl: psl = np.abs(c[k])
        k += 1
    return psl


def labs_search_psl_neighborhood(L, nfes=1e4, seed=1):
    """Search for bit sequences for the sequence with the highest PSL.

    Args: 
        L (int): Length of the search sequence
        nfes (int): Number of allowed fucntion evaluations.
        seed (int): Seed for random generator

    Returns:
        Tuple[np.ndarray, float]:
            1. Best sulution found
            2. Best solution energy
    """
    rng = np.random.RandomState(seed)
    best_x = random_sol(rng, L)
    c = np.zeros(L)
    best_psl = eval_psl(best_x, c)
    n, curr_x, curr_psl = 1, np.copy(best_x), best_psl
    while n < nfes:
        best_neighbor_psl, best_neighbor_i = np.Inf, -1
        for i in range(L):
            if n + 1 >= nfes: return best_x, best_psl
            else: n += 1
            e = neighbor_psl(curr_x, c, i)
            if e < best_neighbor_psl: best_neighbor_i, best_neighbor_psl = i, e
        if best_neighbor_psl >= curr_psl:
            if n + 1 >= nfes: return best_x, best_psl
            else: n += 1
            curr_x = random_sol(rng, L)
            curr_psl = eval_psl(curr_x, c)
        else:
            update_labs(curr_x, c, best_neighbor_i)
            curr_psl = best_neighbor_psl
            if curr_psl < best_psl: best_x, best_psl = np.copy(curr_x), curr_psl
    return best_x, best_psl


if __name__ == "__main__":
    l, s, n = 100, 1, 1e5
    if len(argv) > 1: l = int(argv[1])
    if len(argv) > 2: s = int(argv[2])
    if len(argv) > 3: n = int(argv[3])
    print('---------- Energy ----------')
    t_start_e = time.time()
    x_e, e = labs_search_e_neighborhood(l, n, s)
    t_end_e = time.time()
    print('x:\n', x_e, '\nE: ', e, '\nF: ', eval_mf(x_e), '\nIn %f s.' % (t_end_e - t_start_e), '\n')
    print('---------- PSL ----------')
    t_start_psl = time.time()
    x_psl, psl = labs_search_psl_neighborhood(l, n, s)
    t_end_psl = time.time()
    print('x:\n', x_psl, '\npsl: ', psl, '\nE: ', eval_e_x(x_psl), '\nF: ', eval_mf(x_psl), '\nIn %f s.' % (t_end_psl - t_start_psl))


