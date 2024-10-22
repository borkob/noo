# encoding=utf8


from sys import argv
from threading import Thread
from queue import Queue
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


def update_labs(x, c, i):
    k, lmt = 1, np.max([len(x) - i, i + 1])
    while k < lmt:
        ck = c[k]
        if i + k < len(x): ck -= 2 * x[i] * x[k + i]
        if k <= i: ck -= 2 * x[i - k] * x[i]
        c[k] = ck
        k += 1
    x[i] = - x[i]


def neighbor_e(x, c, i):
    r"""TODO.

    Args:
        x (list[int]): TODO.
        c (list[int]): TODO.
        i (int): TODO.

    Returns:
        int: TODO.
    """
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


def neighbor_thread_e(xq, cq, iq, rq):
    r"""Thread for calculating the energy of neighbors.

    Args:
        xq (Queue[Union[None, list[int]]]): Quene for getting the sequence.
        cq (Queue[list[int]]): Quene for getting the C of sequence.
        iq (Queue[list[int]]): Quene for getting the index of premutation of sequence.
        rq (Queue[tuple[int, int]]): Queue for return data to main thread.
    """
    while True:
        x = xq.get()
        if x == None: return
        c, inds = cq.get(), iq.get()
        rq.put([(i, neighbor_e(x, c, i)) for i in inds])


def labs_search_e_neighborhood(L, nfes=1e4, seed=1, threads_no=1):
    r"""Search for bit sequences for the sequence with the lowers energy.

    Args: 
        L (int): Length of the search sequence.
        nfes (int): Number of allowed fucntion evaluations.
        seed (int): Seed for random generator.
        threads_no (int): Number of threads to run on.

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

    xqs = [Queue(maxsize=1) for i in range(threads_no - 1)]
    cqs = [Queue(maxsize=1) for i in range(threads_no - 1)]
    iqs = [Queue(maxsize=1) for i in range(threads_no - 1)]
    rq = Queue(maxsize=threads_no - 1)
    ts = [Thread(target=neighbor_thread_e, args=(xqs[i], cqs[i], iqs[i], rq)) for i in range(threads_no - 1)]
    for t in ts: t.start()

    while n < nfes:
        best_neighbor_e, best_neighbor_i = np.Inf, -1
        
        # TODO
        for i in range(threads_no - 1):
            # TODO prepare threads work load

        # FIXME delete when proper
        for i in range(L):
            if n + 1 >= nfes: break
            else: n += 1
            e = neighbor_e(curr_x, c, i)
            if e < best_neighbor_e: best_neighbor_i, best_neighbor_e = i, e

        if best_neighbor_e >= curr_e:
            if n + 1 >= nfes: break
            else: n += 1
            curr_x = random_sol(rng, L)
            curr_e = eval_e(curr_x, c)
        else:
            update_labs(curr_x, c, best_neighbor_i)
            curr_e = best_neighbor_e
            if curr_e < best_e: best_x, best_e = np.copy(curr_x), curr_e

    for q in xqs: q.put(None)
    for t in ts: t.join()
    return best_x, best_e


if __name__ == "__main__":
    l, s, n, t = 10, 1, 1e4, 1
    if len(argv) > 1: l = int(argv[1])
    if len(argv) > 2: t = int(argv[1])
    if len(argv) > 3: s = int(argv[2])
    if len(argv) > 4: n = int(argv[3])
    print('---------- Energy ----------')
    x_e, e = labs_search_e_neighborhood(l, n, s, t)
    print('x:\n', x_e, '\nE: ', e, '\nF: ', eval_mf(x_e), '\n')
    #print('---------- PSL ----------')
    #x_psl, psl = labs_search_psl_neighborhood(l, n, s)
    #print('x:\n', x_psl, '\npsl: ', psl, '\nE: ', eval_e_x(x_psl), '\nF: ', eval_mf(x_psl))


