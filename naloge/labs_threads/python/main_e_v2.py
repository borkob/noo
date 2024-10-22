# encoding=utf8


from sys import argv
import time
from multiprocessing import Queue, Process
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


def update_e(x, c, i):
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
        x (list[int]): TODO
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


def labs_search_thread_e_neighborhood(L, rq, nfes=1e4, seed=1):
    r"""Search for bit sequences for the sequence with the lowers energy.

    Args: 
        L (int): Length of the search sequence.
        rq (Queue): return Queue.
        nfes (int): Number of allowed fucntion evaluations.
        seed (int): Seed for random generator.
    """
    rng = np.random.RandomState(seed)
    best_x = random_sol(rng, L)
    c = np.zeros(L)
    best_e = eval_e(best_x, c)
    n, curr_x, curr_e = 1, np.copy(best_x), best_e
    while n < nfes:
        best_neighbor_e, best_neighbor_i = np.Inf, -1
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
            update_e(curr_x, c, best_neighbor_i)
            curr_e = best_neighbor_e
            if curr_e < best_e: best_x, best_e = np.copy(curr_x), curr_e
    rq.put((best_x, best_e), block=False)


if __name__ == "__main__":
    l, s, n, ts_no = 100, 1, 1e5, 4 
    if len(argv) > 1: l = int(argv[1])
    if len(argv) > 2: ts_no = int(argv[2])
    if len(argv) > 3: s = int(argv[3])
    if len(argv) > 4: n = int(argv[4])
    nfes_t = n // ts_no
    print('Number of used threads %d.\nNumber of evaluations per thread %d.\n' % (ts_no, nfes_t))
    q = Queue(maxsize=ts_no)
    ts = [Process(target=labs_search_thread_e_neighborhood, args=(l, q, nfes_t, s + i + 1)) for i in range(ts_no - 1)]
    print('Threads createed')
    t_start = time.time()
    for t in ts: t.start()
    print('Treads started witout main thread')
    labs_search_thread_e_neighborhood(l, q, nfes_t, s)
    for t in ts: t.join()
    t_end = time.time()
    print('Threads ended in %f s.\n\nResults from threads:' % (t_end - t_start))
    best = (None, np.Inf)
    for i in range(ts_no):
        r = q.get()
        print(r)
        if r[1] < best[1]: best = r
    print('\nGlobal best found with energy %d\nSolutions:\n%s' % (best[1], best[0]))

