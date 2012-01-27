import itertools
import math
import random
import cProfile


def param_to_simplex(nn, kk):
    # least sig first
    def to_binary(n):
        xs = []
        for i in range(max(nn, kk)):
            xs.append(n % 2)
            n = n  // 2
        return xs
    vs = [to_binary(n) for n in range(2**nn)]
    hs = [to_binary(n) for n in range(2**kk)]
    def _param_to_simplex(w, b, c):
        def energy(v, h):
            e = 0
            e += sum(w[i][j] 
                    for (i, j) 
                    in itertools.product(range(nn), range(kk))
                    if v[i] and h[j]
                    )
            e += sum(b[i] for i in range(nn) if v[i])
            e += sum(c[j] for j in range(kk) if h[j])
            return e
        def psi(v, h):
            return math.exp(energy(v, h))
        z = sum(psi(v, h) for v in vs for h in hs)
        return [sum(psi(v, h) for h in hs) / z for v in vs]
    return _param_to_simplex

def rand_param(nn, kk):
    f = lambda: random.uniform(-3, 3)
    w = [[f() for j in range(kk)] for i in range(nn)]
    b = [f() for i in range(nn)]
    c = [f() for j in range(kk)]
    return (w, b, c)

def sample(nn, kk):
    f = param_to_simplex(nn, kk)
    def _sample(n):
        xs = []
        for i in range(n):
            p = rand_param(nn, kk)
            xs.append(f(*p))
        return xs
    return _sample


def total_variation_distance(p, q):
    return 0.5 * sum(abs(x - y) for (x, y) in zip(p, q))

def hausdorff(f, xs, ys):
    a = max(min(f(x, y) for y in ys) for x in xs)
    b = max(min(f(x, y) for x in xs) for y in ys)
    return max(a, b)

def main():
    kks = range(0, 4)
    ss = [sample(3, 2)(1000) for kk in kks]

    for i in kks[:-1]:
        for j in kks[i+1:]:
            d = hausdorff(total_variation_distance, ss[i], ss[j])
            print("The hausdorff distance between kk=%d and kk=%d is %f" % (i, j, d))

main()
# cProfile.run('main()')
