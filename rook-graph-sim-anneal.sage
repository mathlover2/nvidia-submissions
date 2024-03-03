# Global constants. These are set at invocation.
SIZE = 0
FACTOR = 0
ITER_PER_TEMP = 0
DIM = [10, 10]
G = Graph()

def simulated_annealing(s_0, energy, neighbor, kmax, temperature, P, verbose=False):
    s, k, e_0 = s_0, 1, energy(s_0)
    T = temperature(0)
    if verbose:
        print(k, e_0, T)
    while k < kmax:
        T = temperature(k)
        s_new = neighbor(s)
        e_new = energy(s_new)
        if e_0 == 0:
            break
        if verbose:
            if k % 1000 == 0:
                print(k, e_0, T)
        if P(e_0, e_new, T) >= random():
            s = s_new
            e_0 = e_new
        k += 1
    return s

    
def count_uncovered(subset):
    elements = subset.elements
    g = subset.g
    acc = set()
    for element in elements:
        acc.update(g.neighbor_dict[element])
    return g.order() - len(acc)

def energy_function(s):
    return count_uncovered(s)

def neighbor_function(s):
    return s.mutate()

def probability_function(e_0, e_new, T):
    if e_new <= e_0:
        return 1
    return math.exp((e_0 - e_new)/T)

def update_temperature_function(k):
    return 2*FACTOR**floor(k/ITER_PER_TEMP)


class GSubset:
    def __init__(self, elements):
        self.g = G
        self.vertices = frozenset([1..self.g.order()])
        self.elements = set(elements)
    def mutate(self):
        element_to_change = choice(list(self.elements))
        try:
            new_element = choice(list(self.g.neighbor_dict[element_to_change]
                                      - set(self.elements)))
        except IndexError:
            new_element = choice(list(self.vertices - set(self.elements)))
        new_set = self.elements - {element_to_change} | {new_element}
        return GSubset(new_set)

    
def generate_initial_state(size):
    return GSubset(sample(G.vertices(), size))

if __name__ == "__main__":
    import argparse
    import signal
    
    parser = argparse.ArgumentParser(description='Process command line arguments.')
    parser.add_argument('size', metavar='k', type=int,
                        help='Size of the covering')
    parser.add_argument('decreasefactor', metavar='chi_t', type=float,
                        help='Temperature decrease factor')
    parser.add_argument('iterationspertemp', metavar='N', type=int,
                        help='Number of iterations between application of '
                        + 'decrese factor')
    parser.add_argument('dim', metavar='d', type=int, nargs='+',
                        help='Dimensions of the rook graph')
    args = parser.parse_args()

    SIZE = args.size
    FACTOR = args.decreasefactor
    ITER_PER_TEMP = args.iterationspertemp
    DIM = args.dim
    
    G = graphs.RookGraph(DIM)
    G.relabel([1..G.order()])
    G.neighbor_dict = dict()
    for n in [1..G.order()]:
        G.neighbor_dict[n] = frozenset(G.neighbors(n, closed=True))

    def sigint_handler(signal, frame):
        exit(0)

    signal.signal(signal.SIGINT, sigint_handler)

    s = generate_initial_state(SIZE)
    res = simulated_annealing(s, energy_function, neighbor_function, Infinity,
                              update_temperature_function, probability_function, True)
    print(res.elements)
