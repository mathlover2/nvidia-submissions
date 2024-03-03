from itertools import product as iprod, combinations as icomb

class ConstraintDict(dict):
    r"""A dict object"""
    def __init__(self, network, *args, **kwds):
        dict.__init__(self, *args, **kwds)
        self.network = network

    def __getitem__(self, index):
        if index in self.keys():
            return dict.__getitem__(self, index)
        else:
            if len(index) == 2:
                i, j = index
                if (j, i) in self.keys():
                    return [(b, a) for (a, b) in dict.__getitem__(self, (j, i))]
            domains = [self.network.domains[name] for name in index]
            return list(iprod(*domains))     
        


class ConstraintNetwork(SageObject):
    r"""A constraint network object, suited for small examples."""
    def __init__(self, doms, cons):
        # `doms` is a dict of name-domain pairs. These are in the form
        # of lists. `cons` is a dict mapping domain name pairs to
        # constraint data. Constraint data can be either a list of
        # allowed tuple values or a function that, applied to elements
        # in the Cartesian product of the domains, returns True if the
        # tuple is allowed.
        self.var_names = list(doms.keys())
        self.domains = dict()
        for name in self.var_names:
            self.domains[name] = doms[name]
        self.constraints = ConstraintDict(self)
        for con in cons.keys():
            if callable(cons[con]):
                f = cons[con]
                domains = [self.domains[c] for c in con]
                self.constraints[con] = [tup for tup in iprod(*domains) if f(*tup)]
            else:
                self.constraints[con] = cons[con]
            

def revise(C, i, j):
    r"""Given the input constraint network C and names i and j, return the
domain of the variable with name i after being reduced wrt the
variable with name j.
    """
    Di, Dj, Rij = C.domains[i], C.domains[j], C.constraints[i, j]
    Di_new = []
    for ai in Di:
        if not any((ai, aj) in Rij for aj in Dj):
            pass
        else:
            Di_new.append(ai)
    return Di_new

def revise_3(C, i, j, k):
    Di, Dj, Dk = C.domains[i], C.domains[j], C.domains[k]
    Rij = C.constraints[i, j]
    Rik, Rjk = C.constraints[i, k], C.constraints[j, k]
    Rij_new = []
    for (a, b) in Rij:
        if not any((a, c) in Rik and (b, c) in Rjk for c in Dk):
            pass
        else:
            Rij_new.append((a, b))
    return Rij_new
                        
def ac_1(C, verbose=True):
    r"""Make C arc-consistent via a naive method. At the moment, this
algorithm only works for binary constraint networks.

    """
    unchanged = False
    while not unchanged:
        unchanged = True
        for (p1, p2) in C.constraints:
            for (i, j) in [(p1, p2), (p2, p1)]:
                Di_old = C.domains[i][:]
                Di_new = revise(C, i, j)
                if Di_old != Di_new:
                    if verbose:
                        print(i, ":", Di_old, " -> ", Di_new)
                    unchanged = False
                    C.domains[i] = Di_new

def pc_2(C, verbose=True):
    Q = [(i,k,j) for i in C.var_names
         for j in C.var_names for k in C.var_names
         if i < j and k != i and k != j]
    while Q:
        i, k, j = Q.pop()
        Rij_old = C.constraints[i, j][:]
        Rij_new = revise_3(C, i, j, k)
        if Rij_old != Rij_new:
            if verbose:
                print((i, j), ":", Rij_old, " -> ", Rij_new)
            C.constraints[i, j] = Rij_new
            Q.extend([(l, i, j) for l in C.var_names if l != i and l != j])
            Q.extend([(l, j, i) for l in C.var_names if l != i and l != j])
