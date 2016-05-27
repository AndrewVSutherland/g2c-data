# Canonizes the number fields obtained in the endomorphism representation;
# the input and the output of these functions are all in string form, as by
# exporting from Magma to Sage we have already "broken" the link with actual
# objects.

def Minimal_Sequence(Ls):
    # Finds the sequence of minimal height in a list of sequences; if there are
    # two, then it return the first in the lexicographical ordering.
    sortf = lambda L : [ max([ QQ(c).height() for c in L ]), L[0] ]
    return sorted(Ls, key = sortf)

def Canonize_Field(frep):
    # Canonizes the number field associated to frep. Second return value gives
    # isomorphism, in the sense that it yields the roots of the canonized
    # polynomial in the original field.
    R.<x> = PolynomialRing(QQ)
    den = frep[len(frep) - 1]
    frepdiv = [ c / den for c in frep]
    f = R(frepdiv)
    fcan = R(str(gp.polredabs(f)))
    L.<r> = NumberField(f)
    froots = fcan.roots(L)
    return [ fcan.list(), [ rt[0].list() for rt in froots ] ]

def Canonize_Subfield(grep, frep):
    # Canonizes subfield given by grep in ambient frep that is supposed to have
    # been polredabs'ed already.
    R.<x> = PolynomialRing(QQ)
    f = R(frep)
    L.<r> = NumberField(f)
    tup = L.subfield(L(grep))
    K = tup[0]
    phi = tup[1]
    can = Canonize_Field(K.0.minpoly().list())
    Kcan = can[0]
    rtreps = [ phi(K(rt)).list() for rt in can[1] ]
    return [ Kcan, Minimal_Sequence(rtreps)[0] ]
