def Elliptic_Curve_From_Lattice(L, frep, h, prec = prec0, epscomp = epscomp0,
        epsLLL = epsLLL0): 
    # Input:    A lattice of complex numbers L, a representation of a number
    #           field frep, and number representing a homomorphism h from that
    #           field into the complex numbers.
    # Output:   An algebraic equation of the corresponding elliptic curve.
    gp.set_precision(prec)
    g4 = gp.elleisnum(L, 4, flag = 1)
    g6 = gp.elleisnum(L, 6, flag = 1)
    eisvals_an = [ 12 * g4 , 216 * g6 ]
    eisvals_an = [ eisval.sage() for eisval in eisvals_an ]
    # TODO: Next line can be replaced by a smaller field; if we just keep track
    # of its hom then that gives descents.
    eisvals_alg = [ magma.AlgebraizeElementInField(magma(eisval), magma(frep),
        magma(h), epscomp = epscomp, epsLLL = epsLLL) for eisval in eisvals_an ]
    return eisvals_alg

# Part of TODO:
# This allows us to find the curves. To get points later, use
# gp.ellwp(w, 0.00001), which indeed has a quadratic pole.
