load('Initialize.sage')

# Ambient ring:
R.<x> = PolynomialRing(QQ)
# Perhaps make these ratios / powers? Anyway, voodoo for now.
prec = 300
epscomp = 10^(-prec + 30)
epsLLL = 5^(-prec + 7)
epsinv = 2^(-4)
Bound = 48

# Curve input: specify g and h in its equation y^2 + h y = g.

# Some interesting curves:
# Many automorphisms, and descent to a smaller field down from degree 12:
# (a personal favorite)
g = R(1)
h = x^3 + 1
# Another degree 12 field:
#g = -2*x^6 + 2*x^3 - 1
#h = x^3 + 1
# Endomorphisms and splittings defined over QQ:
#g = x^4 + x^2
#h = x^3 + 1
# Case where [RR, CC] occurs:
g = x^5 + x^4 + 2*x^3 + x^2 + x
h = x^2 + x
# The big degree 48 case:
#g = x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1
#h = R(0)
# Purported restriction of scalars:
g = x^6 + 2*x^5 + 2*x^4 - 9*x^3 -12*x^2 - 3*x + 26
h = x

# Directly from database:
# Typical (USp(4)) case:
#[g,h] = [[-1,-2,-1,1,2,2,1],[1]]
# Infinite loop:
#[g,h] = [[2,0,-64,0,-16,0,-1],[0,1]]
# Hard to jiggle:
#[g,h] = [[-3,-14,3,35,-53,42,-14],[0,1,1]]

# Coercion if necessary
g = R(g)
h = R(h)

subst = x
#subst = x + 1
#subst = -2/3*x + 1 
#subst = (-1/3*x + 1/3)/(-3/2*x + 2)
den = subst.denominator()
g = R(den^6 * g(subst))
h = R(den^3 * h(subst))

LGEB = magma.GeometricEndomorphismBasis(h, g, prec = prec, epscomp = epscomp,
        epsLLL = epsLLL, epsinv = epsinv, nvals = 3)
As = LGEB[0]
Rs = LGEB[1]
P = LGEB[2]
print As
print Rs
print P

AsPol = [ magma.PolynomializeMatrix(A, epscomp = epscomp, epsLLL = epsLLL) for A
        in As ]
frep = Common_Splitting_Field(AsPol, Bound = Bound)
print frep

LAM = magma.AlgebraizeMatricesInField(As, AsPol, frep, epscomp = epscomp, nvals
        = 2)
AsAlg = LAM[0]
fhom = LAM[1]
print AsAlg

LEL = magma.EndomorphismLatticeG2(AsAlg, As, Rs, Geometric = false, AddTensor =
        true, AddRing = true, AddSatoTate = true, AddSplitting = true, nvals =
        3)
EDs = LEL[0]
fsubgen = LEL[1]
idems = LEL[2]
print idems

EDs_sage = []
for ED in EDs:
    ED1 = Canonize_Subfield(ED[1].sage(), frep)
    ED2 = [ [ Canonize_Field(factorsQQ[1].sage())[0], factorsQQ[2].sage() ] for
            factorsQQ in ED[2] ]
    ED3 = [ repr(factorRR) for factorRR in ED[3] ]
    ED4 = ED[4].sage()
    ED5 = repr(ED[5])
    EDs_sage.append([ED1, ED2, ED3, ED4, ED5])
print EDs_sage

# Optimize:
SplFoD = Canonize_Subfield(fsubgen.sage(), frep)
fsubrep_opt = SplFoD[0]
fsubgen_opt = SplFoD[1]
fsubhom_opt = magma.InducedEmbedding(fsubgen_opt, fhom, frep).sage()
print SplFoD

Lats = magma.LatticesFromIdempotents(idems[2], P, epscomp = epscomp, epsLLL =
        epsLLL, epsinv = epsinv)
ECs = [ Elliptic_Curve_From_Lattice(Lat.sage(), fsubrep_opt,
    fsubhom_opt, prec = prec, epscomp = epscomp, epsLLL = epsLLL) for
    Lat in Lats ]
print ECs

# Optional file printing:
#outputfile = 'Examples_Output.txt'
#with open(outputfile, 'w') as outputstream:
#    outputstream.write(str(frep) + '\n')
#    outputstream.write(str(EDs_sage) + '\n')
#    outputstream.write(str(SplFoD) + '\n')
#    outputstream.write(str(ECs) + '\n')
