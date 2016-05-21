# Generates endomorphism data for gce format

import os, shutil

load('Initialize.sage')

# Ambient ring:
R.<x> = PolynomialRing(QQ)
# Parameters; perhaps make these ratios / powers? Anyway, voodoo for now.
prec = 300
epscomp = 10^(-prec + 30)
epsLLL = 5^(-prec + 7)
epsinv = 2^(-4)
Bound = 48

# Specify first input and output:
#inputfile = 'LMFDB/gce_lmfdb_input.txt'
#outputfile = 'LMFDB/gce_lmfdb_output.txt'
# For debugging:
inputfile = 'LMFDB/gce_delta.txt'
outputfile = 'LMFDB/gce_delta_output.txt'

# Safeguards against Van Wamelen's infinite loops and other errors:
# First substitutions (do not touch or an infinite loop will result):
subst_list = [ x + 1, x ]
# Bound on height of random substitutions after these have been tried:
B = 3
# The maximum number of retries:
maxrun = 2^5
# The counter for the current run and the line:
run = 0
counter = 0
# The curve that is death:
hell = [ 56306 ]

# Automatic determination of line length of input:
with open(inputfile) as inputstream:
    n = min([ len(line.split(':')) for line in inputstream ])

stop = False
exhaust = False
while not stop:
    run += 1
    counter = 0
    exhaust = len(subst_list) == 0
    if not exhaust:
        subst = subst_list.pop()
    else:
        while True:
            num = R([QQ.random_element(B) for i in range (2)])
            den = R([QQ.random_element(B) for i in range (2)])
            if num != 0 and den != 0:
                break
            # We do not have to guard against the constant case, since that
            # bugs in Magma already and will get caught below.
        subst = num/den
    stop = True
    print "Run:", run
    print "Substitution:", subst
    with open(inputfile) as inputstream:
        with open(outputfile, 'w') as outputstream:
            for line in inputstream:
                counter += 1
                linestrip = line.rstrip()
                linesplit = linestrip.split(':')
                linestart = ":".join(linesplit[:n])
                # We have to calculate if there is no information on the line
                # already:
                if len(linesplit) == n:
                    # In the USp(4) case we know everything:
                    if linesplit[8] == "USp(4)":
                        outputstream.write(linestart + ':' +
                                # Modify this if the data display changes
                                "[[[[0,1],[0]],[[[0,1],-1]],['RR'],[1,-1],'USp(4)']]:[[0,1],[0]]:[]"
                                + '\n')
                    # Avoiding a nasty infinite loop:
                    elif (subst == x) and (counter in hell):
                        outputstream.write(line)
                    else:
                        print counter
                        pol_list = eval(linesplit[3])
                        g = R(pol_list[0])
                        h = R(pol_list[1])
                        den = subst.denominator()
                        g = R(den^6 * g(subst))
                        h = R(den^3 * h(subst))
                        try:
                            #See Examples.sage for comments on this code
                            LGEB = magma.GeometricEndomorphismBasis(h, g, prec
                                    = prec, epscomp = epscomp, epsLLL = epsLLL,
                                    epsinv = epsinv, nvals = 3)
                            As = LGEB[0]
                            Rs = LGEB[1]
                            P = LGEB[2]
                            AsPol = [ magma.PolynomializeMatrix(A, epscomp =
                                epscomp, epsLLL = epsLLL) for A in As ]
                            frep = Common_Splitting_Field(AsPol, Bound = Bound)
                            LAM = magma.AlgebraizeMatricesInField(As, AsPol,
                                    frep, epscomp = epscomp, nvals = 2)
                            AsAlg = LAM[0]
                            fhom = LAM[1]
                            LEL = magma.EndomorphismLatticeG2(AsAlg, As, Rs,
                                    Geometric = false, AddTensor = true,
                                    AddRing = true, AddSatoTate = true,
                                    AddSplitting = true, nvals = 3)
                            EDs = LEL[0]
                            fsubgen = LEL[1]
                            idems = LEL[2]
                            EDs_sage = []
                            for ED in EDs:
                                ED1 = Canonize_Subfield(ED[1].sage(), frep)
                                ED2 = [ [ Canonize_Field(factorsQQ[1].sage())[0],
                                    factorsQQ[2].sage() ] for factorsQQ in
                                    ED[2] ]
                                ED3 = [ repr(factorRR) for factorRR in ED[3] ]
                                ED4 = ED[4].sage()
                                ED5 = repr(ED[5])
                                EDs_sage.append([ED1, ED2, ED3, ED4, ED5])
                            SplFoD = Canonize_Subfield(fsubgen.sage(), frep)
                            fsubrep_opt = SplFoD[0]
                            fsubgen_opt = SplFoD[1]
                            fsubhom_opt = magma.InducedEmbedding(fsubgen_opt,
                                    fhom, frep).sage()
                            Lats = magma.LatticesFromIdempotents(idems[2], P,
                                    epscomp = epscomp, epsLLL = epsLLL, epsinv
                                    = epsinv)
                            ECs = [ Elliptic_Curve_From_Lattice(Lat.sage(),
                                fsubrep_opt, fsubhom_opt, prec = prec, epscomp
                                = epscomp, epsLLL = epsLLL) for Lat in Lats ]
                            if EDs_sage[len(EDs_sage) - 1][4] == linesplit[8]:
                                outputstream.write(linestart
                                        + ':' + repr(EDs_sage).replace('\n', '').replace(' ', '')
                                        + ':' + repr(SplFoD).replace('\n', '').replace(' ', '')
                                        + ':' + repr(ECs).replace('\n', '').replace(' ', '')
                                        + '\n')
                            else:
                                # Case where we do not add information because
                                # the Sato-Tate groups do not match:
                                outputstream.write(line)
                        except:
                            # Case where we do not add information because of
                            # an error:
                            outputstream.write(line)
                            stop = False
                else:
                    # Case where we skip because we had already calculated this
                    # line:
                    outputstream.write(line)
    inputfile = 'LMFDB/intermediate.txt'
    shutil.copyfile(outputfile, inputfile)
    if run >= maxrun:
        stop = True
os.remove(inputfile)
