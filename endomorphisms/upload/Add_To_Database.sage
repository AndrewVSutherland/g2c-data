# Main upload file for old format

#inputfile = 'database_output_final.txt'
# For testing and debugging:
inputfile = 'database_output_test.txt'

with open(inputfile) as inputstream:
    num_lines = sum(1 for line in inputstream)
print num_lines

# Ambient needed in what follows:
R.<x> = PolynomialRing(QQ)
counter_done = 0
# For resuming:
#counter_done = 56300
while counter_done < num_lines:
    import pymongo
    lmfdb = pymongo.MongoClient(host = 'localhost', port = int(37010))
    lmfdb.numberfields.authenticate('editor', '282a29103a17fbad')
    lmfdb.elliptic_curves.authenticate('editor', '282a29103a17fbad')
    lmfdb.genus2_curves.authenticate('editor', '282a29103a17fbad')
    nfdb = lmfdb.numberfields.fields
    ecdb_QQ = lmfdb.elliptic_curves.curves
    ecdb_nf = lmfdb.elliptic_curves.nfcurves
    endodb = lmfdb.genus2_curves.endomorphisms
    with open(inputfile) as inputstream:
        try:
            counter_current = 0
            for line in inputstream:
                counter_current += 1
                if counter_current > counter_done:
                    counter_done += 1
                    print counter_done
                    linestrip = line.rstrip()
                    linesplit = linestrip.split(':')

                    # Starting new dictionary for the curve:
                    endodata = {'label': linesplit[0]}

                    # Input endomorphism lattice and prettify:
                    EDs = eval(linesplit[3])
                    EDs_se = sage_eval(linesplit[3])
                    N = len(EDs)
                    for i in range(N):
                        # Use strings in subfield description and add labels:
                        EDs[i][0][1] = [ repr(c) for c in EDs_se[i][0][1] ]
                        coeffs = str(EDs[i][0][0]).replace(' ', '').replace('[','').replace(']', '')
                        nf = nfdb.find_one({'coeffs': coeffs})
                        if nf:
                            EDs[i][0].insert(0, nf['label'])
                        else:
                            EDs[i][0].insert(0, '')
                        for factorQQ in EDs[i][1]:
                            # Add labels for subfields:
                            coeffs = str(factorQQ[0]).replace(' ', '').replace('[','').replace(']', '')
                            nf = nfdb.find_one({'coeffs': coeffs})
                            if nf:
                                factorQQ.insert(0, nf['label'])
                            else:
                                factorQQ.insert(0, '')
                    EDs = sorted(EDs, key = lambda t : len(t[0][1]))
                    endodata['lattice'] = EDs

                    # Field of definition:
                    endodata['fod_label'] = EDs[N - 1][0][0]
                    endodata['fod_coeffs'] = EDs[N - 1][0][1]

                    # Information over QQ and QQbar:
                    endodata['factorsQQ_base'] = EDs[0][1]
                    endodata['factorsQQ_geom'] = EDs[N - 1][1]
                    endodata['factorsRR_base'] = EDs[0][2]
                    endodata['factorsRR_geom'] = EDs[N - 1][2]
                    endodata['ring_base'] = EDs[0][3]
                    endodata['ring_geom'] = EDs[N - 1][3]
                    # Following line is superfluous because of check:
                    endodata['st_group_base'] = EDs[0][4]
                    # This is also superfluous, but I keep it in for symmetry reasons:
                    endodata['st_group_geom'] = EDs[N - 1][4]

                    # Simple if there is one factor of the algebra, not a quaternion
                    # algebra:
                    endodata['is_simple_base'] = (len(EDs[0][1]) == 1 and EDs[0][1][0][2] != 1)
                    endodata['is_simple_geom'] = (len(EDs[N - 1][1]) == 1 and EDs[N - 1][1][0][2] != 1)

                    # Splitting information:
                    subfield_rep = sage_eval(linesplit[4])
                    endodata['spl_fod_coeffs'] = [ int(c) for c in subfield_rep[0] ]
                    endodata['spl_fod_gen'] = [ repr(c) for c in subfield_rep[1] ]
                    factors_se = sage_eval(linesplit[5])
                    endodata['spl_facs_coeffs'] = [ [ [ repr(c) for c in coeff ] for coeff in factor ] for factor in factors_se ]
                    # Try to refine:
                    endodata['spl_facs_labels'] = []
                    endodata['spl_facs_condnorms'] = []
                    coeffs = str(subfield_rep[0]).replace(' ', '').replace('[','').replace(']', '')
                    nf = nfdb.find_one({'coeffs': coeffs})
                    field_label = ''
                    if nf:
                        field_label = nf['label']
                    endodata['spl_fod_label'] = field_label
                    # If the field is available, then we can polish the
                    # factors:
                    K.<r> = NumberField(R(subfield_rep[0]))
                    for factor in factors_se:
                        Ecoeffs0 = [ K(c) for c in factor ]
                        E0 = EllipticCurve([ -Ecoeffs0[0] / 48, -Ecoeffs0[1] / 864 ])
                        # Now look in the database for curves with the same
                        # j-invariant and conductor (norm) and check if they
                        # are isomorphic:
                        if K.degree() == 1:
                            conductor = int(E0.conductor())
                            endodata['spl_facs_condnorms'].append(int(conductor))
                            jinv = repr(E0.j_invariant())
                            ecs = ecdb_QQ.find({'conductor': conductor, 'jinv': jinv})
                            for ec in ecs:
                                E = EllipticCurve([ QQ(a.encode('ascii')) for a in ec['ainvs'] ])
                                if E0.is_isomorphic(E):
                                    endodata['spl_facs_labels'].append(ec['lmfdb_label'])
                        else:
                            conductor_norm = int(E0.conductor().norm())
                            endodata['spl_facs_condnorms'].append(int(conductor_norm))
                            jinv = [ repr(c) for c in E0.j_invariant().list() ]
                            ecs = ecdb_nf.find({'field_label': field_label, 'conductor_norm': conductor_norm, 'jinv': jinv})
                            for ec in ecs:
                                E = EllipticCurve([ K([ sage_eval(c) for c in a ]) for a in ec['ainvs'] ])
                                if E0.is_isomorphic(E):
                                    endodata['spl_facs_labels'].append(ec['label'])

                    # Finally!
                    if endodb.find_one({'label': endodata['label']}):
                        endodb.find_one_and_replace({'label': endodata['label']}, endodata)
                    else:
                        endodb.insert(endodata)
                    #print "Addition done"

        except:
            # Back a bit further, just to be sure
            print "Exception encountered"
            counter_done -= 1

# Check at end:
endodb = lmfdb.genus2_curves.endomorphisms
print endodb.find_one({'label': '169.a.169.1'})
print endodb.count()
