# Real endomorphism algebra over Qbar: End(A_\Qbar) \otimes R, as determined by ST^0(A) = U(1),SU(2),U(1)xSU(2),SU(2)xSU(2),USp(4) -- see FKRS12
real_geom_end_alg = { \
"C_1" : "M_2(C)", \
"C_2" : "M_2(C)", \
"C_3" : "M_2(C)", \
"C_4" : "M_2(C)", \
"C_6" : "M_2(C)", \
"D_2" : "M_2(C)", \
"D_3" : "M_2(C)", \
"D_4" : "M_2(C)", \
"D_6" : "M_2(C)", \
"T" : "M_2(C)", \
"O" : "M_2(C)", \
"J(C_1)" : "M_2(C)", \
"J(C_2)" : "M_2(C)", \
"J(C_3)" : "M_2(C)", \
"J(C_4)" : "M_2(C)", \
"J(C_6)" : "M_2(C)", \
"J(D_2)" : "M_2(C)", \
"J(D_3)" : "M_2(C)", \
"J(D_4)" : "M_2(C)", \
"J(D_6)" : "M_2(C)", \
"J(T)" : "M_2(C)", \
"J(O)" : "M_2(C)", \
"C_{2,1}" : "M_2(C)", \
"C_{4,1}" : "M_2(C)", \
"C_{6,1}" : "M_2(C)", \
"D_{2,1}" : "M_2(C)", \
"D_{4,1}" : "M_2(C)", \
"D_{6,1}" : "M_2(C)", \
"D_{3,2}" : "M_2(C)", \
"D_{4,2}" : "M_2(C)", \
"D_{6,2}" : "M_2(C)", \
"O_1" : "M_2(C)", \
"E_1" : "M_2(R)", \
"E_2" : "M_2(R)", \
"E_3" : "M_2(R)", \
"E_4" : "M_2(R)", \
"E_6" : "M_2(R)", \
"J(E_1)" : "M_2(R)", \
"J(E_2)" : "M_2(R)", \
"J(E_3)" : "M_2(R)", \
"J(E_4)" : "M_2(R)", \
"J(E_6)" : "M_2(R)", \
"F" : "C x C", \
"F_a" : "C x C", \
"F_{ab}" : "C x C", \
"F_{ac}" : "C x C", \
"F_{a,b}" : "C x C", \
"G_{1,3}" : "C x R", \
"N(G_{1,3})" : "C x R", \
"G_{3,3}" : "R x R", \
"N(G_{3,3})" : "R x R", \
"USp(4)": "R" \
}

# Real endomorphism algebra over Q: End(A) \otimes R, taken from Table 8 in FKRS12
real_end_alg = { \
"C_1" : "M_2(C)", \
"C_2" : "C x C", \
"C_3" : "C x C", \
"C_4" : "C x C", \
"C_6" : "C x C", \
"D_2" : "C", \
"D_3" : "C", \
"D_4" : "C", \
"D_6" : "C", \
"T" : "C", \
"O" : "C", \
"J(C_1)" : "H", \
"J(C_2)" : "C", \
"J(C_3)" : "C", \
"J(C_4)" : "C", \
"J(C_6)" : "C", \
"J(D_2)" : "R", \
"J(D_3)" : "R", \
"J(D_4)" : "R", \
"J(D_6)" : "R", \
"J(T)" : "R", \
"J(O)" : "R", \
"C_{2,1}" : "M_2(R)", \
"C_{4,1}" : "C", \
"C_{6,1}" : "C", \
"D_{2,1}" : "R x R", \
"D_{4,1}" : "R", \
"D_{6,1}" : "R", \
"D_{3,2}" : "R x R", \
"D_{4,2}" : "R x R", \
"D_{6,2}" : "R x R", \
"O_1" : "R", \
"E_1" : "M_2(R)", \
"E_2" : "C", \
"E_3" : "C", \
"E_4" : "C", \
"E_6" : "C", \
"J(E_1)" : "R x R", \
"J(E_2)" : "R", \
"J(E_3)" : "R", \
"J(E_4)" : "R", \
"J(E_6)" : "R", \
"F" : "C x C", \
"F_a" : "C x R", \
"F_{ab}" : "R x R", \
"F_{ac}" : "R", \
"F_{a,b}" : "R x R", \
"G_{1,3}" : "C x R", \
"N(G_{1,3})" : "R x R", \
"G_{3,3}" : "R x R", \
"N(G_{3,3})" : "R", \
"USp(4)": "R" \
}

def qdisc(n):
    d = ZZ(n).squarefree_part()
    return d if d%4 == 0 or d%4 == 1 else 4*d

def make_disc_key(D):
    Dz = D.abs()
    D1 = int(Dz.log(10)) if Dz else 0
    return '%03d%s' % (D1, str(Dz))

# construct isogeny class labels -- returns a dictionary indexed by hash
def loadclasslabels(filename):
    conductors = {}
    for r in open(filename):
        cond = int(r.split(":")[1])
        conductors[cond] = {}
    for r in open(filename):
        s = r.split(":")
        cond=int(s[1])
        hash=int(s[2])
        if hash in conductors[cond]:
            continue
        n=len(conductors[cond])
        conductors[cond][hash] = str(cond)+"."+chr(ord('a')+n)
    labels = { h:v[h] for v in conductors.values() for h in v }
    return labels
    
def scanLdata(filename):
    fp = open(filename,"r")
    Ldataoffsets = {}
    while true:
        off = fp.tell()
        x = fp.readline()
        if not x: break
        hash = int(x.split(":")[0])
        Ldataoffsets[hash] = off
    return Ldataoffsets

# to deal with memory issues we process isogeny classes in half-open conductor intervals [a,b)
# this seems inefficient but isn't really, the time to rescan the file is negligible compares to the time to process all the Lfactors
def loadclasses(filename,ecproduct_file,ecquadratic_file,mfproduct_file,mfhilbert_file,mincond,maxcond,class_labels,Ldatafile,Ldataoffsets):
    ecproduct_dict = {}
    for r in open(ecproduct_file):
        s = r.split(":")
        ecproduct_dict[int(s[0])] = [eval(s[1]),eval(s[2])]
    ecquadratic_dict = {}
    for r in open(ecquadratic_file):
        s = r.split(":")
        ecquadratic_dict[int(s[0])] = eval(s[1])
    mfproduct_dict = {}
    for r in open(mfproduct_file):
        s = r.split(":")
        mfproduct_dict[int(s[0])] = eval(s[1])
    mfhilbert_dict = {}
    for r in open(mfhilbert_file):
        s = r.split(":")
        mfhilbert_dict[int(s[0])] = eval(s[1])
    R.<x>=PolynomialRing(QQ)
    hashes = {}
    classes = []
    Lfunctions = []
    for r in open(filename):
        s = r.split(":")
        cond=int(s[1])
        if cond < mincond or cond >= maxcond:
            continue
        hash=int(s[2])
        if hash in hashes:
            continue
        assert hash in class_labels
        hashes[hash] = 1
        badprimes = prime_factors(int(s[1]))
        badpolys = s[7].split(",")
        bad_lfactors = [[int(badprimes[i]),[int(c) for c in R(badpolys[i]).list()]] for i in range(len(badprimes))]
        good_lfactors = [[int(a[0]),int(a[1]),int(a[2])] for a in eval(s[13])]
        stgroup = s[8].strip()
        is_gl2_type = true if real_end_alg[stgroup] == "R x R" or real_end_alg[stgroup] == "C" else false
        rec = {"label":class_labels[hash],"cond":cond,"hash":hash,"root_number":int(s[6]),"st_group":stgroup,"real_end_alg":real_end_alg[stgroup],"real_geom_end_alg":real_geom_end_alg[stgroup],"is_gl2_type":is_gl2_type,"bad_lfactors":bad_lfactors,"good_lfactors":good_lfactors}
        rec["rat_end_alg"] = "Q" if real_end_alg[stgroup] == "R" else ""
        rec["rat_geom_end_alg"] = "Q" if real_geom_end_alg[stgroup] == "R" else ""
        rec["geom_end_field"] = "1.1.1.1" if real_geom_end_alg[stgroup] == "R" or stgroup == "G_{3,3}" or stgroup == "E_1" else ""
        if stgroup == "USp(4)":
            rec["is_simple"] = true
            rec["is_geom_simple"] = true
        if real_geom_end_alg[stgroup] == "M_2(C)":
            rec["is_geom_simple"] = false
        if real_geom_end_alg[stgroup] == "C x R":
            rec["is_geom_simple"] = false
        if stgroup == "F_{ac}":
            rec["is_simple"] = true
            rec["is_geom_simple"] = true
        if hash in ecproduct_dict:
            rec["ecproduct"] = ecproduct_dict[hash][0]
            rec["is_simple"] = false
            rec["is_geom_simple"] = false
            if stgroup == "G_{3,3}":
                rec["rat_end_alg"] = "Q x Q"
                rec["rat_geom_end_alg"] = "Q x Q"
            if stgroup == "E_1":
                rec["rat_end_alg"] = "M_2(Q)"
                rec["rat_geom_end_alg"] = "M_2(Q)"
            if stgroup == "N(G_{1,3})":
                rec["rat_end_alg"] = "Q x Q"
                n = min(ecproduct_dict[hash][1])
                assert n < 0
                n = qdisc(n)
                rec["rat_geom_end_alg"] = "Q x Qsqrt%d"%(n)
                rec["geom_end_field"] = "2.0.%d.1"%(-n)
        if not "is_simple" in rec:
            for a in good_lfactors:
                if a[0] > 100: break
                Lpoly = R([1,-a[1],a[2],-a[0]*a[1],a[0]*a[0]])
                if Lpoly.is_irreducible():
                    rec["is_simple"] = true
                    if rec["geom_end_field"] == "1.1.1.1":
                        rec["is_geom_simple"] = true
                    break
        if hash in ecquadratic_dict:
            rec["ecquadratic"] = ecquadratic_dict[hash]
        if hash in mfproduct_dict:
            rec["mfproduct"] = mfproduct_dict[hash]
        if hash in mfhilbert_dict:
            rec["mfhilbert"] = mfhilbert_dict[hash]
        if hash in Ldataoffsets:
            lfunc = {"hash":hash}
            fp = open(Ldatafile)
            fp.seek(Ldataoffsets[hash])
            lfuncdata = fp.readline()
            fp.close()
            t = lfuncdata.strip().split(":")
            assert t[0] == str(hash)
            assert int(t[1]) == cond
            assert int(t[2]) == int(s[6])
            lfactors = {}
            for x in bad_lfactors:
                lfactors[x[0]] = x[1]
            for a in good_lfactors:
                if a[0] >= 100:
                    break
                p = a[0]; a1 = a[1]; a2 = a[2]
                lfactors[p] = [int(1),int(-a1),int(a2),int(-p*a1),int(p*p)]
            lfunc["euler_factors"] = str([lfactors[p] for p in prime_range(0,100)])
            lfunc["bad_lfactors"] = bad_lfactors
            lfunc["special_values"] = t[3] # [[float(x[0]),float(x[1])] for x in pari(t[3])] 
            lfunc["zeros"] = t[4] # [float(x) for x in pari(t[4])]
            rank = len([x for x in eval(t[4]) if x == 0])
            lfunc["order_of_vanishing"] = int(rank)
            rec["analytic_rank"] = int(rank)
            lfunc["plot"] = t[5] # [[float(x[0]),float(x[1])] for x in pari(t[5])]
            lfunc["root_number"] = t[2]
            lfunc["analytic_normalization"] = "1/2"
            lfunc["degree"] = int(4)
            lfunc["gamma_factors"] = "[[],[0,0]]" # [[],[float(0),float(0)]] 
            lfunc["conductor"] = cond
            lfunc["algebraic"] = true
            lfunc["motivic_weight"] = int(1)
            lfunc["primitive"] = true if real_end_alg[stgroup] == "R" else false
            lfunc["instances"] = [ "/L/Genus2Curve/Q/%d/%s"%(cond,class_labels[hash].split(".")[-1]) ]
            lfunc["central_character"] = "%d.1"%(cond)
            Lfunctions.append(lfunc)
        classes.append(rec)
    return classes,Lfunctions

def loadcurves(filename,iso_classes):
    R.<x>=PolynomialRing(QQ)
    cdlabels = {}
    start = walltime()
    for r in open(filename):
        s = r.split(":")
        cond=int(s[1])
        hash=int(s[2])
        stgroup = s[8].strip()
        iso_class = iso_classes.find_one({"hash":hash})
        clabel = iso_class["label"]
        cdlabel = clabel+"."+s[0]
        if not cdlabel in cdlabels:
            cdlabels[cdlabel] = {}
        n = len(cdlabels[cdlabel])
        label=cdlabel+"."+chr(ord('1')+n)
        polys = s[3].split("[")[1].split("]")[0].split(",")
        f = R(polys[0])
        h = R(polys[1])
        g = 4*f + h^2
        assert g.discriminant() != 0 and (g.degree() == 5 or g.degree() == 6)
        num_rat_wpts = int(len(g.roots()) + 6-g.degree())
        disc_key = make_disc_key(ZZ(s[0]))
        disc_sign = int(s[4])
        inv = [str(i) for i in eval(s[5])]
        aut_grp = [int(n) for n in eval(s[9])]
        geom_aut_grp = [int(n) for n in eval(s[10])]
        torsion = [int(n) for n in eval(s[11])]
        torsion_order = int(prod(torsion))
        two_selmer_rank = int(eval(s[12]))
        min_eqn = [[int(a) for a in f.list()],[int(a) for a in h.list()]]
        rec = {"label":label,"cond":cond,"class":clabel,"disc_key":disc_key,"disc_sign":disc_sign,"min_eqn":min_eqn,"igusa_clebsch":inv,"aut_grp":aut_grp,"geom_aut_grp":geom_aut_grp,"torsion":torsion,"torsion_order":torsion_order,"num_rat_wpts":num_rat_wpts,"two_selmer_rank":two_selmer_rank}
        if iso_class["rat_end_alg"] == "Q":
            rec["end_ring"] = "Z"
        if iso_class["rat_geom_end_alg"] == "Q":
            rec["geom_end_ring"] = "Z"
        # duplicate certain attributes from isogeny class for search purposes
        rec["st_group"] = stgroup
        rec["real_geom_end_alg"] = real_geom_end_alg[stgroup]
        rec["is_gl2_type"] = iso_class["is_gl2_type"]
        if "is_simple" in iso_class.keys():
            rec["is_simple"] = iso_class["is_simple"]
        if "is_geom_simple" in iso_class.keys():
            rec["is_geom_simple"] = iso_class["is_geom_simple"]
        cdlabels[cdlabel][label] = rec
    curves = [ rec for v in cdlabels.values() for rec in v.values() ]
    return curves

    
# convert products of elliptic curves specified by cremona labels to a pair of lmfdb isogeny class labels
def convert_ecproduct_labels(infile,outfile,lmfdb):
    outfp = open(outfile,"w")
    for r in open(infile):
        s = r.split(":")
        labels = s[2].strip().split("*")
        e1 = lmfdb.elliptic_curves.curves.find_one({"label":labels[0]})
        label1 = e1["lmfdb_iso"]
        cm1 = e1["cm"]
        e2 = lmfdb.elliptic_curves.curves.find_one({"label":labels[1]})
        label2 = e2["lmfdb_iso"]
        cm2 = e2["cm"]
        outfp.write("%s:['%s','%s']:[%d,%d]\n"%(s[1],label1,label2,cm1,cm2))
    outfp.close()
