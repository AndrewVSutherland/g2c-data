import ast

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

# Map automorphism GAP ids to string labels (eventually these should be labels of automorphism groups in the LMFDB), storing pairs of ints in mongo-DB is unhelpful
aut_grp_ids = {
"[2,1]" : "2G1",\
"[4,1]" : "4G1",\
"[4,2]" : "4G2",\
"[6,2]" : "6G2",\
"[8,3]" : "8G3",\
"[10,2]" : "10G2",\
"[12,4]" : "12G4",\
"[24,8]" : "24G8",\
"[48,29]" : "48G29",\
}

def igusa_clebsch_to_igusa(I):
    # Conversion from Igusa-Clebsch to Igusa
    J2 = I[0]//8
    J4 = (4*J2**2 - I[1])//96
    J6 = (8*J2**3 - 160*J2*J4 - I[2])//576
    J8 = (J2*J6 - J4**2)//4
    J10 = I[3]//4096
    return [J2,J4,J6,J8,J10]

def igusa_to_g2(J):
    # Conversion from Igusa to G2
    if J[0] != 0:
        return [J[0]**5/J[4], J[0]**3*J[1]/J[4], J[0]**2*J[2]/J[4]]
    elif J[1] != 0:
        return [0, J[1]**5/J[4]**2, J[1]*J[2]/J[4]]
    else:
        return [0,0,J[2]**5/J[4]**3]

def qdisc(n):
    d = ZZ(n).squarefree_part()
    return d if d%4 == 0 or d%4 == 1 else 4*d

def make_disc_key(D):
    Dz = D.abs()
    D1 = int(Dz.log(10)) if Dz else 0
    return '%03d%s' % (D1, str(Dz))
    
def list2string(li):
    li2 = [str(x) for x in li]
    return ','.join(li2)

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
def loadclasses(filename,ecproduct_file,ecquadratic_file,mfproduct_file,mfhilbert_file,mincond,maxcond,class_labels,Ldatafile,Ldataoffsets,load_plot=1):
    ecproduct_dict = {}
    for r in open(ecproduct_file):
        s = r.split(":")
        ecproduct_dict[int(s[0])] = [ast.literal_eval(s[1]),ast.literal_eval(s[2])]
    ecquadratic_dict = {}
    for r in open(ecquadratic_file):
        s = r.split(":")
        ecquadratic_dict[int(s[0])] = ast.literal_eval(s[1])
    mfproduct_dict = {}
    for r in open(mfproduct_file):
        s = r.split(":")
        mfproduct_dict[int(s[0])] = ast.literal_eval(s[1])
    mfhilbert_dict = {}
    for r in open(mfhilbert_file):
        s = r.split(":")
        mfhilbert_dict[int(s[0])] = ast.literal_eval(s[1])
    R.<x>=PolynomialRing(QQ)
    hashes = {}
    classes = []
    Lfunctions = []
    Linstances = []
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
        good_lfactors = [[int(a[0]),int(a[1]),int(a[2])] for a in ast.literal_eval(s[16])]
        stgroup = s[8].strip()
        is_gl2_type = true if real_end_alg[stgroup] == "R x R" or real_end_alg[stgroup] == "C" else false
        rec = {"label":class_labels[hash],"cond":cond,"hash":hash,"Lhash":str(hash),"root_number":int(s[6]),"st_group":stgroup,"real_geom_end_alg":real_geom_end_alg[stgroup],"is_gl2_type":is_gl2_type,"bad_lfactors":bad_lfactors}
        if hash in ecquadratic_dict:
            rec["ecquadratic"] = ecquadratic_dict[hash]
        if hash in mfproduct_dict:
            rec["mfproduct"] = mfproduct_dict[hash]
        if hash in mfhilbert_dict:
            rec["mfhilbert"] = mfhilbert_dict[hash]
        if hash in Ldataoffsets:
            lfunc = {"Lhash":str(hash)}
            lfunc["load_key"] = "G2Q"
            fp = open(Ldatafile)
            fp.seek(Ldataoffsets[hash])
            lfuncdata = fp.readline()
            fp.close()
            t = lfuncdata.strip().split(":")
            # L-function data format
            # hash (string) : cond (int) : root_number (int) : analytic_rank (int) : dirichlet_coeffs (list of 9 ints a2,...,a10) : leading_term (float): positive_zeros (list of floats) : plot_delta : plot_values (list of floats) 
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
                lfactors[p] = [int(1),int(a1),int(a2),int(p*a1),int(p*p)]
            lfunc["euler_factors"] = [[int(a) for a in lfactors[p]] for p in prime_range(0,100)]
            lfunc["bad_lfactors"] = bad_lfactors
            lfunc["root_number"] = int(t[2])
            lfunc["order_of_vanishing"] = int(t[3])
            rec["analytic_rank"] = int(t[3])
            As = [int(x) for x in pari(t[4])]
            lfunc["A2"] = int(As[0])
            lfunc["A3"] = int(As[1])
            lfunc["A4"] = int(As[2])
            lfunc["A5"] = int(As[3])
            lfunc["A6"] = int(As[4])
            lfunc["A7"] = int(As[5])
            lfunc["A8"] = int(As[6])
            lfunc["A9"] = int(As[7])
            lfunc["A10"] = int(As[8])
            lfunc["a2"] = [float(As[0])/float(p),float(0)]
            lfunc["a3"] = [float(As[1])/float(p),float(0)]
            lfunc["a4"] = [float(As[2])/float(p),float(0)]
            lfunc["a5"] = [float(As[3])/float(p),float(0)]
            lfunc["a6"] = [float(As[4])/float(p),float(0)]
            lfunc["a7"] = [float(As[5])/float(p),float(0)]
            lfunc["a8"] = [float(As[6])/float(p),float(0)]
            lfunc["a9"] = [float(As[7])/float(p),float(0)]
            lfunc["a10"] = [float(As[8])/float(p),float(0)]
            lfunc["leading_term"] = float(t[5])
            lfunc["positive_zeros"] = [float(x) for x in pari(t[6])]
            lfunc["z1"] = lfunc["positive_zeros"][0]
            lfunc["z2"] = lfunc["positive_zeros"][1]
            lfunc["z3"] = lfunc["positive_zeros"][2]
            if load_plot:
                lfunc["plot_delta"] = float(t[7])
                lfunc["plot_values"] = [float(x) for x in pari(t[8])]
            lfunc["analytic_normalization"] = "1/2"
            lfunc["degree"] = int(4)
            lfunc["gamma_factors"] = [[],[int(0),int(0)]]
            lfunc["conductor"] = cond
            lfunc["algebraic"] = true
            lfunc["motivic_weight"] = int(1)
            lfunc["primitive"] = true if real_end_alg[stgroup] == "R" else false
            lfunc["origin"] = "Genus2Curve/Q/%d/%s"%(cond,class_labels[hash].split(".")[-1])
            lfunc["central_character"] = "%d.1"%(cond)
            lfunc["self_dual"] = true
            lfunc["coefficient_field"] = "1.1.1.1"
            lfunc["st_group"] = stgroup
            Lfunctions.append(lfunc)
            linst={"Lhash":str(hash)}
            linst["url"]=lfunc["origin"]
            linst["type"]="G2Q"
            Linstances.append(linst)
        classes.append(rec)
    return classes,Lfunctions,Linstances

def loadcurves(filename,iso_classes):
    R.<x>=PolynomialRing(QQ)
    cdlabels = {}
    for r in open(filename):
        s = r.split(":")
        cond=int(s[1])
        hash=int(s[2])
        stgroup = s[8].strip()
        iso_class = iso_classes.find_one({"Lhash":str(hash)})
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
        abs_disc = int(ZZ(s[0]))
        disc_key = make_disc_key(ZZ(s[0]))
        disc_sign = int(s[4])
        igusa_clebsch = [str(i) for i in ast.literal_eval(s[5])]
        igusa = [str(i) for i in igusa_clebsch_to_igusa(ast.literal_eval(s[5]))]
        g2inv = [str(i) for i in igusa_to_g2 (igusa_clebsch_to_igusa(ast.literal_eval(s[5])))]
        aut_grp = [int(n) for n in ast.literal_eval(s[9])]
        geom_aut_grp = [int(n) for n in ast.literal_eval(s[10])]
        aut_grp_id = aut_grp_ids["".join(s[9].split())] # use join of null split to remove whitespace
        geom_aut_grp_id = aut_grp_ids["".join(s[10].split())]
        torsion = [int(n) for n in ast.literal_eval(s[11])]
        torsion_order = int(prod(torsion))
        two_selmer_rank = int(ast.literal_eval(s[12]))
        has_square_sha = int(ast.literal_eval(s[13]))
        assert has_square_sha in [0,1]
        has_square_sha = true if has_square_sha == 1 else false
        locally_solvable = int(ast.literal_eval(s[14]))
        assert locally_solvable in [0,1]
        locally_solvable = true if locally_solvable == 1 else false
        globally_solvable = int(ast.literal_eval(s[15]))
        assert globally_solvable in [-1,0,1]
        min_eqn = [[int(a) for a in f.list()],[int(a) for a in h.list()]]
        analytic_rank = iso_class["analytic_rank"]
        rec = {"label":label,"cond":cond,"class":clabel,"abs_disc":abs_disc,"disc_key":disc_key,"disc_sign":disc_sign,"min_eqn":min_eqn,"igusa_clebsch":igusa_clebsch,"igusa":igusa,"g2inv":g2inv,
                    "aut_grp":aut_grp,"geom_aut_grp":geom_aut_grp, # remove this line once these are no longer needed
                    "aut_grp_id":aut_grp_id,"geom_aut_grp_id":geom_aut_grp_id,
                   "torsion":torsion,"torsion_order":torsion_order,"num_rat_wpts":num_rat_wpts,"analytic_rank":analytic_rank,"two_selmer_rank":two_selmer_rank,"has_square_sha":has_square_sha,"locally_solvable":locally_solvable,"globally_solvable":globally_solvable}
        # duplicate certain attributes from isogeny class for search purposes
        rec["st_group"] = stgroup
        rec["real_geom_end_alg"] = real_geom_end_alg[stgroup]
        rec["is_gl2_type"] = iso_class["is_gl2_type"]
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

def lookup_numberfield(min_eqn,lmfdb):
    def from_polynomial(cls, pol):
        pol = PolynomialRing(QQ, 'x')(str(pol))
        pol *= pol.denominator()
        R = pol.parent()
        pol = R(pari(pol).polredabs())
        return cls.from_coeffs([int(c) for c in pol.coefficients(sparse=False)])

