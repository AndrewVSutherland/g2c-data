password = open("password.txt", "r").readlines()[0].strip()
from pymongo import MongoClient
update_Lfunctions = True
maxcond = 1000000
delta = 10000
curvefile = "gce_1000000_lmfdb.txt"
ecproduct = "gce_1000000_ecproduct_lmfdb.txt"
ecquadratic = "gce_1000000_ecquadratic_lmfdb.txt"
mfhilbert = "gce_1000000_mfhilbert_lmfdb.txt"
mfproduct = "gce_1000000_mfproduct_lmfdb.txt"
Ldatafile = "gce_1000000_ldata.txt"
conn = MongoClient(host='localhost',port=int(37010))  # warwick  # warwick
conn.admin.authenticate("lmfdb","lmfdb") # get read access to everything (so we can lookup labels)
print "Connected to LMFDB"
db = conn.genus2_curves
Ldb = conn.Lfunctions
db.authenticate('editor', password) 
Ldb.authenticate('editor', password)
print "Removing L-function data"
if update_Lfunctions:
    Ldb.Lfunctions.remove({"load_key":"G2Q"})
    Ldb.instances.remove({"type":"G2Q"})
begin=walltime()
load("loadcurves.sage")
print "Dropping new curves files (if present)"
db.drop_collection("new_curves")
db.drop_collection("new_isogeny_classes")
print "determining isogeny class labels in %s"%(curvefile)
start=walltime()
class_labels=loadclasslabels(curvefile)
print "constructed %d isogeny class labels in %.3f secs"%(len(class_labels),walltime()-start)
print "loading L-function data from %s"%(Ldatafile)
start=walltime()
Ldataoffsets=scanLdata(Ldatafile)
print "found data for %d L-functions in %.3f secs"%(len(Ldataoffsets),walltime()-start)
for mincond in range(0,maxcond+1,delta):
    print "loading isogeny classes from %s"%(curvefile)
    start=walltime()
    isogeny_classes, Lfunctions,Linstances = loadclasses(curvefile,ecproduct,ecquadratic,mfproduct,mfhilbert,mincond,mincond+delta,class_labels,Ldatafile,Ldataoffsets,load_plot=1) # change load_plot to 1 to get plots
    print "loaded %d isogeny classes and %d L-functions in conductor range [%d,%d) in %.3f secs"%(len(isogeny_classes),len(Lfunctions),mincond,mincond+delta,walltime()-start)
    if len(isogeny_classes) == 0:
        continue
    print "inserting isogeny class data into LMFDB"
    start = walltime()
    if isogeny_classes:
        db.new_isogeny_classes.insert(isogeny_classes)
    print "inserted %d isogeny classes in conductor range [%d,%d) into LMFDB in %.3f secs"%(len(isogeny_classes),mincond,mincond+delta,walltime()-start)
    start = walltime()
    if update_Lfunctions and Lfunctions:
        Ldb.Lfunctions.insert(Lfunctions)
        Ldb.instances.insert(Linstances)
        print "inserted %d L-functions in conductor range [%d,%d) into LMFDB in %.3f secs"%(len(Lfunctions),mincond,mincond+delta,walltime()-start)
    print "inserted a total of %d isogeny classes < %d into LMFDB in %.3f secs"%(db.new_isogeny_classes.count(),mincond+delta,walltime()-begin)
db.new_isogeny_classes.ensure_index("Lhash")
db.new_isogeny_classes.ensure_index("hash")
db.new_isogeny_classes.ensure_index("label")
db.new_isogeny_classes.ensure_index("cond")
db.new_isogeny_classes.ensure_index("analytic_rank")
db.new_isogeny_classes.ensure_index("st_group")
print "loading curves from %s"%(curvefile)
start=walltime()
curves=loadcurves(curvefile,db.new_isogeny_classes)
print "loaded %d curves in %.3f secs"%(len(curves),walltime()-start)
print "inserting curve data into LMFDB"
start=walltime()
if curves:
    db.new_curves.insert(curves)
print "inserted %d curves in %.3f secs"%(db.new_curves.count(),walltime()-start)
db.new_curves.ensure_index("label")
db.new_curves.ensure_index("cond")
db.new_curves.ensure_index("class")
db.new_curves.ensure_index("abs_disc")
db.new_curves.ensure_index("disc_key")
db.new_curves.ensure_index("is_gl2_type")
db.new_curves.ensure_index("st_group")
db.new_curves.ensure_index("aut_grp")
db.new_curves.ensure_index("aut_grp_id")
db.new_curves.ensure_index("geom_aut_grp")
db.new_curves.ensure_index("geom_aut_grp_id")
db.new_curves.ensure_index("igusa_clebsch")
db.new_curves.ensure_index("g2inv")
db.new_curves.ensure_index("analytic_rank")
db.new_curves.ensure_index("real_geom_end_alg")
db.new_curves.ensure_index("torsion")
db.new_curves.ensure_index("torsion_order")
db.new_curves.ensure_index("locally_solvable")
db.new_curves.ensure_index("has_square_sha")
db.new_curves.ensure_index([("cond",int(1)),("class",int(1)),("disc_key",int(1)),("label",int(1))])
print "Total time %.3f secs"%(walltime()-begin)
print "Checking data..."

def check_labels(oldcurves,newcurves):
    i = 0
    for c1 in oldcurves.find():
        c2 = newcurves.find_one({"label":c1["label"]})
        for k in c1.keys():
            if k != "_id" and k != "min_eqn":
                if c1[k] != c2[k]:
                    print "old:",c1
                    print "new:", c2
                    print "key mismatch",k,c1[k],c2[k]
                    assert c1[k]==c2[k]
        i = i+1
    print "verified", i, "records"

check_labels(db.curves,db.new_curves)

db.drop_collection("old_curves")
db.drop_collection("old_isogeny_classes")
db.curves.rename("old_curves")
db.isogeny_classes.rename("old_isogeny_classes")
db.new_curves.rename("curves")
db.new_isogeny_classes.rename("isogeny_classes")
print "New data copied in"


