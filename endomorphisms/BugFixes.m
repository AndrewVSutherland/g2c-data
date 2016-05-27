// Bug fixes for automorphism groups and fixgroups of fields


function AutomorphismGroupFixed(L);
// Returns the automorphism group without duplicates.

Gf := Automorphisms(L);

function aut_comp(phi1, phi2);
// Returns the composition of two automorphisms
return map<L -> L | x:->phi2(phi1(x)) >;
end function;

function aut_eq(phi1, phi2);
//Assumes that L is simple
return phi1(L.1) eq phi2(L.1);
end function;

Gpres, m1 := GenericGroup(Gf : Mult:=aut_comp, Eq:=aut_eq, Id:=Gf[1]);
m2, Gp := CosetAction(Gpres, sub<Gpres | >); // get permutation action
Gphi := Inverse(m2)*m1;

return Gp, Gf, Gphi;

end function;


function FixedGroupFixed(L, K, Gp, Gphi);
// Assumes that K is simple (to skip DefiningPolynomial step)

HpList := [ g : g in Gp | L ! (Gphi(g))(K.1) eq L ! K.1 ];
if not #HpList * Degree(K) eq Degree(L) then
    error Error("Magma took an incorrect FixedGroup... again");
end if;

return sub<Gp | HpList >;

end function;
