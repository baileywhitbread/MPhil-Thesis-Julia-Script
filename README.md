# MPhil Thesis Julia Script

This is the Julia script appearing in my MPhil thesis. 

It runs on Julia version 1.10.2 and the current version of Chevie as of 13/06/24.

This script computes representation-theoretic data needed to compute counting polynomials of character varieties.

## 14/06/24
I found a problem with the <tt>`orbit`</tt> function. 

If <tt>`G=rootdatum(:gl,4)`</tt> and <tt>`pseudo_levi_orbit_reps = reflection_subgroup.(Ref(G),sscentralizer_reps(G))`</tt> then <tt>`pseudo_levi_orbits = orbits(G,pseudo_levi_orbit_reps)`</tt> is incorrect. 

There are five copies of the $A_2$ root system when there should only be four copies. 

Replacing <tt>`G=rootdatum(:gl,4)`</tt> with <tt>`G=coxgroup(:A,3)`</tt> does not resolve this issue.