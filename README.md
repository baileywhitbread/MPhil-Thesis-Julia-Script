# MPhil Thesis Julia Script

This is the Julia script appearing in my MPhil thesis. 

It runs on Julia version 1.10.2 and the current version of Chevie as of 13/06/24.

This script computes representation-theoretic data needed to compute counting polynomials of character varieties.

## 14/06/24
I found a problem with the `orbit` function. 

If `G=rootdatum(:gl,4)` and `pseudo_levi_orbit_reps = reflection_subgroup.(Ref(G),sscentralizer_reps(G))` then `pseudo_levi_orbits = orbits(G,pseudo_levi_orbit_reps)` is incorrect. 

There are five copies of the $A_2$ root system when there should only be four copies. 

Replacing `G=rootdatum(:gl,4)` with `G=coxgroup(:A,3)` does not resolve this issue.