# MPhil Thesis Julia Script

This is the Julia script appearing in my MPhil thesis. 

It runs on Julia version 1.10.2 and the current version of Chevie as of 13/06/24.

This script computes representation-theoretic data needed to compute counting polynomials of character varieties.

## 14/06/24
I found a problem with the `orbit` function. 

If `G=rootdatum(:gl,4)` and `pseudo_levi_orbit_reps = reflection_subgroup.(Ref(G),sscentralizer_reps(G))` then `pseudo_levi_orbits = orbits(G,pseudo_levi_orbit_reps)` is incorrect. 

There are five copies of the $A_2$ root system when there should only be four copies. 

Replacing `G=rootdatum(:gl,4)` with `G=coxgroup(:A,3)` does not resolve this issue.

## 25/06/24
This duplication does not occur if $G=GL_1,GL_2,GL_3,SO_5, G_2, F_4, E_6, E_7$. I have not checked $E_8$. 

It does occur if $G=GL_4,\ldots,GL_{10},SO6,SO7,\ldots,SO_{12}$.
