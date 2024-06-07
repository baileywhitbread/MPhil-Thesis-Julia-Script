# The purpose of this script is to compute E-polynomials associated of character varieties
# The important function is EX(G,g,n) which computes the E-polynomial E(X;q) = |X(Fq)|

# Load Jean Michel's package
using Chevie

# Helper functions
function orderpol(G)
	return PermRoot.generic_order(G,Pol(:q))
end

function pi0(L)
	return length(algebraic_center(L).AZ)
end

function subset(L,M)
	return issubset(inclusion(L),inclusion(M))
end

function equal(L,M)
	return subset(L,M) && subset(M,L)
end

function mob(A,B,poset)
	if equal(A,B)
		return 1
	elseif subset(A,B)
		mob_value = 0
		for element in poset
			if subset(A,element) && subset(element,B) && !equal(element,B)
				mob_value += mob(A,element,poset)
			end
		end
		return (-1)*mob_value
	else
		error("First argument must be a subset of the second argument")
	end
end

function nu(L,iplevis,plevis)
	nu_value = 0
	for iplevi in iplevis
		if subset(L,iplevi)
			nu_value += mob(L,iplevi,plevis)*pi0(iplevi)
		end
	end
	return nu_value
end

# Functions for computing E-polynomials

# Compute the types of G and the relevant data
function gptypes(G)
	# This returns a matrix (called gptypes) of data required to compute E-polynomials
	
	# gptypes[i,:][1] = the type [L,ρ] as a tuple
	# gptypes[i,:][2] = |Φ(L)+| of ith type
	# gptypes[i,:][3] = ρ(1) of ith type (unipotent character degree)
	# gptypes[i,:][4] = |L(Fq)| of ith type
	# gptypes[i,:][5] = χᵨ(1) of ith type (Weyl gp character degree)
	# gptypes[i,:][6] = |W(L)| of ith type
	# gptypes[i,:][7] = |[L]| of ith type
	# gptypes[i,:][8] = ν(L) of ith type
	
	# Gather required info about G
	G_dual = rootdatum(simplecoroots(G),simpleroots(G));
	G_positive_root_size = Int(length(roots(G))/2);
	
	# Compute pseudo Levi orbit representatives and pseudo Levi orbits
	plorbit_reps = reflection_subgroup.(Ref(G_dual),sscentralizer_reps(G_dual)); 
	plorbits = orbits(G_dual,plorbit_reps); 
	# A problem: |[G]| > 1 iff G=SO5
	# Therefore we kill the duplicate that occurs
	if G==rootdatum(:so,5)
		plorbits[4] = [plorbits[4][1]];
	end
	
	# Compute all pseudo Levis and all isolated pseudo Levis	
	plevis = [];
	iplevis = [];
	for plorbit in plorbits
		for plevi in plorbit
			append!(plevis,[plevi])
			# Check is pseudo Levi is isolated
			# This happens iff it has the same number of simple roots as G
			if length(gens(plevi)) == length(gens(G))
				append!(iplevis,[plevi])
			end
		end
	end
	
	# Create gptypes
	gptypes = Array{Any}(nothing,0,8);
	# Pick L in the type [L,ρ]
	for plevi in plorbit_reps
		# Compute relevant data for L
		plevi_order = orderpol(plevi);
		plevi_positive_root_size = Int(length(roots(plevi))/2);
		plevi_orbit_size = length(orbit(G_dual,plevi)); 
		plevi_weyl_size = length(plevi);
		plevi_nu = nu(plevi,iplevis,plevis);
		plevi_uc = UnipotentCharacters(plevi);
		plevi_uc_names = charnames(plevi_uc,limit=true);
		plevi_uc_degs = degrees(plevi_uc);
		# Pick unipotent character ρ in the type [L,ρ]
		for i in 1:length(plevi_uc)
			# Check if unipotent character is principal
			# This happens iff q-1 does not divide its degree
			if Int(plevi_uc_degs[i](1))!=0
				# Add row of data for the type [L,ρ]
				gptypes = vcat(gptypes,
				[(plevi,plevi_uc_names[i]) plevi_positive_root_size plevi_uc_degs[i] plevi_order Int(plevi_uc_degs[i](1)) plevi_weyl_size plevi_orbit_size plevi_nu]
				);
			end
		end
	end
	return gptypes
end

function gp_row_term(G,i,genus_num,puncture_num)
	# Returns (|Z(Fq)|/|T(Fq)|^n) * ||tau||(q)^(2g-2+n) * S_tau(q)
	# where tau is the ith type, g = genus_num and n = puncture_num
	
	# Grab necessary data
	G_rank = rank(G);
	G_ssrank = semisimplerank(G);
	T = torus(G_rank); # Split maximal torus of G
	Z = torus(G_rank-G_ssrank); # Centre of G
	G_positive_root_size = Int(length(roots(G))/2); # |Φ(G)+| 
	G_weyl_group_size = length(G); # |W|
	gptype_data = gptypes(G) 
	
	# For readability:
	# gptypes[i,:][1] = the ith type, ie. the pair [L,ρ]
	# gptypes[i,:][2] = |Φ(L)+| for ith type
	# gptypes[i,:][3] = ρ(1) for ith type (unipotent character degree)
	# gptypes[i,:][4] = |L(Fq)| for ith type
	# gptypes[i,:][5] = χᵨ(1) for ith type (Weyl gp character degree)
	# gptypes[i,:][6] = |W(L)| for ith type
	# gptypes[i,:][7] = |[L]| for ith type
	# gptypes[i,:][8] = ν(L) for ith type
	
	coeff = orderpol(Z)//(orderpol(T)^puncture_num);
	tau_q = ((Pol(:q)^(G_positive_root_size-gptype_data[i,:][2]))*gptype_data[i,:][4])//gptype_data[i,:][3];
	S_tau_q = orderpol(Z)*(gptype_data[i,:][5]^puncture_num)*gptype_data[i,:][7]*((G_weyl_group_size//gptype_data[i,:][6])^(puncture_num-1))*gptype_data[i,:][8]
	
	return coeff*(tau_q)^(2*genus_num-2+puncture_num)*S_tau_q
end

function EX(G,genus_num,puncture_num)
	# Returns the E-polynomial E(X;q) = |X(Fq)| associated to (G,g,n)
	gptype_data = gptypes(G);
	epol = Pol([0]);
	for row_number in 1:size(gptype_data)[1]
		epol += gp_row_term(G,row_number,genus_num,puncture_num)
	end
	return Pol{Int64}(epol)
end

# Display human-readable table
function gptable(G)
	gptype_data = gptypes(G)
	clabels = ["|Φ(L)+|","|L(Fq)|","ρ(1)","χᵨ(1)","|W(L)|","|[L]|","ν(L)"];
	rlabels = xrepr.(Ref(rio()),gptype_data[:,1]); # xrepr(rio(),?) is a string of ? when printed on the REPR
	repr_gptype_data = xrepr.(Ref(rio()),gptype_data[:,2:size(gptype_data)[2]]);
	println("A type is a W-orbit [L,ρ] where ")
	println("L is an endoscopy of G")
	println("ρ is a principal unipotent of L(Fq)")
	println("Φ(L)+ is the set of positive roots of L")
	println("|L(Fq)| is the size of L(Fq)")
	println("ρ(1) is the degree of the unipotent character ρ")
	println("χᵨ(1) is the degree of the Weyl group character associated to ρ")
	println("W(L) is the Weyl group of L")
	println("[L] is the orbit of L under the natural W-action")
	println("ν(L) is an integer only depending on L")
	println("")
	return showtable(repr_gptype_data;col_labels=clabels,rows_label="Types [L,ρ]",row_labels=rlabels)
end

# Testing counting functions
function is_palindromic(f)
	return (f(0)!= 0) && (f.c == f.c[end:-1:1])
end

function euler_zero(G,genus_max,puncture_max)
	# Checks if χ(X) is zero or non-zero from g=1,2,..,genus_max and n=1,2,...,puncture_max
	for g in 1:genus_max
		for n in 1:puncture_max
			try 
				if EX(G,g,n)(1)==0
					println("χ(X)=0 when g=",g," and n=",n)
				elseif EX(G,g,n)(1)!=0
					println("χ(X) non-zero when g=",g," and n=",n)
				end
			catch err
				if isa(err,OverflowError)
					println("Overflow error when g=",g," and n=",n)
				else
					println(err," when g=",g," and n=",n)
				end
			end
		end
	end
end