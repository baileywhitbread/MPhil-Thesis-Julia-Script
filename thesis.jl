using Chevie # Load Jean Michel's package
G=coxgroup(:G,2); 
uc=UnipotentCharacters(G);
################################################################################################
G_dual=coxgroup(:G,2); 
pseudo_levi_orbit_reps = reflection_subgroup.(Ref(G_dual),sscentralizer_reps(G_dual));
pseudo_levi_orbits = orbits(G_dual,pseudo_levi_orbit_reps);
################################################################################################
pseudo_levis = [];
isolated_pseudo_levis = [];
for pseudo_levi_orbit in pseudo_levi_orbits
	for pseudo_levi in pseudo_levi_orbit
		append!(pseudo_levis,[pseudo_levi])
		if length(gens(pseudo_levi)) == length(gens(G)) # Isolated iff no. of simples equal
			append!(isolated_pseudo_levis,[pseudo_levi])
		end
	end
end
################################################################################################
function subset(L,M)
	return issubset(inclusion(L),inclusion(M))
end

function equal(L,M)
	return inclusion(L) == inclusion(M)
end
################################################################################################
function pi0(L) 
	return length(algebraic_center(L).AZ) 
end
################################################################################################
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
################################################################################################
function nu(L)
	nu_value = 0
	for isolated_pseudo_levi in isolated_pseudo_levis
		if subset(L,isolated_pseudo_levi)
			nu_value += mob(L,isolated_pseudo_levi,pseudo_levis)*pi0(isolated_pseudo_levi)
		end
	end
	return nu_value
end
################################################################################################
type_data = Array{Any}(nothing,0,8);
for pseudo_levi in pseudo_levi_orbit_reps
	pseudo_levi_order_poly = PermRoot.generic_order(pseudo_levi,Pol(:q));
	pseudo_levi_positive_root_size = Int(length(roots(pseudo_levi))/2);
	pseudo_levi_orbit_size = length(orbit(G_dual,pseudo_levi)); 
	pseudo_levi_weyl_size = length(pseudo_levi);
	pseudo_levi_nu = nu(pseudo_levi);
	pseudo_levi_uc = UnipotentCharacters(pseudo_levi);
	pseudo_levi_uc_names = charnames(pseudo_levi_uc,limit=true);
	pseudo_levi_uc_degree_polys = degrees(pseudo_levi_uc);
	for i in 1:length(pseudo_levi_uc)
		if Int(pseudo_levi_uc_degree_polys[i](1))!=0 # Check unipotent character principal
			type_row = Array{Any}(nothing,1,0);
			global type_row = hcat(type_row,[(pseudo_levi,pseudo_levi_uc_names[i])]);
			global type_row = hcat(type_row,[pseudo_levi_positive_root_size]);
			global type_row = hcat(type_row,[pseudo_levi_uc_degree_polys[i]]);
			global type_row = hcat(type_row,[pseudo_levi_order_poly]);
			global type_row = hcat(type_row,[Int(pseudo_levi_uc_degree_polys[i](1))]);
			global type_row = hcat(type_row,[pseudo_levi_weyl_size]);
			global type_row = hcat(type_row,[pseudo_levi_orbit_size]);
			global type_row = hcat(type_row,[pseudo_levi_nu]);
			global type_data = vcat(type_data,type_row);
		end
	end
end

