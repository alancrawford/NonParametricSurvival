
# Data Setup 

get_bin_counts(df::DataFrame, tvar::Symbol, failvar::Symbol) =
	sort!( combine(groupby(df, [tvar, failvar]), nrow ), tvar )


get_bin_counts(df::DataFrame, tvar::Symbol, failvar::Symbol, wgt_var::Symbol) =
	sort!(combine(groupby(df, [tvar, failvar]), wgt_var => sum => :nrow), tvar)



function np_survival(df::DataFrame, tvar::Symbol, failvar::Symbol; kwargs...)
	kwd = Dict(kwargs)
	d = zeros(maximum(df[!, tvar]))
	m = zeros(maximum(df[!, tvar]))

	if !haskey(kwd, :wgt_var)
		dfc = get_bin_counts(df, tvar, failvar)
	else 
		dfc = get_bin_counts(df, tvar, failvar, kwd[:wgt_var])
	end
	for dfr in eachrow(dfc)
		if dfr[failvar].==1
			d[dfr[tvar]] = dfr[:nrow] 
		else 
			m[dfr[tvar]] = dfr[:nrow] 
		end
	end

	dm = d + m
	cumu_dm = accumulate(+, dm)

	r = [cumu_dm[end]]
	for (n, exits) in enumerate(dm)
		push!(r, r[n] - exits)
	end
	pop!(r)

	t = findall(x->x!=0. , dm)

	return np_survival(t,d[t],m[t],r[t])
end
