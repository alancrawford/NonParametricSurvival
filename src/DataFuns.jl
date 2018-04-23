
# Data Setup 

function get_bin_counts(df::DataFrame, tvar::Symbol, failvar::Symbol)
	df[:ones] = 1
	cts = sort(
			@based_on(DataFrames.groupby(df, [tvar, failvar]), 
			counts = sum(:ones)),
			cols = order(:time)
		)
end

function np_survival(df::DataFrame, tvar::Symbol, failvar::Symbol)
	d = zeros(maximum(df[tvar]))
	m = zeros(maximum(df[tvar]))

	dfc = get_bin_counts(df, tvar, failvar)
	for dfr in eachrow(dfc)
		if dfr[failvar].==1
			d[dfr[tvar]] = dfr[:counts] 
		else 
			m[dfr[tvar]] = dfr[:counts] 
		end
	end

	dm = d + m
	cumu_dm = accumulate(+, dm)
	desc_dm = maximum(cumu_dm) - cumu_dm

	r = [desc_dm[1]]
	for (n, exits) in enumerate(dm)
		push!(r, r[n] - exits)
	end
	pop!(r)

	t = find(dm)

	return np_survival(t,d[t],m[t],r[t])
end
