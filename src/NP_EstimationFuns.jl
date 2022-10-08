
# Hazard, Cumulative Hazards and 

hazard(d,r) = d./r

function nelsonaalen(d::Vector,r::Vector) 
	CH = accumulate(+, d./r)
	V = accumulate(+, d./(r.^2))
	return cumu_hazard(CH, V, sqrt.(V), sqrt.(V./(CH.*CH)))
end

function kaplanmeier(d::Vector,r::Vector) 
	S = accumulate(*, (1 .- d./r))
	top =	accumulate(+, (d./r)./(r .- d))
	V = top.*(S.*S)
	term = log.( (r .- d)./r )
	for (i,term_i) in enumerate(term)
		isinf(term_i) ? term[i] = 0. : nothing
	end
	bottom = accumulate(+, term)   # sqrt(1 / log(S)²) = 1 / ∑ln(1 .- d./r)
	sigma = sqrt.(top)./bottom
	return survivor(S, V, sqrt.(V), sigma)
end

hazard(npd::np_survival) = hazard(npd.d, npd.r)
kaplanmeier(npd::np_survival) = kaplanmeier(npd.d, npd.r)
nelsonaalen(npd::np_survival) = nelsonaalen(npd.d, npd.r) 

hazard(df::DataFrame, tvar::Symbol, failvar::Symbol) = hazard(np_survival(df,tvar,failvar))
kaplanmeier(df::DataFrame, tvar::Symbol, failvar::Symbol)  = kaplanmeier(np_survival(df,tvar,failvar))
nelsonaalen(df::DataFrame, tvar::Symbol, failvar::Symbol) = nelsonaalen(np_survival(df,tvar,failvar))

get_ci(km::survivor, alpha::Float64) = [km.surv.*exp.(-quantile(Normal(), (1-alpha/2 )).*km.sigma) km.surv.*exp.(quantile(Normal(), (1-alpha/2 )).*km.sigma)]



