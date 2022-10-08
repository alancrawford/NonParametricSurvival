using NonParametricSurvival
using CSV, DataFrames, CategoricalArrays

df = CSV.read("/Users/alancrawford/Documents/catheter.csv", DataFrame);
df.female = categorical(df.female);

df[!, :ones] .= 1;
npd = np_survival(df, :time, :infect);

sv = kaplanmeier(npd);
ch = nelsonaalen(npd);

# Print Out Survivor Function
[ npd.t npd.r npd.m npd.d sv.surv sv.se ]

# Print Out CH Function
[ npd.t npd.r npd.m npd.d ch.ch ch.se ]

using Plots
scatter(npd.t, ch.ch, line=:step, title="Cumulative Hazard", legend=false)
scatter(npd.t, sv.surv, line=:step, title="Survivor fn.", legend=false)


df[!, :rd] = rand(nrow(df));
npdr = np_survival(df, :time, :infect; wgt_var=:rd);
svr = kaplanmeier(npdr);
[npdr.t npdr.r npdr.m npdr.d svr.surv svr.se ]


# ------------------------------------------- #
