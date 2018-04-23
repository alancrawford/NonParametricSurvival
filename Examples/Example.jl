using NonParametricSurvival
using CSV

df = CSV.read("Documents/catheter.csv");
DataFrames.categorical!(df, :female);


npd = np_survival(df, :time, :infect);
sv = kaplanmeier(npd);
ch = nelsonaalen(npd);

# Print Out Survivor Function
[map(x->Int.(x), [npd.t npd.r npd.m npd.d]) sv.surv sv.se ]

# Print Out CH Function
[map(x->Int.(x), [npd.t npd.r npd.m npd.d]) ch.ch ch.se ]

using Plots
scatter(npd.t, ch.ch, line=:step, title="Cumulative Hazard", legend=false)
scatter(npd.t, sv.surv, line=:step, title="Survivor fn.", legend=false)