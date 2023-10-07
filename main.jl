using LinearAlgebra,DataFrames,CSV,Distributions,StatsBase,SparseArrays
using JuMP,Mosek,MosekTools
using StatsPlots,LaTeXStrings,ColorSchemes

# Include helpers
include("src/helpers.jl")
include("src/load_data.jl")
include("src/eval_econ_props.jl")

# Include models
include("src/models/socp_mc.jl")
include("src/models/lp_mc_twostage.jl")
include("src/models/lp_mc_scenarios.jl")


settings = Dict(:Nsteps         => 20,      #Steps for piecewise linearization of quadratic cost
                :ε              => 0.05,    #chance constraints joint violation probability
                :T              => 24,      #market clearing periods
                :S              => 50,      #:S = zero => scenario numbers from CC feasibility requirement, S = num scenarios
                :det            => false,   #det = true => risk parameter = 0
                :drcc           => false,   #drcc = true => distributionally robust cc modeling. ε must be raised for feasibility.
                :histdata       => false,   #histdata = true => estimation of moments from historical data
                :Wcap           => 52*5,    #installed wind capacity for each wind farm
                :ESRs           => true,    #Whether storage is included or not
                :ESRoverride    => true,    #setting this flag to true sizes ESRs as a function of wind penetration
                :ramplimits     => true,    #modeling ramping rate limits for generators
                :σ              => 0.1,     #std. deviation for wind farm forecast errors
                :MRR            => 0.00,    #minimum reserve requirement factor: 0.15 for uncongested, 0.5 for congested
                :bottleneck     => true)   #setting this flag true reduces some line limits to simulate congestion

networkdata         = load_data(settings)
winddata,loaddata   = prep_wind_load_data(settings,networkdata)
# Uncomment the following to compute the number of scnearios to provide CC-like guarantees
#settings[:S]       = compute_num_scenarios(settings,networkdata)

##  FOR SINGLE RUN ONLY
#Deterministic Problem
settings[:det]      = true
sol_det_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
println("Deterministic SOCP Cost => ", sol_det_SOCP[:cost]/1e3)

# # Chance-constrained stochastic SOCP
settings[:det]      = false
sol_sto_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
println("Stochastic SOCP Cost => ", sol_sto_SOCP[:cost]/1e3)

# Compute LMPs and flexibility payoff
sol_sto_SOCP[:Π]     = compute_lmps(settings,networkdata,sol_sto_SOCP)
flex_gens_payoff     = analyze_gen_flex_payoffs(settings,networkdata,sol_det_SOCP,sol_sto_SOCP)
flex_esrs_payoff     = analyze_esr_flex_payoffs(settings,networkdata,sol_det_SOCP,sol_sto_SOCP)

#Compute market outcomes for players
gen_outcomes_SOCP   = process_gen_outcomes(settings,networkdata,sol_sto_SOCP)
esr_outcomes_SOCP   = process_esr_outcomes(settings,networkdata,sol_sto_SOCP)

#Ex-ante OOS constraint violations evaluation
ex_ante_res         = compute_ex_ante_violprob(settings,networkdata,winddata,sol_sto_SOCP)


#Deterministic LP market clearing -- Option 1: Two-stage/Sequential, Deterministic
settings[:bottleneck] == true ? settings[:MRR] = 0.5 : settings[:MRR] = 0.15
sol_det_LP_DA       = det_lp_da_mc(settings,networkdata,winddata,loaddata)
res_realtime_det_LP = DataFrame(iter = Int[], status = Any[], solvetime=Float64[], cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
for scen in 1:winddata[:S]
    sol_det_LP_RT   = det_lp_rt_mc(settings,networkdata,winddata,loaddata,sol_det_LP_DA,scen)
    # Uncomment the lines below to show scenario-wise real-time operation outcomes
    # println("Scen =>", scen)
    # println("LP Realtime Cost => ", round(sol_det_LP_RT[:cost],digits=2))
    # println("Load shed =>",         round(sum(sol_det_LP_RT[:p_shed]),digits=2))
    # println("Wind spilled =>",  round(sum(sol_det_LP_RT[:w_spill]),digits=2))
    # println("==========")
    push!(res_realtime_det_LP, [scen, sol_det_LP_RT[:status], sol_det_LP_RT[:solvetime], round(sol_det_LP_RT[:cost],digits=2), round(sum(sol_det_LP_RT[:p_shed]),digits=2), round(sum(sol_det_LP_RT[:w_spill]), digits=2)])
end
@show res_realtime_det_LP
println("Deterministic LP Cost => ", (mean(res_realtime_det_LP[!, :cost_act]) + sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr])/1e3)
gen_outcomes_det_LP = process_gen_outcomes(settings,networkdata,sol_det_LP_DA)
esr_outcomes_det_LP = process_esr_outcomes(settings,networkdata,sol_det_LP_DA)

#LP market clearing -- Option 2: Stochastic LP
sol_sto_LP          = sto_lp_mc(settings,networkdata,winddata,loaddata)
println("Stochastic LP Cost => ", sol_sto_LP[:cost]/1e3)
sol_sto_LP[:Π]      = compute_lmps(settings,networkdata,sol_sto_LP)
gen_outcomes_LP     = process_gen_outcomes(settings,networkdata,sol_sto_LP)
esr_outcomes_LP     = process_esr_outcomes(settings,networkdata,sol_sto_LP)