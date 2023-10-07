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


# Setup results folders and subfolders and create them
R_FILES_DIR_NAME = "results/r_files/"
mkpath(R_FILES_DIR_NAME)
RESULTS_DATA_DIR_NAME = "results/data/"
mkpath(RESULTS_DATA_DIR_NAME)
RESULTS_PLOTS_DIR_NAME = "results/plots/"
mkpath(RESULTS_PLOTS_DIR_NAME)


settings = Dict(:Nsteps         => 20,      #Steps for piecewise linearization of quadratic cost
                :ε              => 0.05,    #chance constraints joint violation probability
                :T              => 24,      #market clearing periods
                :S              => 2,     #:S = zero => scenario numbers from CC feasibility requirement, S = num scenarios
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
    # println("Scen =>", scen)
    # println("LP Realtime Cost => ", round(sol_det_LP_RT[:cost],digits=2))
    # println("Load shed =>",         round(sum(sol_det_LP_RT[:p_shed]),digits=2))
    # println("Wind spilled =>",  round(sum(sol_det_LP_RT[:w_spill]),digits=2))
    # println("==========")
    push!(res_realtime_det_LP, [scen, sol_det_LP_RT[:status], sol_det_LP_RT[:solvetime], round(sol_det_LP_RT[:cost],digits=2), round(sum(sol_det_LP_RT[:p_shed]),digits=2), round(sum(sol_det_LP_RT[:w_spill]), digits=2)])
end
#@show res_realtime_det_LP
println("Deterministic LP Cost => ", (mean(res_realtime_det_LP[!, :cost_act]) + sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr])/1e3)
gen_outcomes_det_LP = process_gen_outcomes(settings,networkdata,sol_det_LP_DA)
esr_outcomes_det_LP = process_esr_outcomes(settings,networkdata,sol_det_LP_DA)

#LP market clearing -- Option 2: Stochastic LP
sol_sto_LP          = sto_lp_mc(settings,networkdata,winddata,loaddata)
println("Stochastic LP Cost => ", sol_sto_LP[:cost]/1e3)
sol_sto_LP[:Π]      = compute_lmps(settings,networkdata,sol_sto_LP)
gen_outcomes_LP     = process_gen_outcomes(settings,networkdata,sol_sto_LP)
esr_outcomes_LP     = process_esr_outcomes(settings,networkdata,sol_sto_LP)


## Payment Rate Evaluation experiments --- ##
rescapvals = 0:5:30
paymentrate=[]
for i in 1:length(rescapvals)
    #SOLVE DETERMINISTIC PROBLEM
    settings[:det]      = true
    sol_det_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
    println("Deterministic SOCP Cost => ", sol_det_SOCP[:cost]/1e3)

    #SOLVE STOCHASTIC PROBLEM
    networkdata[:gens][6].r̅ = rescapvals[i]
    settings[:det]      = false
    sol_sto_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
    println("Stochastic SOCP Cost => ", sol_sto_SOCP[:cost]/1e3)

    #COMPUTE LMPS AND FLEXIBILITY PAYOFF
    sol_sto_SOCP[:Π]     = compute_lmps(settings,networkdata,sol_sto_SOCP)
    flex_gens_payoff     = analyze_gen_flex_payoffs(settings,networkdata,sol_det_SOCP,sol_sto_SOCP)
    push!(paymentrate, flex_gens_payoff[:flex_pay_rate][6,11])
end

### Network price plots
xy_data = DataFrame(CSV.File("data/nodes_coords.csv"))
out_data = DataFrame(node=Int[],x=Any[],y=Any[],ntype=Any[],el_price=Any[])
for n in 1:length(networkdata[:buses])
    hour=23
    push!(out_data, [n,xy_data[!,:x][n],xy_data[!,:y][n], xy_data[!, :ntype][n], sol_sto_SOCP[:Π][:Π_E][n,hour]])
end
CSV.write(joinpath(R_FILES_DIR_NAME, "out_data_24node.csv"), out_data)

kde_settings = Dict(:gridsize => 0.1, :bias => 1)
kerData = prepare_kde_data(kde_settings)
CSV.write(joinpath(R_FILES_DIR_NAME, "price_kernel_plot.csv"), kerData)



## - Data for In-sample Day-ahead Table --- ##
#No bottleneck case
tabledata = DataFrame(Model=Any[],ExpCost=Any[], CompTime=Any[], CostEner=Any[], CostFlex=Any[], GensMW=Any[], ESRChgMW=Any[], ESRDisMW=Any[])
push!(tabledata,["Mcc",
                sol_sto_SOCP[:cost],
                sol_sto_SOCP[:solvetime],
                sol_sto_SOCP[:cost_ener],
                sol_sto_SOCP[:cost_flex],
                sum(sol_sto_SOCP[:p]),
                sum(sol_sto_SOCP[:b][sol_sto_SOCP[:b] .> 0]),
                sum(sol_sto_SOCP[:b][sol_sto_SOCP[:b] .< 0])])

push!(tabledata,["R1",
                sol_det_LP_DA[:cost],
                sol_det_LP_DA[:solvetime],
                sol_det_LP_DA[:cost_res_gen],
                sol_det_LP_DA[:cost_res_esr],
                sum(sol_det_LP_DA[:p]),
                sum(sol_det_LP_DA[:b][sol_det_LP_DA[:b] .> 0]),
                sum(sol_det_LP_DA[:b][sol_det_LP_DA[:b] .< 0])])

push!(tabledata,["R2",
                sol_sto_LP[:cost],
                sol_sto_LP[:solvetime],
                sol_sto_LP[:cost_ener],
                sol_sto_LP[:cost_flex],
                sum(sol_sto_LP[:p_da]),
                sum(sol_sto_LP[:b_da][sol_sto_LP[:b_da] .> 0]),
                sum(sol_sto_LP[:b_da][sol_sto_LP[:b_da] .< 0])])
@show tabledata

## Export Data for Mcc -> compare bottlenecks and no bottlenecks
st = 4
en = 12
settings[:bottleneck]   = false
settings[:det]          = false
networkdata         = load_data(settings)
winddata,loaddata   = prep_wind_load_data(settings,networkdata)

settings[:det]      = true
sol_det_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
println("Deterministic SOCP Cost => ", sol_det_SOCP[:cost]/1e3)

settings[:det]      = false
sol_sto_SOCP            = sto_socp_mc(settings,networkdata,winddata,loaddata)
println("Stochastic SOCP Cost => ", sol_sto_SOCP[:cost]/1e3)

sol_sto_SOCP[:Π]     = compute_lmps(settings,networkdata,sol_sto_SOCP)
flex_gens_payoff     = analyze_gen_flex_payoffs(settings,networkdata,sol_det_SOCP,sol_sto_SOCP)
flex_esrs_payoff     = analyze_esr_flex_payoffs(settings,networkdata,sol_det_SOCP,sol_sto_SOCP)


socp_gen_adj_poli_no_bottlenecks  = DataFrame([st:en sol_sto_SOCP[:α][:,st:en]'], [:hours, :F1, :F2, :F3, :F4, :F5, :F6, :F7, :F8, :F9, :F10, :F11, :F12])
socp_gen_nom_disp_no_bottlenecks  = DataFrame([st:en sol_sto_SOCP[:p][:,st:en]'], [:hours, :F1, :F2, :F3, :F4, :F5, :F6, :F7, :F8, :F9, :F10, :F11, :F12])
socp_esr_adj_poli_no_bottlenecks  = DataFrame([st:en sol_sto_SOCP[:γ][:,st:en]'], [:hours, :S1, :S2, :S3])
socp_esr_nom_disp_no_bottlenecks  = DataFrame([st:en sol_sto_SOCP[:b][:,st:en]'], [:hours, :S1, :S2, :S3])
CSV.write(joinpath(RESULTS_DATA_DIR_NAME,"socp_gen_adj_poli_no_bottlenecks.csv"), socp_gen_adj_poli_no_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_gen_nom_disp_no_bottlenecks.csv"), socp_gen_nom_disp_no_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_esr_adj_poli_no_bottlenecks.csv"), socp_esr_adj_poli_no_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_esr_nom_disp_no_bottlenecks.csv"), socp_esr_nom_disp_no_bottlenecks)

socp_gen_flex_payment_no_bottlenecks = DataFrame([st:en flex_gens_payoff[:flex_pay][:,st:en]'],         [:hours, :F1, :F2, :F3, :F4, :F5, :F6, :F7, :F8, :F9, :F10, :F11, :F12])
socp_gen_flex_payrate_no_bottlenecks = DataFrame([st:en flex_gens_payoff[:flex_pay_rate][:,st:en]'],    [:hours, :F1, :F2, :F3, :F4, :F5, :F6, :F7, :F8, :F9, :F10, :F11, :F12])
socp_esr_flex_payment_no_bottlenecks = DataFrame([st:en flex_esrs_payoff[:flex_pay][:,st:en]'],         [:hours, :S1, :S2, :S3])
socp_esr_flex_payrate_no_bottlenecks = DataFrame([st:en flex_esrs_payoff[:flex_pay_rate][:,st:en]'],    [:hours, :S1, :S2, :S3])
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_gen_flex_payment_no_bottlenecks.csv"), socp_gen_flex_payment_no_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_gen_flex_payrate_no_bottlenecks.csv"), socp_gen_flex_payrate_no_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_esr_flex_payment_no_bottlenecks.csv"), socp_esr_flex_payment_no_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_esr_flex_payrate_no_bottlenecks.csv"), socp_esr_flex_payrate_no_bottlenecks)



settings[:bottleneck]   = true
settings[:det]      = true
sol_det_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
println("Deterministic SOCP Cost => ", sol_det_SOCP[:cost]/1e3)


settings[:det]          = false
networkdata         = load_data(settings)
winddata,loaddata   = prep_wind_load_data(settings,networkdata)
sol_sto_SOCP            = sto_socp_mc(settings,networkdata,winddata,loaddata)
println("Stochastic SOCP Cost => ", sol_sto_SOCP[:cost]/1e3)

sol_sto_SOCP[:Π]     = compute_lmps(settings,networkdata,sol_sto_SOCP)
flex_gens_payoff     = analyze_gen_flex_payoffs(settings,networkdata,sol_det_SOCP,sol_sto_SOCP)
flex_esrs_payoff     = analyze_esr_flex_payoffs(settings,networkdata,sol_det_SOCP,sol_sto_SOCP)

socp_gen_adj_poli_with_bottlenecks  = DataFrame([st:en sol_sto_SOCP[:α][:,st:en]'], [:hours, :F1, :F2, :F3, :F4, :F5, :F6, :F7, :F8, :F9, :F10, :F11, :F12])
socp_gen_nom_disp_with_bottlenecks  = DataFrame([st:en sol_sto_SOCP[:p][:,st:en]'], [:hours, :F1, :F2, :F3, :F4, :F5, :F6, :F7, :F8, :F9, :F10, :F11, :F12])
socp_esr_adj_poli_with_bottlenecks  = DataFrame([st:en sol_sto_SOCP[:γ][:,st:en]'], [:hours, :S1, :S2, :S3])
socp_esr_nom_disp_with_bottlenecks  = DataFrame([st:en sol_sto_SOCP[:b][:,st:en]'], [:hours, :S1, :S2, :S3])
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_gen_adj_poli_with_bottlenecks.csv"), socp_gen_adj_poli_with_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_gen_nom_disp_with_bottlenecks.csv"), socp_gen_nom_disp_with_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_esr_adj_poli_with_bottlenecks.csv"), socp_esr_adj_poli_with_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_esr_nom_disp_with_bottlenecks.csv"), socp_esr_nom_disp_with_bottlenecks)


socp_gen_flex_payment_with_bottlenecks = DataFrame([st:en flex_gens_payoff[:flex_pay][:,st:en]'],         [:hours, :F1, :F2, :F3, :F4, :F5, :F6, :F7, :F8, :F9, :F10, :F11, :F12])
socp_gen_flex_payrate_with_bottlenecks = DataFrame([st:en flex_gens_payoff[:flex_pay_rate][:,st:en]'],    [:hours, :F1, :F2, :F3, :F4, :F5, :F6, :F7, :F8, :F9, :F10, :F11, :F12])
socp_esr_flex_payment_with_bottlenecks = DataFrame([st:en flex_esrs_payoff[:flex_pay][:,st:en]'],         [:hours, :S1, :S2, :S3])
socp_esr_flex_payrate_with_bottlenecks = DataFrame([st:en flex_esrs_payoff[:flex_pay_rate][:,st:en]'],    [:hours, :S1, :S2, :S3])
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_gen_flex_payment_with_bottlenecks.csv"), socp_gen_flex_payment_with_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_gen_flex_payrate_with_bottlenecks.csv"), socp_gen_flex_payrate_with_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_esr_flex_payment_with_bottlenecks.csv"), socp_esr_flex_payment_with_bottlenecks)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "socp_esr_flex_payrate_with_bottlenecks.csv"), socp_esr_flex_payrate_with_bottlenecks)



### Flexibility payments plot
R1data = []
for t in 1:settings[:T]
    push!(R1data, sum(networkdata[:gens][g].c_r * sol_det_LP_DA[:r_p][g,t] for g in 1:length(networkdata[:gens]))
                 +sum(networkdata[:esrs][s].c_r * sol_det_LP_DA[:r_b][s,t] for s in 1:length(networkdata[:esrs])))
    # push!(R2data, ((1/winddata[:S])*(sum(sum(sum(pwl_cost_approx(settings,networkdata,g,0)[:slopes] .* sol_sto_LP[:p_rt][g,:,t,scen]) for g in 1:length(networkdata[:gens])) for scen in 1:winddata[:S]))
    #               + (1/winddata[:S])*sum(sum(networkdata[:esrs][s].c_l * sol_sto_LP[:b_rt][s,t,scen] for s in 1:length(networkdata[:esrs])) for scen in 1:winddata[:S])))
end
flex_payments_data = DataFrame([1:24 sol_sto_SOCP[:λ_R] R1data], [:hours, :Mcc, :R1])
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "compare_socp_R1.csv"), flex_payments_data)

## --- BASIC PLOTS FOR PRICING ANALYSIS --- ###
#Plotting Flexibility provision - alpha
pyplot()
p1 = groupedbar(sol_sto_SOCP[:α]',
                legend=:none,
                bar_position=:stack,
                #label = ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9" "F10" "F11" "F12"],
                palette =:Paired_12,
                title = "Adjustment policies",
                box=:on)

p2 = groupedbar(sol_sto_SOCP[:p]',
                legend=:outerleft,
                bar_position=:stack,
                label = ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9" "F10" "F11" "F12"],
                palette =:Paired_12,
                title = "Nominal power",
                box=:on)

p3 = plot(sum(winddata[:ŵ],dims=1)',
            lw=4,
            label="Wind",
            legend=:topleft,
            #ylim=(900,2700)
            )
p3 = plot!(twinx(),loaddata[:D],
            lw=4,
            colour=[RGB(0,0,0)],
            label="Load",
            box = :on,
            grid = :off,
            legend=:topright,
            #ylim=(900,2700)
            )
plot(p1,p2,p3, layout = (3,1), size=(1000,800))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "with_storage_with_bottleneck_50pcRES_dispatch.pdf"))

#Compare deterministic and stochastic solution in Dispatch
pyplot()
p1 = groupedbar(sol_det_SOCP[:p]',
                legend=:none,
                bar_position=:stack,
                #label = ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9" "F10" "F11" "F12"],
                palette =:Paired_12,
                title = "Deterministic Solution",
                box=:on)

p2 = groupedbar(sol_sto_SOCP[:p]',
                legend=:outerleft,
                bar_position=:stack,
                label = ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9" "F10" "F11" "F12"],
                palette =:Paired_12,
                title = "Stochastic Solution",
                box=:on)
plot(p1,p2, layout = (2,1), size=(1000,800))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "with_storage_with_bottleneck_50pcRES_compare_det_sto.pdf"))
#Obtain flexibility payments and flexibility pay rates
pyplot()
p1 = groupedbar(flex_gens_payoff[:flex_pay]',
                legend=:outerleft,
                bar_position=:stack,
                label = ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9" "F10" "F11" "F12"],
                palette =:Paired_12,
                title="Hourly Flexibility Payment",
                box=:on)
p1 = plot!(sol_sto_SOCP[:λ_R],
            lw=4,
            ls=:dash,
            color= RGB(1,87/255,34/255),
            label="Total Payment")
p2 = groupedbar(flex_gens_payoff[:flex_pay_rate]',
                legend=:outerleft,
                bar_position=:stack,
                label = ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9" "F10" "F11" "F12"],
                palette =:Paired_12,
                title="Flexibility Payment Rate",
                box=:on)
plot(p1,p2, layout = (2,1), size=(1000,800))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "with_storage_with_bottleneck_50pcRES_flex_rate_payments.pdf"))


pyplot()
p1 = groupedbar(flex_esrs_payoff[:flex_pay]',
                legend=:top,
                bar_position=:stack,
                label = ["S1" "S2" "S3"],
                palette =:Paired_12,
                title="Hourly Flexibility Payment",
                box=:on)
p1 = plot!(sol_sto_SOCP[:λ_R],
            lw=4,
            ls=:dash,
            color= RGB(1,87/255,34/255),
            label="Total Payment")
p2 = groupedbar(flex_esrs_payoff[:flex_pay_rate]',
                legend=:none,
                bar_position=:stack,
                label = ["S1" "S2" "S3"],
                palette =:Paired_12,
                title="Flexibility Payment Rate",
                box=:on)
p3 = groupedbar(sol_sto_SOCP[:γ]',
                legend=:none,
                bar_position=:stack,
                ylim=(0,1),
                #label = ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9" "F10" "F11" "F12"],
                palette =:Paired_12,
                title = "Adjustment policies",
                box=:on)

p4 = groupedbar(sol_sto_SOCP[:b]',
                legend=:none,
                bar_position=:stack,
                #label = ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9" "F10" "F11" "F12"],
                palette =:Paired_12,
                title = "Nominal power",
                box=:on)

plot(p1,p2,p3,p4,layout = (2,2), size=(1000,800))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "ESRS_with_bottleneck_50pcRES_flex_rate_payments_dispatch.pdf"))




## -- Graph 1: Dispatch comparison of three models -- ##
#System price
pyplot()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
        "font.size" => 10,
        "axes.labelweight" => "bold",
        "axes.labelsize" => 10,
        "xtick.labelsize" => 8,
        "ytick.labelsize" => 8,
        "legend.fontsize" => 10,
)
merge!(rcParams, font0)
p1 = plot([sol_sto_SOCP[:λ_E] sol_det_LP_DA[:λ_E] sol_sto_LP[:λ_da_sys]], lw=4,
        label=["Mcc" "R1" "R2"],
        xlabel="Hours",
        ylabel="EL Price (\$/MWh)",
        colour= [RGB(7/255,98/255,86/255) RGB(252/255,69/255,45/255) RGB(185/255,199/255,204/255)],
        legend=:bottom, box=:on)
p2 = plot(loaddata[:D] .- sum(winddata[:ŵ],dims=1)',lw=4,label="Net demand",legend=:topleft)
p2 = plot!(twinx(), loaddata[:D],lw=4,colour=[RGB(0,0,0)],label="Load",box = :on, grid = :off)
#Flexible generators
p3 = groupedbar([sum(sol_sto_SOCP[:p],dims=1); sum(sum(sol_det_LP_DA[:p][:,seg,:] for seg in 1:settings[:Nsteps]), dims=1) ; sum(sum(sol_sto_LP[:p_da][:,seg,:] for seg in 1:settings[:Nsteps]),dims=1) + sum(mean(sum(sol_sto_LP[:p_rt][:,seg,:,:] for seg in 1:settings[:Nsteps]),dims=3)[:,:,1],dims=1)]',
    xlim=[5,15],
    label=["Mcc" "R1" "R2"],
    xlabel="Hours",
    ylabel="Generators (MW)",
    colour= [RGB(7/255,98/255,86/255) RGB(252/255,69/255,45/255) RGB(185/255,199/255,204/255)],
    legend=:topleft)

    ticklabel = ["10%" "20%" "30%" "40%" "50%" "60%"]
    p3 = groupedbar([flex_α flex_γ]./24,
                    bar_position = :stack,
                    barwidth=0.7,
                    xticks=(1:6, ticklabel),
                    label=["Flex Gens" "ESRs"],
                    title="Adjustment policies allocation",
                    legend= :outerright)
    savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "fig_flexibility_dispatch.pdf"))

#Energy storage
p4 = groupedbar([sum(sol_sto_SOCP[:b],dims=1); sum(sol_det_LP_DA[:b], dims=1); sum(sol_sto_LP[:b_da],dims=1) + sum(mean(sol_sto_LP[:b_rt],dims=3)[:,:,1],dims=1)]',
    xlim=[5,15],
    label=["Mcc" "R1" "R2"],
    xlabel="Hours",
    ylabel="Energy storage (MW)",
    colour= [RGB(7/255,98/255,86/255) RGB(252/255,69/255,45/255) RGB(185/255,199/255,204/255)],
    legend=:none)
plot(p2,p3,p4,p1, layout = (1,4), size=(1000,200))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "uncongested_dispatch_comparison_LP_SOCP.pdf"))

#Export output for plotting in pgfplots
ELpricesOutput  = DataFrame([1:24 sol_sto_SOCP[:λ_E] sol_det_LP_DA[:λ_E] sol_sto_LP[:λ_da_sys]], [:hours, :Mcc,  :R1, :R2])
GensOutput      = DataFrame([(1:24)'; sum(sol_sto_SOCP[:p],dims=1); sum(sum(sol_det_LP_DA[:p][:,seg,:] for seg in 1:settings[:Nsteps]), dims=1) ; sum(sum(sol_sto_LP[:p_da][:,seg,:] for seg in 1:settings[:Nsteps]),dims=1) + sum(mean(sum(sol_sto_LP[:p_rt][:,seg,:,:] for seg in 1:settings[:Nsteps]),dims=3)[:,:,1],dims=1)]',  [:hours, :Mcc,  :R1, :R2])
ESROutput       = DataFrame([(1:24)'; sum(sol_sto_SOCP[:b],dims=1); sum(sol_det_LP_DA[:b], dims=1); sum(sol_sto_LP[:b_da],dims=1) + sum(mean(sol_sto_LP[:b_rt],dims=3)[:,:,1],dims=1)]',  [:hours, :Mcc,  :R1, :R2])
NetDemandOutput = DataFrame([1:24 loaddata[:D] .- sum(winddata[:ŵ],dims=1)'], [:hours, :netdemand])
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "elprices_output.csv"), ELpricesOutput)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "gens_output.csv"), GensOutput)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "esrs_output.csv"), ESROutput)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "netdemand_output.csv"), NetDemandOutput)


## -- Graph 2: Cost comparison of three models over different renewable penetration shares -- ##
Wcap_vals = 52:52:52*6  #Range of wind penetration shares to evaluate
settings[:S] = 100       #Number of in-sample scenarios
cost_wcap = DataFrame(model = Any[],congestion=Any[],resShare=Float64[], costTot=Float64[],costEne=Float64[],costFlex=Float64[],compTime=Float64[])
for i in 1:length(Wcap_vals)
    #-- Specific settings for this run --#
    settings[:Wcap] = Wcap_vals[i]  #Setting wind capacity value to reach desired res penetration
    settings[:bottleneck] = false   #Network is uncongested
    networkdata         = load_data(settings)   #reload data to size the ESR properly
    winddata,loaddata   = prep_wind_load_data(settings,networkdata) #Get wind forecast and scenarios
    #SOCP MC : Mcc
    sol_sto_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
    push!(cost_wcap, ["Mcc", "acongested",
                      100*round(sum(winddata[:ŵ])/sum(loaddata[:D]),digits=2),
                      sol_sto_SOCP[:cost]/1e3,
                      sol_sto_SOCP[:cost_ener]/1e3,
                      sol_sto_SOCP[:cost_flex]/1e3,
                      sol_sto_SOCP[:solvetime]])
    #Deterministic with MRR : R1
    #settings[:MRR]      = determine_MRR(settings,networkdata,winddata,loaddata)
    settings[:bottleneck] == true ? settings[:MRR] = 0.5 : settings[:MRR] = 0.15
    sol_det_LP_DA       = det_lp_da_mc(settings,networkdata,winddata,loaddata)
    res_realtime_det_LP = DataFrame(iter = Int[], status = Any[], solvetime=Float64[], cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
    for scen = 1:winddata[:S]
        sol_det_LP_RT    = det_lp_rt_mc(settings,networkdata,winddata,loaddata,sol_det_LP_DA,scen)
        push!(res_realtime_det_LP, [scen, sol_det_LP_RT[:status], sol_det_LP_RT[:solvetime], round(sol_det_LP_RT[:cost],digits=2), round(sum(sol_det_LP_RT[:p_shed]),digits=2), round(sum(sol_det_LP_RT[:w_spill]), digits=2)])
    end
    push!(cost_wcap, ["R1","acongested",
                      100*round(sum(winddata[:ŵ])/sum(loaddata[:D]),digits=2),
                      (mean(res_realtime_det_LP[!, :cost_act]) + sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr])/1e3,
                      mean(res_realtime_det_LP[!, :cost_act])/1e3,
                      (sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr])/1e3,
                      sol_det_LP_DA[:solvetime] + mean(res_realtime_det_LP[!, :solvetime])])
    #LP Stochastic : R2
    sol_sto_LP          = sto_lp_mc(settings,networkdata,winddata,loaddata)
    push!(cost_wcap, ["R2","acongested",
                      100*round(sum(winddata[:ŵ])/sum(loaddata[:D]),digits=2),
                      sol_sto_LP[:cost]/1e3,
                      sol_sto_LP[:cost_ener]/1e3,
                      sol_sto_LP[:cost_flex]/1e3,
                      sol_sto_LP[:solvetime]])
end
@show cost_wcap

#CONGESTED LINES:
for i in 1:length(Wcap_vals)
    #-- Specific settings for this run --#
    settings[:Wcap] = Wcap_vals[i]  #Setting wind capacity value to reach desired res penetration
    settings[:bottleneck] = true   #Network is congested
    networkdata         = load_data(settings)   #reload data to size the ESR properly
    winddata,loaddata   = prep_wind_load_data(settings,networkdata) #Get wind forecast and scenarios
    #SOCP MC : Mcc
    sol_sto_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
    push!(cost_wcap, ["Mcc", "congested",
                      100*round(sum(winddata[:ŵ])/sum(loaddata[:D]),digits=2),
                      sol_sto_SOCP[:cost]/1e3,
                      sol_sto_SOCP[:cost_ener]/1e3,
                      sol_sto_SOCP[:cost_flex]/1e3,
                      sol_sto_SOCP[:solvetime]])
    #Deterministic with MRR : R1
    settings[:bottleneck] == true ? settings[:MRR] = 0.5 : settings[:MRR] = 0.15
    sol_det_LP_DA       = det_lp_da_mc(settings,networkdata,winddata,loaddata)
    res_realtime_det_LP = DataFrame(iter = Int[], status = Any[], solvetime=Float64[], cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
    for scen = 1:winddata[:S]
        sol_det_LP_RT    = det_lp_rt_mc(settings,networkdata,winddata,loaddata,sol_det_LP_DA,scen)
        push!(res_realtime_det_LP, [scen, sol_det_LP_RT[:status], sol_det_LP_RT[:solvetime], round(sol_det_LP_RT[:cost],digits=2), round(sum(sol_det_LP_RT[:p_shed]),digits=2), round(sum(sol_det_LP_RT[:w_spill]), digits=2)])
    end
    push!(cost_wcap, ["R1","congested",
                      100*round(sum(winddata[:ŵ])/sum(loaddata[:D]),digits=2),
                      (mean(res_realtime_det_LP[!, :cost_act]) + sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr])/1e3,
                      mean(res_realtime_det_LP[!, :cost_act])/1e3,
                      (sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr])/1e3,
                      sol_det_LP_DA[:solvetime] + mean(res_realtime_det_LP[!, :solvetime])])
    #LP Stochastic : R2
    sol_sto_LP          = sto_lp_mc(settings,networkdata,winddata,loaddata)
    push!(cost_wcap, ["R2","congested",
                      100*round(sum(winddata[:ŵ])/sum(loaddata[:D]),digits=2),
                      sol_sto_LP[:cost]/1e3,
                      sol_sto_LP[:cost_ener]/1e3,
                      sol_sto_LP[:cost_flex]/1e3,
                      sol_sto_LP[:solvetime]])

end
@show cost_wcap
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "cost_vs_wcap_100scen.csv"), cost_wcap)


## -- Graph: 3 -- SOCP Prices with RES penetration --- ##
Wcap_vals = 52:52:52*6
settings[:S] = 1000
λ_E_vals = []
λ_R_vals = []
flex_α = zeros(length(Wcap_vals))
flex_γ = zeros(length(Wcap_vals))
for i in 1:length(Wcap_vals)
    settings[:Wcap] = Wcap_vals[i]
    networkdata         = load_data(settings)
    networkdata[:ncon]  = (settings[:T]*length(networkdata[:gens])*4 +          #Flex gens reserve limits and production Limits
                            (settings[:T]-1)*length(networkdata[:gens])*2 +     #Flex gens ramping limits
                                settings[:T]*length(networkdata[:lines])*2 +    #Line flow limits
                                    settings[:T]*length(networkdata[:esrs])*4 + #Storage operation power limits
                                        2*length(networkdata[:esrs]))           #End of day storage bounds
    networkdata[:nvar] = settings[:T]*(length(networkdata[:gens]) + length(networkdata[:esrs]))*2
    winddata,loaddata   = prep_wind_load_data(settings,networkdata)
    sol_sto_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
    push!(λ_E_vals, sol_sto_SOCP[:λ_E])
    push!(λ_R_vals, sol_sto_SOCP[:λ_R])
    flex_α[i] = sum(sol_sto_SOCP[:α])
    flex_γ[i] = sum(sol_sto_SOCP[:γ])
end
#Export results
ELPricesMccWithRES = DataFrame([1:24 λ_E_vals[1] λ_E_vals[2] λ_E_vals[3] λ_E_vals[4] λ_E_vals[5] λ_E_vals[6]],
                                [:hours, :res10, :res20, :res30, :res40, :res50, :res60])
REPricesMccWithRES = DataFrame([1:24 λ_R_vals[1] λ_R_vals[2] λ_R_vals[3] λ_R_vals[4] λ_R_vals[5] λ_R_vals[6]],
                                [:hours, :res10, :res20, :res30, :res40, :res50, :res60])
FlexAllocationWithRES = DataFrame([10:10:60 flex_α flex_γ], [:res, :gens, :esrs])
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "el_prices_with_res.csv"), ELPricesMccWithRES)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "re_prices_with_res.csv"), REPricesMccWithRES)
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "flex_alloc_with_res.csv"), FlexAllocationWithRES)

p1 = plot(λ_E_vals, lw=4,
        label=["10%" "20%" "30%" "40%" "50%" "60%"],
        legendtitle="RES share",
        xlabel="Hours",
        ylabel="EL Price (\$/MWh)",
        legend=:outerright)

p2 = plot(λ_R_vals, lw=4,
        label=["10%" "20%" "30%" "40%" "50%" "60%"],
        legendtitle="RES share",
        xlabel="Hours",
        ylabel="Flexibility Payment (\$)",
        legend=:outerright)

plot(p1,p2, layout=(2,1))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "fig_price_SOCP_w_windcap.pdf"))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "fig_flexibility_dispatch.pdf"))

## --- Cost Recovery Plots --- ##
pyplot()
p1 = groupedbar(["F1","F1","F2","F2","F3","F3","F4","F4","F5","F5","F6","F6","F7","F7","F8","F8","F9","F9","F10","F10","F11","F11","F12","F12"],
            [gen_outcomes_SOCP[:Rev_gen] gen_outcomes_LP[:Rev_gen]],
            label=["Mcc" "R2"],
            xlabel = "Generator #",
            ylabel="Revenues",
            title="Revenues: Generators",
            colour= [RGB(0,121/255,107/255) RGB(199/255,211/255,215/255)],
            legend=:outerright)

p2 = groupedbar(["F1","F1","F2","F2","F3","F3","F4","F4","F5","F5","F6","F6","F7","F7","F8","F8","F9","F9","F10","F10","F11","F11","F12","F12"],
            [gen_outcomes_SOCP[:Cos_gen] gen_outcomes_LP[:Cos_gen]],
            label=["Mcc" "R2"],
            xlabel = "Generator #",
            ylabel="Costs",
            title="Costs: Generators",
            colour= [RGB(0,121/255,107/255) RGB(199/255,211/255,215/255)],
            legend=:outerright)

p3 = groupedbar(["F1","F1","F2","F2","F3","F3","F4","F4","F5","F5","F6","F6","F7","F7","F8","F8","F9","F9","F10","F10","F11","F11","F12","F12"],
            [gen_outcomes_SOCP[:Pro_gen] gen_outcomes_LP[:Pro_gen]],
            label=["Mcc" "R2"],
            xlabel = "Generator #",
            ylabel="Profits",
            title="Profits: Generators",
            colour= [RGB(0,121/255,107/255) RGB(199/255,211/255,215/255)],
            legend=:outerright)
plot(p1,p2,p3, layout = (3,1))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "fig_cost_recovery_generators.pdf"))


## Expt 1: ε vs SOCP cost
ε_vals = 0.01:0.02:0.4
cost_violprob = DataFrame(SOCP=Float64[], viol_prob = Float64[])
for i in 1:length(ε_vals)
    settings[:ε] = ε_vals[i]
    @show settings
    sol_sto_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)
    oos_results         = compute_ex_ante_violprob(settings,networkdata,winddata,sol_sto_SOCP)
    push!(cost_violprob, [sol_sto_SOCP[:cost], oos_results[:ε_stat]])
end
p1 = plot((1 .- ε_vals),
         cost_violprob[!,:SOCP]./1e3,
         lw=4,
         xlim=[0.6,1.0],
         xlabel="Desired confidence level, (1 - ε)",
         ylabel="Expected cost in k\$",
         label="cost",
         markershape=:circle,
         colour=RGB(0,121/255,107/255))
p1 = plot!(twinx(),
          (1 .- ε_vals),
          xlim=[0.6,1.0],
          [cost_violprob[!,:viol_prob]],
          ylabel="Violation probability",
          label="viol. prob.",
          colour=RGBA(0,121/255,107/255,0.5),
          box=:on, grid=:off)
          
p1 = plot((1 .- ε_vals), cost_violprob[!,:SOCP]./1e3, title="Cost vs. Risk trade-off (σ = $(settings[:σ]))", label="", lw=4, xlabel="Desired confidence level, (1 - ε)", ylabel="Expected cost in k\$")
p2 = plot((1 .- ε_vals), cost_violprob[!,:viol_prob], title="Ex-ante constraint violation probability (σ = $(settings[:σ]))", label="", lw=4, xlabel="Desired confidence level, (1 - ε)", ylabel="Violation probability")
plot(p1,p2, layout = (2,1))
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "fig_cost_violation_prob_sigma.pdf"))


## Expt 3: SOCP vs LP : Computation time
num_scens = [10 50 100 200]
comp_cost = DataFrame(SOCPcost = Float64[], SOCPtime = Float64[], LPR1cost = Float64[], LPR1time = Float64[], LPR2cost=Float64[], LPR2time=Float64[])
for i in 1:length(num_scens)
    settings[:S] = num_scens[i]
    winddata, loaddata = prep_wind_load_data(settings,networkdata)
    sol_sto_SOCP        = sto_socp_mc(settings,networkdata,winddata,loaddata)

    settings[:bottleneck] == true ? settings[:MRR] = 0.5 : settings[:MRR] = 0.15
    sol_det_LP_DA       = det_lp_da_mc(settings,networkdata,winddata,loaddata)
    res_realtime_det_LP = DataFrame(iter = Int[], status = Any[], solvetime=Float64[], cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
    for scen = 1:winddata[:S]
        sol_det_LP_RT    = det_lp_rt_mc(settings,networkdata,winddata,loaddata,sol_det_LP_DA,scen)
        push!(res_realtime_det_LP, [scen, sol_det_LP_RT[:status], sol_det_LP_RT[:solvetime], round(sol_det_LP_RT[:cost],digits=2), round(sum(sol_det_LP_RT[:p_shed]),digits=2), round(sum(sol_det_LP_RT[:w_spill]), digits=2)])
    end

    sol_sto_LP          = sto_lp_mc(settings,networkdata,winddata,loaddata)
    push!(comp_cost, [sol_sto_SOCP[:cost]./1e3,
                      sol_sto_SOCP[:solvetime],
                      (sol_det_LP_DA[:cost] + mean(res_realtime_det_LP[!, :cost_act]))/1e3,
                      sol_det_LP_DA[:solvetime],
                      sol_sto_LP[:cost]./1e3,
                      sol_sto_LP[:solvetime]])
end
p1 = groupedbar(["10","10","10","50","50","50","100","100","100","200","200","200"],
                [comp_cost[!,:SOCPtime] comp_cost[!, :LPR1time] comp_cost[!, :LPR2time]],
                label=["Mcc" "R1" "R2"],
                colour= [RGB(0,121/255,107/255) RGB(1,87/255,34/255) RGB(199/255,211/255,215/255)],
                xlabel = "# of scenarios",
                ylabel="Solver time (s)",
                title="Computation time")
p2 = plot(["10","50","100","200"],
        [comp_cost[!,:SOCPcost] comp_cost[!, :LPR1cost] comp_cost[!, :LPR2cost]],
        label=["Mcc" "R1" "R2"],
        colour= [RGB(0,121/255,107/255) RGB(1,87/255,34/255) RGB(199/255,211/255,215/255)],
        title="Exp. cost",
        xlabel="# of scenarios", ylabel = "Expected cost (k\$)",
        lw=4)
plot(p1,p2, layout = (1,2))
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "cost_computation_vs_scenarios.csv"), comp_cost)
savefig(joinpath(RESULTS_PLOTS_DIR_NAME, "fig_cost_computation_vs_scenarios.pdf"))



## Out of Sample Simulations ##
# Augment settings with OOS data
settings[:TDsize]     = 500
# settings[:σ_oos]      = 0.1
# rt_premium            = 1.1 #real-time premium for adjustments
# oosdata = prep_outofsample_data(settings,winddata)
# #Compute SOCP model violation probability for oos
# oos_res = compute_ex_post_violprob(settings,networkdata,winddata,sol_sto_SOCP,oosdata)

#oos_res_df = DataFrame(model=Any[], networkconfig=Any[], rt_premium = Float64[], scen = Int[], oos_cost=Float64[])

settings[:bottleneck] == true ? netconfig = "congested" : netconfig = "acongested"
sol_realtime_SOCP_res = DataFrame(iter = Int[], status = Any[], comp_time=Float64[], cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
oos_infeas    = 0
cost_rec_flag_SOCP = 0
for scen = 1:oosdata[:S]
    sol_realtime_SOCP       = run_sto_SOCP_oos_reopt(settings,networkdata,winddata,sol_sto_SOCP,oosdata,scen)
    cost_recovery_check     = check_cost_rec_scen(settings,networkdata,sol_sto_SOCP,sol_realtime_SOCP)
    if(length(findall(x->x<-0.1, cost_recovery_check[:Pro_gen])) > 0)
        # @warn("Cost recovery failed for scenario -->$(scen)")
        cost_rec_flag_SOCP += 1
    end
    if (sol_realtime_SOCP[:status] ==  MOI.OPTIMAL)
        push!(sol_realtime_SOCP_res, [scen, sol_realtime_SOCP[:status], sol_realtime_SOCP[:solvetime], round(sol_realtime_SOCP[:cost],digits=4), round(sum(sol_realtime_SOCP[:p_shed]),digits=4), round(sum(sol_realtime_SOCP[:w_spill]), digits=4)])
        push!(oos_res_df,["Mcc",netconfig, 0, scen, round(sol_realtime_SOCP[:cost],digits=4)])
    else
        oos_infeas +=1
    end
end
@show sol_realtime_SOCP_res
println("OOS Infeasible ----> $(oos_infeas)")
println("Mean OOS Cost SOCP => ", mean(sol_realtime_SOCP_res[!, :cost_act])/1e3)
println("Increase over in-sample cost SOCP => $(100*(mean(sol_realtime_SOCP_res[!, :cost_act]) - sol_sto_SOCP[:cost])/sol_sto_SOCP[:cost])%")
println("Probability of load shed => $(100*(sum(sol_realtime_SOCP_res[!, :p_shed] .> 0)/settings[:TDsize]))%, Average load shed => $(100*sum(sol_realtime_SOCP_res[!, :p_shed])/(500*sum(loaddata[:D])))%")
println("Probability of wind spilled => $(100*(sum(sol_realtime_SOCP_res[!, :w_spill] .> 0)/settings[:TDsize]))%, Average wind spilled => $(100*sum(sol_realtime_SOCP_res[!, :w_spill])/(500*sum(winddata[:ŵ])))%")

#-- OOS Sample for R1 ---#
sol_realtime_R1_LP_scen_res = DataFrame(iter = Int[], status=Any[], comp_time = Float64[],  cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
oos_infeas = 0
cost_rec_flag_LP = 0
for scen in 1:oosdata[:S]
    sol_realtime_R1_LP_scen          = run_det_LP_oos_reopt(settings,networkdata,winddata,loaddata,sol_det_LP_DA,oosdata,scen,rt_premium)
    cost_recovery_check           = check_cost_rec_scen(settings,networkdata,sol_sto_LP,sol_realtime_R1_LP_scen)
    if(length(findall(x->x<-0.1, cost_recovery_check[:Pro_gen])) > 0)
         @warn("Cost recovery failed for scenario -->$(scen)")
         cost_rec_flag_LP +=1
    end
    if (sol_realtime_R1_LP_scen[:status] ==  MOI.OPTIMAL)
        push!(sol_realtime_R1_LP_scen_res, [scen, sol_realtime_R1_LP_scen[:status], sol_realtime_R1_LP_scen[:solvetime], round(sol_realtime_R1_LP_scen[:cost],digits=4) + sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr], round(sum(sol_realtime_R1_LP_scen[:p_shed]),digits=4), round(sum(sol_realtime_R1_LP_scen[:w_spill]),digits=4)])
        push!(oos_res_df,["R1", netconfig, rt_premium, scen, round(sol_realtime_R1_LP_scen[:cost],digits=4) + sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr]])
    else
        oos_infeas +=1
    end
end
@show sol_realtime_R1_LP_scen_res
println("OOS Infeasible ----> $(oos_infeas)")
println("Mean OOS Cost LP => ", mean(sol_realtime_R1_LP_scen_res[!, :cost_act])/1e3)
insample_expected_cost  = mean(res_realtime_det_LP[!, :cost_act]) + sol_det_LP_DA[:cost_res_gen] + sol_det_LP_DA[:cost_res_esr]
outsample_expected_cost = mean(sol_realtime_R1_LP_scen_res[!, :cost_act])
println("Increase over in-sample cost LP => $(100*(outsample_expected_cost-insample_expected_cost)/insample_expected_cost)%")
println("Probability of load shed => $(100*(sum(sol_realtime_R1_LP_scen_res[!, :p_shed] .> 0)/settings[:TDsize]))%, Average load shed => $(100*sum(sol_realtime_R1_LP_scen_res[!, :p_shed])/(500*sum(loaddata[:D])))%")
println("Probability of wind spilled => $(100*(sum(sol_realtime_R1_LP_scen_res[!, :w_spill] .> 0)/settings[:TDsize]))%, Average wind spilled => $(100*sum(sol_realtime_R1_LP_scen_res[!, :w_spill])/(500*sum(winddata[:ŵ])))%")

#-- OOS Sample for R2 --#
sol_realtime_R2_LP_scen_res = DataFrame(iter = Int[], status=Any[], comp_time = Float64[],  cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
oos_infeas = 0
cost_rec_flag_LP = 0
for scen in 1:oosdata[:S]
    sol_realtime_LP_scen          = run_sto_LP_oos_reopt(settings,networkdata,winddata,loaddata,sol_sto_LP,oosdata,scen,rt_premium)
    cost_recovery_check           = check_cost_rec_scen(settings,networkdata,sol_sto_LP,sol_realtime_LP_scen)
    if(length(findall(x->x<-0.1, cost_recovery_check[:Pro_gen])) > 0)
         @warn("Cost recovery failed for scenario -->$(scen)")
         cost_rec_flag_LP +=1
    end
    if (sol_realtime_LP_scen[:status] ==  MOI.OPTIMAL)
        push!(sol_realtime_R2_LP_scen_res, [scen, sol_realtime_LP_scen[:status], sol_realtime_LP_scen[:solvetime], round(sol_realtime_LP_scen[:cost],digits=4), round(sum(sol_realtime_LP_scen[:p_shed]),digits=4), round(sum(sol_realtime_LP_scen[:w_spill]),digits=4)])
        push!(oos_res_df,["R2", netconfig, rt_premium, scen, round(sol_realtime_LP_scen[:cost],digits=4)])
    else
        oos_infeas +=1
    end
end
@show sol_realtime_R2_LP_scen_res
println("OOS Infeasible ----> $(oos_infeas)")
println("Mean OOS Cost LP => ", mean(sol_realtime_R2_LP_scen_res[!, :cost_act])/1e3)
println("Increase over in-sample cost LP => $(100*(mean(sol_realtime_R2_LP_scen_res[!, :cost_act]) - sol_sto_LP[:cost])/sol_sto_LP[:cost])%")
println("Probability of load shed => $(100*(sum(sol_realtime_R2_LP_scen_res[!, :p_shed] .> 0)/settings[:TDsize]))%, Average load shed => $(100*sum(sol_realtime_R2_LP_scen_res[!, :p_shed])/(500*sum(loaddata[:D])))%")
println("Probability of wind spilled => $(100*(sum(sol_realtime_R2_LP_scen_res[!, :w_spill] .> 0)/settings[:TDsize]))%, Average wind spilled => $(100*sum(sol_realtime_R2_LP_scen_res[!, :w_spill])/(500*sum(winddata[:ŵ])))%")

@show oos_res_df
CSV.write(joinpath(RESULTS_DATA_DIR_NAME, "new_oos_cost_comparison_500scens.csv"), oos_res_df)
