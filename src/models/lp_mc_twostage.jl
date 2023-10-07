function det_lp_da_mc(settings,networkdata,winddata,loaddata)
    """ Deterministic Dayahead Two-stage Market Clearing, Mdet_DA """
    #Agents
    Ng = length(networkdata[:gens])
    Nw = size(winddata[:ŵ],1)
    Ns = length(networkdata[:esrs])

    m = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0, "QUIET" => true))
    #m = Model(optimizer_with_attributes(Ipopt.Optimizer))

    # ---- DEFINE Primal VARIABLES ------ #
    #For power and participation factor
    @variable(m, p[1:Ng,    1:settings[:Nsteps], 1:settings[:T]])       #nominal power
    @variable(m, r_p[1:Ng,  1:settings[:T]])                            #reserve capacity
    @variable(m, b[1:Ns,    1:settings[:T]])       #nominal storage operation
    @variable(m, r_b[1:Ns,  1:settings[:T]])                            #reserve capacity storage

    @objective(m, Min, sum(sum(sum(pwl_cost_approx(settings,networkdata,g,0)[:slopes] .* p[g,:,t]) for g in 1:length(networkdata[:gens])) for t in 1:settings[:T])
                       + sum((ones(settings[:T]) * networkdata[:esrs][s].c_l)'b[s,:] for s in 1:length(networkdata[:esrs]))
                       + sum((ones(settings[:T])*networkdata[:gens][g].c_r)'r_p[g,:] for g in 1:length(networkdata[:gens]))
                       + sum((ones(settings[:T])*networkdata[:esrs][s].c_r)'r_b[s,:] for s in 1:length(networkdata[:esrs])))

    #--Generation Constraints--#
    #Generator cost segments limits
    @constraint(m, p_seg_lim[g=1:Ng, seg=1:settings[:Nsteps], t=1:settings[:T]], 0 <= p[g,seg,t] <= (networkdata[:gens][g].p̅ - networkdata[:gens][g].p̲)/settings[:Nsteps])
    #Generator production limits
    @constraint(m, p_gen_lim_mx[g=1:Ng, t=1:settings[:T]], sum(p[g,seg,t] for seg in 1:settings[:Nsteps]) >= networkdata[:gens][g].p̲)
    @constraint(m, p_gen_lim_mn[g=1:Ng, t=1:settings[:T]], sum(p[g,seg,t] for seg in 1:settings[:Nsteps]) <= networkdata[:gens][g].p̅)
    #Generator production and reserve limits
    @constraint(m, p_res_lim_mx[g=1:Ng, t=1:settings[:T]], sum(p[g,seg,t] for seg in 1:settings[:Nsteps]) + r_p[g,t] <= networkdata[:gens][g].p̅)
    @constraint(m, p_res_lim_mn[g=1:Ng, t=1:settings[:T]], sum(p[g,seg,t] for seg in 1:settings[:Nsteps]) - r_p[g,t] >= networkdata[:gens][g].p̲)
    #Generator reserve limits
    @constraint(m, rese_lim_min[g=1:Ng, t=1:settings[:T]], 0 <= r_p[g,t])
    @constraint(m, rese_lim_max[g=1:Ng, t=1:settings[:T]], r_p[g,t] <= networkdata[:gens][g].r̅)
    #Generator ramping limits
    @constraint(m, p_ram_lim_max[g=1:Ng, t=2:settings[:T]], (sum(p[g,seg,t] for seg in 1:settings[:Nsteps]) + r_p[g,t]) - (sum(p[g,seg,t-1] for seg in 1:settings[:Nsteps]) + r_p[g,t-1]) <=    networkdata[:gens][g].Ra̅)
    @constraint(m, p_ram_lim_min[g=1:Ng, t=2:settings[:T]], (sum(p[g,seg,t] for seg in 1:settings[:Nsteps]) + r_p[g,t]) - (sum(p[g,seg,t-1] for seg in 1:settings[:Nsteps]) + r_p[g,t-1]) >=  - networkdata[:gens][g].Ra̅)

    #--Energy Storage Resource Constraints--#
    #Storage charging and discharging limits
    @constraint(m, esr_pd_lim[s=1:Ns,  t=1:settings[:T]], b[s,t] <= networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d)
    @constraint(m, esr_pc_lim[s=1:Ns,  t=1:settings[:T]], b[s,t] >= -round(1/networkdata[:esrs][s].η_c,digits=4)*networkdata[:esrs][s].P̅c)
    #Storage production and reserve limits
    @constraint(m, esr_chg_lim[s=1:Ns, t=1:settings[:T]], b[s,t] + r_b[s,t] <= networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d)
    @constraint(m, esr_dis_lim[s=1:Ns, t=1:settings[:T]], b[s,t] - r_b[s,t] >= - round(1/networkdata[:esrs][s].η_c,digits=4)*networkdata[:esrs][s].P̅c)
    #Storage reserve limit - non-negative
    @constraint(m, esr_lim_min[s=1:Ns, t=1:settings[:T]], 0 <= r_b[s,t])

    #Storage : SOC minimum and maximum energy limits
    @constraint(m, esr_spill_up[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(b[s,t̂] + r_b[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].s̅)
    @constraint(m, esr_spill_dn[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(b[s,t̂] + r_b[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].s̲)
    #Storage : end-of-day energy neutrality limits
    @constraint(m, esr_eod_max[s=1:Ns,  t=settings[:T]],   networkdata[:esrs][s].E𝟶 - sum(b[s,t̂] + r_b[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].E𝟶 + networkdata[:esrs][s].B_s)
    @constraint(m, esr_eod_min[s=1:Ns,  t=settings[:T]],   networkdata[:esrs][s].E𝟶 - sum(b[s,t̂] + r_b[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].B_s)

    #Power flow in the lines - PTDF formulation
    @constraint(m, plims_max[l=1:length(networkdata[:lines]),t=1:settings[:T]], -(networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
                                                                                +(networkdata[:Ψ]networkdata[:Cgens]*(sum(p[:,seg,t] for seg in 1:settings[:Nsteps])))[l]
                                                                                +(networkdata[:Ψ]networkdata[:Cesr] *(b[:,t]))[l]
                                                                                +(networkdata[:Ψ]winddata[:Cwind]winddata[:ŵ][:,t])[l] <=  networkdata[:f̅][l])
    @constraint(m, plims_min[l=1:length(networkdata[:lines]),t=1:settings[:T]], -(networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
                                                                                +(networkdata[:Ψ]networkdata[:Cgens]*(sum(p[:,seg,t] for seg in 1:settings[:Nsteps])))[l]
                                                                                +(networkdata[:Ψ]networkdata[:Cesr] *(b[:,t]))[l]
                                                                                +(networkdata[:Ψ]winddata[:Cwind]winddata[:ŵ][:,t])[l] >= -networkdata[:f̅][l])

    # balance constraints
    @constraint(m, λ_e[t=1:settings[:T]], sum(sum(p[g,seg,t] for seg in 1:settings[:Nsteps]) for g=1:Ng)
                                        + sum(winddata[:ŵ][k,t] for k=1:Nw)
                                        + sum(b[s,t] for s=1:Ns)
                                        - loaddata[:D][t] == 0)
    @constraint(m, λ_r[t=1:settings[:T]], sum(r_p[g,t] for g=1:Ng)
                                        + sum(r_b[s,t] for s=1:Ns) >= settings[:MRR]*sum(winddata[:ŵ][k,t] for k=1:Nw))
    # solve
    @time optimize!(m)
    status = termination_status(m)
    @info("Primal Spatial LP status ---> $(status)")

    nomflows = zeros(length(networkdata[:lines]), settings[:T])
    for t in 1:settings[:T]
        nomflows[:,t] = (networkdata[:Ψ]*(loaddata[:node_loads][:,t]
                       - networkdata[:Cgens]*JuMP.value.(sum(p[:,seg,t] for seg in 1:settings[:Nsteps]))
                       - networkdata[:Cesr] *JuMP.value.(b[:,t])
                       - winddata[:Cwind]*winddata[:ŵ][:,t]))
    end

    # return solution
    solution = Dict(
    :p              => JuMP.value.(p),
    :r_p            => JuMP.value.(r_p),
    :b              => JuMP.value.(b),
    :r_b            => JuMP.value.(r_b),
    :λ_E            => JuMP.dual.(λ_e),
    :λ_R            => JuMP.dual.(λ_r),
    :ŵ              => winddata[:ŵ],
    :nomflows       => nomflows,
    :cost           => JuMP.objective_value.(m),
    :cost_res_gen   => sum((ones(settings[:T])*networkdata[:gens][g].c_r)' * JuMP.value.(r_p[g,:]) for g in 1:length(networkdata[:gens])),
    :cost_res_esr   => sum((ones(settings[:T])*networkdata[:esrs][s].c_r)' * JuMP.value.(r_b[s,:]) for s in 1:length(networkdata[:esrs])),
    :solvetime      => MOI.get(m, MOI.SolveTime()),
    :model          => m)
    return solution
end

function det_lp_rt_mc(settings,networkdata,winddata,loaddata,sol_det_LP_DA,scen)
    """ Deterministic Reserve Activation, Mdet_RT """
    #Agents
    Ng = length(networkdata[:gens])
    Nw = size(winddata[:ŵ],1)
    Ns = length(networkdata[:esrs])
    C_shed  = 500
    C_spill = 500
    Δ = sum(winddata[:ξ][w,scen,:] for w = 1:Nw)

    m = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0, "QUIET" => true))
    # ---- DEFINE Primal VARIABLES ------ #
    #For power and participation factor
    @variable(m, r_p[1:Ng, 1:settings[:Nsteps], 1:settings[:T]]) #generator reserves activated
    @variable(m, r_b[1:Ns, 1:settings[:T]])                      #storage reserves activated
    @variable(m, p_shed[1:length(networkdata[:buses]), 1:settings[:T]]) #load_shedding to ensure feasibility
    @variable(m, w_spill[1:Nw,1:settings[:T]])      #windspillage to ensure feasibility

    @objective(m, Min,  sum(sum(sum(pwl_cost_approx(settings,networkdata,g,0)[:slopes] .* (sol_det_LP_DA[:p][g,:,t] .+ r_p[g,:,t])) for g in 1:length(networkdata[:gens])) for t in 1:settings[:T])
                       + sum((ones(settings[:T]) * networkdata[:esrs][s].c_l)'sol_det_LP_DA[:b][s,:] .+ (ones(settings[:T]) * networkdata[:esrs][s].c_l)'r_b[s,:] for s in 1:length(networkdata[:esrs]))
                       + sum((C_shed*ones(settings[:T]))'p_shed[d,:] for d in 1:length(networkdata[:buses]))
                       + sum((C_spill*ones(settings[:T]))'w_spill[w,:] for w in 1:Nw))

    #Reserve activation segments limits
    @constraint(m, r_p_lims[g=1:Ng, seg=1:settings[:Nsteps], t=1:settings[:T]], -sol_det_LP_DA[:r_p][g,t]/settings[:Nsteps] <= r_p[g,seg,t] <= sol_det_LP_DA[:r_p][g,t]/settings[:Nsteps])
    #Balance of reserve activation
    @constraint(m, r_p_act[t=1:settings[:T]], sum(sum(r_p[g,seg,t] for seg in 1:settings[:Nsteps]) for g in 1:Ng) + sum(r_b[s,t] for s in 1:Ns) + sum(p_shed[b,t] for b in 1:length(networkdata[:buses])) == Δ[t] + sum(w_spill[k,t] for k in 1:Nw))
    #Reserve activation limits within the day-ahead set values
    @constraint(m, p_res_act[g=1:Ng, t=1:settings[:T]], -sol_det_LP_DA[:r_p][g,t] <= sum(r_p[g,seg,t] for seg in 1:settings[:Nsteps]) <= sol_det_LP_DA[:r_p][g,t])
    @constraint(m, b_res_act[s=1:Ns, t=1:settings[:T]], -sol_det_LP_DA[:r_b][s,t] <= r_b[s,t] <= sol_det_LP_DA[:r_b][s,t])
    #load shedding and wind spillage limits
    @constraint(m, shed_lims[d=1:length(networkdata[:buses]), t=1:settings[:T]], 0 <=  p_shed[d,t]  <= loaddata[:node_loads][d,t])
    @constraint(m, spillage[w=1:Nw, t=1:settings[:T]],  0 <=  w_spill[w,t] <= winddata[:ŵ][w,t] - winddata[:ξ][:,scen,:][w,t])

    #Power flow in the lines - PTDF formulation
    @constraint(m, plims_max[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]*(loaddata[:node_loads][:,t] - p_shed[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cgens]*(sum(sol_det_LP_DA[:p][:,seg,t] for seg in 1:settings[:Nsteps]) + sum(r_p[:,seg,t] for seg in 1:settings[:Nsteps])))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cesr] *(sol_det_LP_DA[:b][:,t]+ r_b[:,t]))[l]
                                                                                 +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - winddata[:ξ][:,scen,:][:,t]))[l] <=  networkdata[:f̅][l])
    @constraint(m, plims_min[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]*(loaddata[:node_loads][:,t] - p_shed[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cgens]*(sum(sol_det_LP_DA[:p][:,seg,t] for seg in 1:settings[:Nsteps]) + sum(r_p[:,seg,t] for seg in 1:settings[:Nsteps])))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cesr] *(sol_det_LP_DA[:b][:,t] + r_b[:,t]))[l]
                                                                                 +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - winddata[:ξ][:,scen,:][:,t]))[l] >= -networkdata[:f̅][l])

    optimize!(m)


    solution = Dict(
    :cost    => JuMP.objective_value(m),
    :r_p        => JuMP.value.(r_p),
    :r_b        => JuMP.value.(r_b),
    :p_shed     => JuMP.value.(p_shed),
    :w_spill    => JuMP.value.(w_spill),
    :status     => termination_status(m),
    :solvetime  => MOI.get(m, MOI.SolveTime())
    )
    return solution
end

function run_det_LP_oos_reopt(settings,networkdata,winddata,loaddata,sol_det_LP_DA,oosdata,scen,rt_pen)
    """Scenario-based stochastic market clearing OOS simulations"""
    #Agents
    Ng = length(networkdata[:gens])
    Ns = length(networkdata[:esrs])
    Nw = size(winddata[:ŵ],1)
    C_shed  = 500
    C_spill = 500
    Δ = sum(oosdata[:ξ_oos][w,scen,:] for w in 1:Nw)

    m = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0, "QUIET" => true))
    # ---- Define primal variables ---- #
    # For adjustments to power
    @variable(m, r_p[1:Ng, 1:settings[:T]])    #real-time adjustment
    @variable(m, r_b[1:Ns, 1:settings[:T]])    #real-time adjustment
    @variable(m, p_shed[1:length(networkdata[:buses]), 1:settings[:T]]) #load_shedding to ensure feasibility
    @variable(m, w_spill[1:Nw,1:settings[:T]])      #windspillage to ensure feasibility

    @objective(m, Min,  sum((sum(sol_det_LP_DA[:p][g,seg,:] for seg in 1:settings[:Nsteps]))'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*(sum(sol_det_LP_DA[:p][g,seg,:] for seg in 1:settings[:Nsteps])) for g in 1:length(networkdata[:gens]))
                        + sum((r_p[g,:])'*diagm(ones(settings[:T])*networkdata[:gens][g].q*rt_pen)*(r_p[g,:]) for g in 1:length(networkdata[:gens]))
                        + sum((ones(settings[:T])*networkdata[:gens][g].c)' * (sum(sol_det_LP_DA[:p][g,seg,:] for seg in 1:settings[:Nsteps]))  for g in 1:length(networkdata[:gens]))
                        + sum((ones(settings[:T])*networkdata[:gens][g].c*rt_pen)' * (r_p[g,:])  for g in 1:length(networkdata[:gens]))
                        + sum((ones(settings[:T])*networkdata[:esrs][s].c_l)' * (sol_det_LP_DA[:b][s,:]) for s in 1:length(networkdata[:esrs]))
                        + sum((ones(settings[:T])*networkdata[:esrs][s].c_l*rt_pen)' * (r_b[s,:]) for s in 1:length(networkdata[:esrs]))
                        + sum((C_shed*ones(settings[:T]))'p_shed[d,:] for d in 1:length(networkdata[:buses]))
                        + sum((C_spill*ones(settings[:T]))'w_spill[w,:] for w in 1:Nw))

    #generator constraints rt
    @constraint(m, p_res_lim_mx[g=1:Ng, t=1:settings[:T]], sum(sol_det_LP_DA[:p][:,seg,:] for seg in 1:settings[:Nsteps])[g,t] + r_p[g,t] <= networkdata[:gens][g].p̅)
    @constraint(m, p_res_lim_mn[g=1:Ng, t=1:settings[:T]], sum(sol_det_LP_DA[:p][:,seg,:] for seg in 1:settings[:Nsteps])[g,t] + r_p[g,t] >= networkdata[:gens][g].p̲)
    @constraint(m, rese_lim_max[g=1:Ng, t=1:settings[:T]], -networkdata[:gens][g].r̅ <= r_p[g,t] <= networkdata[:gens][g].r̅)
    @constraint(m, p_res_ramp_a[g=1:Ng, t=2:settings[:T]], -networkdata[:gens][g].Ra̅ <= (sum(sol_det_LP_DA[:p][:,seg,:] for seg in 1:settings[:Nsteps])[g,t] + r_p[g,t]) - (sum(sol_det_LP_DA[:p][:,seg,:] for seg in 1:settings[:Nsteps])[g,t-1] + r_p[g,t-1]) <=  networkdata[:gens][g].Ra̅)

    #storage constraints
    @constraint(m, esr_chg_lim[s=1:Ns, t=1:settings[:T]],  sol_det_LP_DA[:b][s,t] + r_b[s,t] <= networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d)
    @constraint(m, esr_dis_lim[s=1:Ns, t=1:settings[:T]],  sol_det_LP_DA[:b][s,t] + r_b[s,t] >= - 1/networkdata[:esrs][s].η_c*networkdata[:esrs][s].P̅c)
    @constraint(m, esr_spill_up[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(sol_det_LP_DA[:b][s,t̂] + r_b[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].s̅)
    @constraint(m, esr_spill_dn[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(sol_det_LP_DA[:b][s,t̂] + r_b[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].s̲)
    @constraint(m, esr_eod_max[s=1:Ns, t=settings[:T]],    networkdata[:esrs][s].E𝟶 - sum(sol_det_LP_DA[:b][s,t̂] + r_b[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].E𝟶 + networkdata[:esrs][s].B_s)
    @constraint(m, esr_eod_min[s=1:Ns, t=settings[:T]],    networkdata[:esrs][s].E𝟶 - sum(sol_det_LP_DA[:b][s,t̂] + r_b[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].B_s)

    #load shedding and wind spillage limits
    @constraint(m, shed_lims[d=1:length(networkdata[:buses]), t=1:settings[:T]], 0 <=  p_shed[d,t]  <= loaddata[:node_loads][d,t])
   @constraint(m, spillage[w=1:Nw, t=1:settings[:T]],  0 <=  w_spill[w,t] <= winddata[:ŵ][w,t] - oosdata[:ξ_oos][:,scen,:][w,t])

    # #Power flow in the lines - PTDF formulation
    @constraint(m, plims_max[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]*(loaddata[:node_loads][:,t] - p_shed[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cgens]*(sum(sol_det_LP_DA[:p][:,seg,:] for seg in 1:settings[:Nsteps])[:,t] + r_p[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cesr]*(sol_det_LP_DA[:b][:,t] + r_b[:,t]))[l]
                                                                                 +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - oosdata[:ξ_oos][:,scen,:][:,t] - w_spill[:,t]))[l] <=  networkdata[:f̅][l])
    @constraint(m, plims_min[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]*(loaddata[:node_loads][:,t] - p_shed[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cgens]*(sum(sol_det_LP_DA[:p][:,seg,:] for seg in 1:settings[:Nsteps])[:,t] + r_p[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cesr]*(sol_det_LP_DA[:b][:,t] + r_b[:,t]))[l]
                                                                                 +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - oosdata[:ξ_oos][:,scen,:][:,t] - w_spill[:,t]))[l] >= -networkdata[:f̅][l])

    # # balance constraints
    @constraint(m, r_p_activation_balance[t=1:settings[:T]], sum(r_p[g,t] for g in 1:Ng)
                                                           + sum(r_b[s,t] for s in 1:Ns)
                                                           + sum(p_shed[b,t] for b in 1:length(networkdata[:buses]))
                                                           + sum((winddata[:ŵ][w,t] - oosdata[:ξ_oos][:,scen,:][w,t]) - sol_det_LP_DA[:ŵ][w,t] - w_spill[w,t] for w in 1:Nw)
                                                           == 0)

    # solve
    optimize!(m)

    solution = Dict(
    :cost       => JuMP.objective_value(m),
    :r_p        => JuMP.value.(r_p),
    :r_b        => JuMP.value.(r_b),
    :p_shed     => JuMP.value.(p_shed),
    :w_spill    => JuMP.value.(w_spill),
    :status     => termination_status(m),
    :solvetime => 0.0,
    :da_model   => "R1",
    )
    return solution
end



function determine_MRR(settings,networkdata,winddata,loaddata)
    """Endogenously determine MRR value for R1 based on LOLP"""
    @info("Finding optimal MRR")
    @info("---------------")
    MRR = settings[:MRR]
    @info("MRR --> $(MRR)")
    sol_det_LP_DA = det_lp_da_mc(settings,networkdata,winddata,loaddata)
    res_realtime_det_LP = DataFrame(iter = Int[], status = Any[], solvetime=Float64[], cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
    for scen = 1:winddata[:S]
        sol_det_LP_RT    = det_lp_rt_mc(settings,networkdata,winddata,loaddata,sol_det_LP_DA,scen)
        push!(res_realtime_det_LP, [scen, sol_det_LP_RT[:status], sol_det_LP_RT[:solvetime], round(sol_det_LP_RT[:cost],digits=2), round(sum(sol_det_LP_RT[:p_shed]),digits=2), round(sum(sol_det_LP_RT[:w_spill]), digits=2)])
    end
    @info("P Shed --> $(mean(res_realtime_det_LP[!,:p_shed]))")
    while ((count(x->x>0, res_realtime_det_LP[!, :p_shed])/winddata[:S] >= 0.05))
        MRR = MRR + 0.01
        settings[:MRR] = MRR
        sol_det_LP_DA = det_lp_da_mc(settings,networkdata,winddata,loaddata)
        res_realtime_det_LP = DataFrame(iter = Int[], status = Any[], solvetime=Float64[], cost_act = Float64[], p_shed = Float64[], w_spill = Float64[])
        for scen = 1:winddata[:S]
            sol_det_LP_RT    = det_lp_rt_mc(settings,networkdata,winddata,loaddata,sol_det_LP_DA,scen)
            push!(res_realtime_det_LP, [scen, sol_det_LP_RT[:status], sol_det_LP_RT[:solvetime], round(sol_det_LP_RT[:cost],digits=2), round(sum(sol_det_LP_RT[:p_shed]),digits=2), round(sum(sol_det_LP_RT[:w_spill]), digits=2)])
        end
        @show res_realtime_det_LP
        @info("MRR --> $(MRR)")
        @info("P Shed --> $(mean(res_realtime_det_LP[!,:p_shed]))")
    end
    @info("Optimal MRR --> $(MRR)")
    @info("---------------")
    return MRR
end
