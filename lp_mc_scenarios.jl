function sto_lp_mc(settings,networkdata,winddata,loaddata)
    """Scenario-based Stochastic Market Clearing with LP"""
    #Agents
    Ng = length(networkdata[:gens])
    Nw = size(winddata[:ŵ],1)
    Ns = length(networkdata[:esrs])
    C_shed  = 3000
    m = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0, "QUIET" => true))
    #DA variables
    @variable(m, p_da[1:Ng, 1:settings[:Nsteps], 1:settings[:T]])       #DA power production from generators
    @variable(m, b_da[1:Ns, 1:settings[:T]])       #DA power production from ESRs
    @variable(m, w_da[1:Nw, 1:settings[:T]])                            #DA wind power
    #RT variables
    @variable(m, p_rt[1:Ng, 1:settings[:Nsteps], 1:settings[:T], 1:winddata[:S]])   #RT power production
    @variable(m, b_rt[1:Ns, 1:settings[:T], 1:winddata[:S]])   #RT production from ESRs
    @variable(m, w_rt[1:Nw, 1:settings[:T], 1:winddata[:S]])                        #Wind RT production
    @variable(m, p_shed[1:length(networkdata[:buses]), 1:settings[:T], 1:winddata[:S]]) #Load shed

    @objective(m, Min, sum(sum(sum(pwl_cost_approx(settings,networkdata,g,0)[:slopes] .* p_da[g,:,t]) for g in 1:length(networkdata[:gens])) for t in 1:settings[:T])
                      + sum((ones(settings[:T]) * networkdata[:esrs][s].c_l)'b_da[s,:] for s in 1:length(networkdata[:esrs]))
                      + sum((1/winddata[:S])*(sum(sum(sum(pwl_cost_approx(settings,networkdata,g,0)[:slopes] .* p_rt[g,:,t,scen]) for g in 1:length(networkdata[:gens])) for t in 1:settings[:T])
                      + sum((ones(settings[:T]) * networkdata[:esrs][s].c_l)'b_rt[s,:,scen] for s in 1:length(networkdata[:esrs]))
                                + sum((C_shed*ones(settings[:T]))'p_shed[b,:,scen] for b in 1:length(networkdata[:buses]))) for scen in 1:winddata[:S]))
    ##----Generator constraints - DA---##
    #Limits on cost segments
    @constraint(m, p_seg_lim[g=1:Ng, seg=1:settings[:Nsteps], t=1:settings[:T]], 0 <= p_da[g,seg,t] <= (networkdata[:gens][g].p̅ - networkdata[:gens][g].p̲)/settings[:Nsteps])
    #Limits on production max and min
    @constraint(m, p_da_lim_mx[g=1:Ng, t=1:settings[:T]], sum(p_da[g,seg,t] for seg in 1:settings[:Nsteps]) >= networkdata[:gens][g].p̲)
    @constraint(m, p_da_lim_mn[g=1:Ng, t=1:settings[:T]], sum(p_da[g,seg,t] for seg in 1:settings[:Nsteps]) <= networkdata[:gens][g].p̅)
    #limits on ramping - max and min
    @constraint(m, p_da_ramp_lim_max[g=1:Ng, t=2:settings[:T]], sum(p_da[g,seg,t] for seg in 1:settings[:Nsteps]) - sum(p_da[g,seg,t-1] for seg in 1:settings[:Nsteps]) <= networkdata[:gens][g].Ra̅)
    @constraint(m, p_da_ramp_lim_min[g=1:Ng, t=2:settings[:T]], sum(p_da[g,seg,t] for seg in 1:settings[:Nsteps]) - sum(p_da[g,seg,t-1] for seg in 1:settings[:Nsteps]) >= -networkdata[:gens][g].Ra̅)
    ##----ESR constraints - DA---##
    #limits on charging and discharging
    @constraint(m, b_da_lim_ch[s=1:Ns, t=1:settings[:T]], b_da[s,t] <= networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d)
    @constraint(m, b_da_lim_ds[s=1:Ns, t=1:settings[:T]], b_da[s,t] >= -1/networkdata[:esrs][s].η_c*networkdata[:esrs][s].P̅c)
    #ESR SOC limits minimum and maximum constraints
    @constraint(m, esr_spill_up_da[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(b_da[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].s̅)
    @constraint(m, esr_spill_dn_da[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(b_da[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].s̲)
    #ESR EOD content limits
    @constraint(m, esr_eod_max_da[s=1:Ns,  t = settings[:T]], networkdata[:esrs][s].E𝟶 - sum(b_da[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].E𝟶 + networkdata[:esrs][s].B_s)
    @constraint(m, esr_eod_min_da[s=1:Ns,  t = settings[:T]], networkdata[:esrs][s].E𝟶 - sum(b_da[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].B_s)
    ##----WPP constraints - DA---##
    @constraint(m, w_da_lim_mx[k=1:Nw, t=1:settings[:T]], w_da[k,t] <= winddata[:Wcap])
    @constraint(m, w_da_lim_mn[k=1:Nw, t=1:settings[:T]], w_da[k,t] >= 0)
    ##----System Constraints - DA---##
    #DA Power Balance
    @constraint(m, λ_da_sys[t=1:settings[:T]], sum(sum(p_da[g,seg,t] for seg in 1:settings[:Nsteps]) for g=1:Ng) + sum(w_da[k,t] for k=1:Nw) + sum(b_da[s,t] for s=1:Ns) - loaddata[:D][t] == 0)
    #DA Flow limits
    @constraint(m, ρ̅_da[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
                                                                            +(networkdata[:Ψ]networkdata[:Cgens]*(sum(p_da[:,seg,t] for seg in 1:settings[:Nsteps])))[l]
                                                                            +(networkdata[:Ψ]networkdata[:Cesr]*(b_da[:,t]))[l]
                                                                            +(networkdata[:Ψ]winddata[:Cwind]*(w_da[:,t]))[l] <=  networkdata[:f̅][l])
    @constraint(m, ρ̲_da[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
                                                                            +(networkdata[:Ψ]networkdata[:Cgens]*(sum(p_da[:,seg,t] for seg in 1:settings[:Nsteps])))[l]
                                                                            +(networkdata[:Ψ]networkdata[:Cesr]*(b_da[:,t]))[l]
                                                                            +(networkdata[:Ψ]winddata[:Cwind]*(w_da[:,t]))[l] >= -networkdata[:f̅][l])

    ##----Generator Constraints - Scenario-wise---##
    #Generator cost RT segments limits
    @constraint(m, p_rt_sg_lim[g=1:Ng, seg=1:settings[:Nsteps], t=1:settings[:T], scen=1:winddata[:S]], 0 <= p_rt[g,seg,t,scen] <= (networkdata[:gens][g].r̅ + networkdata[:gens][g].r̅)/settings[:Nsteps])
    #Generator adjustment limits RT
    @constraint(m, p_rt_lim_mx[g=1:Ng, t=1:settings[:T], scen=1:winddata[:S]], sum(p_rt[g,seg,t,scen] for seg in 1:settings[:Nsteps]) <=  networkdata[:gens][g].r̅)
    @constraint(m, p_rt_lim_mn[g=1:Ng, t=1:settings[:T], scen=1:winddata[:S]], sum(p_rt[g,seg,t,scen] for seg in 1:settings[:Nsteps]) >= -networkdata[:gens][g].r̅)
    #Generator total production DA + RT
    @constraint(m, p_tt_lim_mx[g=1:Ng, t=1:settings[:T], scen=1:winddata[:S]], sum(p_da[g,seg,t] for seg in 1:settings[:Nsteps]) + sum(p_rt[g,seg,t,scen] for seg in 1:settings[:Nsteps]) <= networkdata[:gens][g].p̅)
    @constraint(m, p_tt_lim_mn[g=1:Ng, t=1:settings[:T], scen=1:winddata[:S]], sum(p_da[g,seg,t] for seg in 1:settings[:Nsteps]) + sum(p_rt[g,seg,t,scen] for seg in 1:settings[:Nsteps]) >= networkdata[:gens][g].p̲)
    #Generator ramping constraints in RT
    @constraint(m, p_tt_ram_up[g=1:Ng, t=2:settings[:T], scen=1:winddata[:S]],
                                                    (sum(p_da[g,seg,t]   + p_rt[g,seg,t,scen]    for seg in 1:settings[:Nsteps]))
                                                  - (sum(p_da[g,seg,t-1] + p_rt[g,seg,t-1,scen]  for seg in 1:settings[:Nsteps]))
                                                    <=  networkdata[:gens][g].Ra̅)
    @constraint(m, p_tt_ram_dn[g=1:Ng, t=2:settings[:T], scen=1:winddata[:S]],
                                                    (sum(p_da[g,seg,t]   + p_rt[g,seg,t,scen]    for seg in 1:settings[:Nsteps]))
                                                  - (sum(p_da[g,seg,t-1] + p_rt[g,seg,t-1,scen]  for seg in 1:settings[:Nsteps]))
                                                    >= -networkdata[:gens][g].Ra̅)
    ##----ESR Constraints - Scenario-wise---##
    #ESR total production limits DA + RT
    @constraint(m, b_tt_lim_ch[s=1:Ns, t=1:settings[:T], scen=1:winddata[:S]], b_da[s,t] + b_rt[s,t,scen] <= networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d)
    @constraint(m, b_tt_lim_ds[s=1:Ns, t=1:settings[:T], scen=1:winddata[:S]], b_da[s,t] + b_rt[s,t,scen] >= -1/networkdata[:esrs][s].η_c*networkdata[:esrs][s].P̅c)
    # ESR SOC limits minimum and maximum constraints DA + RT
    @constraint(m, esr_spill_up_rt[s=1:Ns, t=1:settings[:T], scen=1:winddata[:S]], networkdata[:esrs][s].E𝟶 - sum(b_da[s,t̂] + b_rt[s,t̂,scen] for t̂ in 1:t) <= networkdata[:esrs][s].s̅)
    @constraint(m, esr_spill_dn_rt[s=1:Ns, t=1:settings[:T], scen=1:winddata[:S]], networkdata[:esrs][s].E𝟶 - sum(b_da[s,t̂] + b_rt[s,t̂,scen] for t̂ in 1:t) >= networkdata[:esrs][s].s̲)
    #ESR EDO min and max DA + RT
    @constraint(m, esr_eod_max_rt[s=1:Ns, t=settings[:T],    scen=1:winddata[:S]], networkdata[:esrs][s].E𝟶 - sum(b_da[s,t̂] + b_rt[s,t̂,scen] for t̂ in 1:t) <= networkdata[:esrs][s].E𝟶 + networkdata[:esrs][s].B_s)
    @constraint(m, esr_eod_min_rt[s=1:Ns, t=settings[:T],    scen=1:winddata[:S]], networkdata[:esrs][s].E𝟶 - sum(b_da[s,t̂] + b_rt[s,t̂,scen] for t̂ in 1:t) >= networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].B_s)

    #wind constraints - RT
    @constraint(m, w_rt_mx[k=1:Nw, t=1:settings[:T], scen=1:winddata[:S]], w_da[k,t] + w_rt[k,t,scen] <= winddata[:ŵ][k,t] - winddata[:ξ][k,scen,t])
    @constraint(m, w_rt_mn[k=1:Nw, t=1:settings[:T], scen=1:winddata[:S]], w_da[k,t] + w_rt[k,t,scen] >= 0)
    #load shedding - RT
    @constraint(m, pshed_rt_mx[b=1:length(networkdata[:buses]), t=1:settings[:T], scen=1:winddata[:S]], p_shed[b,t,scen] <= loaddata[:node_loads][b,t])
    @constraint(m, pshed_rt_mn[b=1:length(networkdata[:buses]), t=1:settings[:T], scen=1:winddata[:S]], p_shed[b,t,scen] >= 0)
    # balance constraint - RT
    @constraint(m, λ_rt_sys[t=1:settings[:T], scen=1:winddata[:S]], sum(sum(p_rt[g,seg,t,scen] for seg in 1:settings[:Nsteps]) for g=1:Ng) + sum(b_rt[s,t,scen] for s=1:Ns) + sum(w_rt[k,t,scen] for k=1:Nw) + sum(p_shed[b,t,scen] for b=1:length(networkdata[:buses])) == 0)
    #Power flows real-time
    @constraint(m, ρ̅_rt[l=1:length(networkdata[:lines]), t=1:settings[:T], scen=1:winddata[:S]],    - (networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
                                                                                                    + (networkdata[:Ψ]p_shed[:,t,scen])[l]
                                                                                                    + (networkdata[:Ψ]networkdata[:Cgens]*(sum(p_da[:,seg,t] for seg in 1:settings[:Nsteps]) .+ sum(p_rt[:,seg,t,scen] for seg in 1:settings[:Nsteps])))[l]
                                                                                                    + (networkdata[:Ψ]networkdata[:Cesr] *(b_da[:,t] .+ b_rt[:,t,scen]))[l]
                                                                                                    + (networkdata[:Ψ]winddata[:Cwind]*(w_da[:,t] .+ w_rt[:,t,scen]))[l] <=  networkdata[:f̅][l])
    @constraint(m, ρ̲_rt[l=1:length(networkdata[:lines]), t=1:settings[:T], scen=1:winddata[:S]],    - (networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
                                                                                                    + (networkdata[:Ψ]p_shed[:,t,scen])[l]
                                                                                                    + (networkdata[:Ψ]networkdata[:Cgens]*(sum(p_da[:,seg,t] for seg in 1:settings[:Nsteps]) .+ sum(p_rt[:,seg,t,scen] for seg in 1:settings[:Nsteps])))[l]
                                                                                                    + (networkdata[:Ψ]networkdata[:Cesr] *(b_da[:,t].+ b_rt[:,t,scen]))[l]
                                                                                                    + (networkdata[:Ψ]winddata[:Cwind]*(w_da[:,t] .+ w_rt[:,t,scen]))[l] >= -networkdata[:f̅][l])
    # solve
    @time optimize!(m)
    status = termination_status(m)
    @info("Stochastic LP status ---> $(status)")

    # return solution
    solution = Dict(
    :p_da          => JuMP.value.(p_da),
    :b_da          => JuMP.value.(b_da),
    :w_da          => JuMP.value.(w_da),
    :p_rt          => JuMP.value.(p_rt),
    :b_rt          => JuMP.value.(b_rt),
    :w_rt          => JuMP.value.(w_rt),
    :p_shed        => JuMP.value.(p_shed),
    :λ_da_sys      => JuMP.dual.(λ_da_sys),
    :λ_rt_sys      => JuMP.dual.(λ_rt_sys),
    :ρ̅_da          => JuMP.dual.(ρ̅_da),
    :ρ̲_da          => JuMP.dual.(ρ̲_da),
    :ρ̅_rt          => JuMP.dual.(ρ̅_rt),
    :ρ̲_rt          => JuMP.dual.(ρ̲_rt),
    :cost          => JuMP.objective_value.(m),
    :cost_ener     => (sum(sum(sum(pwl_cost_approx(settings,networkdata,g,0)[:slopes] .* JuMP.value.(p_da[g,:,t])) for g in 1:length(networkdata[:gens])) for t in 1:settings[:T])
                      + sum((ones(settings[:T]) * networkdata[:esrs][s].c_l)'JuMP.value.(b_da[s,:]) for s in 1:length(networkdata[:esrs]))),
    :cost_flex     => (sum((1/winddata[:S])*(sum(sum(sum(pwl_cost_approx(settings,networkdata,g,0)[:slopes] .* JuMP.value.(p_rt[g,:,t,scen])) for g in 1:length(networkdata[:gens])) for t in 1:settings[:T])
                            + sum((ones(settings[:T]) * networkdata[:esrs][s].c_l)'JuMP.value.(b_rt[s,:,scen]) for s in 1:length(networkdata[:esrs]))
                                    + sum((C_shed*ones(settings[:T]))'JuMP.value.(p_shed[b,:,scen]) for b in 1:length(networkdata[:buses]))) for scen in 1:winddata[:S])),
    :solvetime     => MOI.get(m, MOI.SolveTime()),
    :model         => m)
    return solution
end


function run_sto_LP_oos_reopt(settings,networkdata,winddata,loaddata,sol_sto_LP,oosdata,scen,rt_pen)
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

    @objective(m, Min,  sum((sum(sol_sto_LP[:p_da][g,seg,:] for seg in 1:settings[:Nsteps]))'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*(sum(sol_sto_LP[:p_da][g,seg,:] for seg in 1:settings[:Nsteps])) for g in 1:length(networkdata[:gens]))
                        + sum((r_p[g,:])'*diagm(ones(settings[:T])*networkdata[:gens][g].q*rt_pen)*(r_p[g,:]) for g in 1:length(networkdata[:gens]))
                        + sum((ones(settings[:T])*networkdata[:gens][g].c)' * (sum(sol_sto_LP[:p_da][g,seg,:] for seg in 1:settings[:Nsteps]))  for g in 1:length(networkdata[:gens]))
                        + sum((ones(settings[:T])*networkdata[:gens][g].c*rt_pen)' * (r_p[g,:])  for g in 1:length(networkdata[:gens]))
                       + sum((ones(settings[:T])*networkdata[:esrs][s].c_l)' * (sol_sto_LP[:b_da][s,:]) for s in 1:length(networkdata[:esrs]))
                       + sum((ones(settings[:T])*networkdata[:esrs][s].c_l*rt_pen)' * (r_b[s,:]) for s in 1:length(networkdata[:esrs]))
                       + sum((C_shed*ones(settings[:T]))'p_shed[d,:] for d in 1:length(networkdata[:buses]))
                       + sum((C_spill*ones(settings[:T]))'w_spill[w,:] for w in 1:Nw))

    #generator constraints rt
    @constraint(m, p_res_lim_mx[g=1:Ng, t=1:settings[:T]], sum(sol_sto_LP[:p_da][:,seg,:] for seg in 1:settings[:Nsteps])[g,t] + r_p[g,t] <= networkdata[:gens][g].p̅)
    @constraint(m, p_res_lim_mn[g=1:Ng, t=1:settings[:T]], sum(sol_sto_LP[:p_da][:,seg,:] for seg in 1:settings[:Nsteps])[g,t] + r_p[g,t] >= networkdata[:gens][g].p̲)
    @constraint(m, rese_lim_max[g=1:Ng, t=1:settings[:T]], -networkdata[:gens][g].r̅ <= r_p[g,t] <= networkdata[:gens][g].r̅)
    @constraint(m, p_res_ramp_a[g=1:Ng, t=2:settings[:T]], -networkdata[:gens][g].Ra̅ <= (sum(sol_sto_LP[:p_da][:,seg,:] for seg in 1:settings[:Nsteps])[g,t] + r_p[g,t]) - (sum(sol_sto_LP[:p_da][:,seg,:] for seg in 1:settings[:Nsteps])[g,t-1] + r_p[g,t-1]) <=  networkdata[:gens][g].Ra̅)

    #storage constraints
    @constraint(m, esr_chg_lim[s=1:Ns, t=1:settings[:T]],  sol_sto_LP[:b_da][s,t] + r_b[s,t] <= networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d)
    @constraint(m, esr_dis_lim[s=1:Ns, t=1:settings[:T]],  sol_sto_LP[:b_da][s,t] + r_b[s,t] >= - 1/networkdata[:esrs][s].η_c*networkdata[:esrs][s].P̅c)
    @constraint(m, esr_spill_up[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(sol_sto_LP[:b_da][s,t̂] + r_b[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].s̅)
    @constraint(m, esr_spill_dn[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(sol_sto_LP[:b_da][s,t̂] + r_b[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].s̲)
    @constraint(m, esr_eod_max[s=1:Ns, t=settings[:T]],    networkdata[:esrs][s].E𝟶 - sum(sol_sto_LP[:b_da][s,t̂] + r_b[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].E𝟶 + networkdata[:esrs][s].B_s)
    @constraint(m, esr_eod_min[s=1:Ns, t=settings[:T]],    networkdata[:esrs][s].E𝟶 - sum(sol_sto_LP[:b_da][s,t̂] + r_b[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].B_s)

    #load shedding and wind spillage limits
    @constraint(m, shed_lims[d=1:length(networkdata[:buses]), t=1:settings[:T]], 0 <=  p_shed[d,t]  <= loaddata[:node_loads][d,t])
    @constraint(m, spillage[w=1:Nw, t=1:settings[:T]],  0 <=  w_spill[w,t] <= winddata[:ŵ][w,t] - oosdata[:ξ_oos][:,scen,:][w,t])

    # #Power flow in the lines - PTDF formulation
    @constraint(m, plims_max[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]*(loaddata[:node_loads][:,t] - p_shed[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cgens]*(sum(sol_sto_LP[:p_da][:,seg,:] for seg in 1:settings[:Nsteps])[:,t] + r_p[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cesr]*(sol_sto_LP[:b_da][:,t] + r_b[:,t]))[l]
                                                                                 +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - oosdata[:ξ_oos][:,scen,:][:,t] - w_spill[:,t]))[l] <=  networkdata[:f̅][l])
    @constraint(m, plims_min[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]*(loaddata[:node_loads][:,t] - p_shed[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cgens]*(sum(sol_sto_LP[:p_da][:,seg,:] for seg in 1:settings[:Nsteps])[:,t] + r_p[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cesr]*(sol_sto_LP[:b_da][:,t] + r_b[:,t]))[l]
                                                                                 +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - oosdata[:ξ_oos][:,scen,:][:,t] - w_spill[:,t]))[l] >= -networkdata[:f̅][l])

    # # balance constraints
    @constraint(m, r_p_activation_balance[t=1:settings[:T]], sum(r_p[g,t] for g in 1:Ng)
                                                           + sum(r_b[s,t] for s in 1:Ns)
                                                           + sum(p_shed[b,t] for b in 1:length(networkdata[:buses]))
                                                           + sum((winddata[:ŵ][w,t] - oosdata[:ξ_oos][:,scen,:][w,t]) - sol_sto_LP[:w_da][w,t] - w_spill[w,t] for w in 1:Nw)
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
    :da_model   => "R2",
    )
    return solution
end
