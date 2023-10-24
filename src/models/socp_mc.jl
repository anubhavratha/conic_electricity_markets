function sto_socp_mc(settings,networkdata,winddata,loaddata)
    """ Chance Constrained Market Clearing, Mcc: SOCP problem """
    #Agents
    Ng = length(networkdata[:gens])
    Ns = length(networkdata[:esrs])
    Nw = size(winddata[:ŵ],1)
    m = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0, "QUIET" => true, "MSK_IPAR_INTPNT_SOLVE_FORM" => MSK_SOLVE_PRIMAL))
    # ---- DEFINE Primal VARIABLES ------ #
    #For power and participation factor
    @variable(m, p[1:Ng, 1:settings[:T]])
    @variable(m, α[1:Ng, 1:settings[:T]])
    @variable(m, b[1:Ns, 1:settings[:T]])
    @variable(m, γ[1:Ns, 1:settings[:T]])
    #Auxiliary Variables for converting QCP to SOCP
    @variable(m, z_p[1:Ng, 1:settings[:T]])
    @variable(m, z_α[1:Ng, 1:settings[:T]])

    # minimize primal objective
    @objective(m, Min, sum(sum((z_p[g,t]  + z_α[g,t] + networkdata[:gens][g].c  * p[g,t]) for g=1:Ng) for t=1:settings[:T])
                     + sum(sum(networkdata[:esrs][s].c_l * b[s,t] for s=1:Ns) for t=1:settings[:T]))

    @constraint(m, obj_ref_p[g=1:Ng, t=1:settings[:T]], [1/2; z_p[g,t]; sqrt(networkdata[:gens][g].q)*p[g,t]] in RotatedSecondOrderCone())
    @constraint(m, obj_ref_α[g=1:Ng, t=1:settings[:T]], [1/2; z_α[g,t]; sqrt(networkdata[:gens][g].q)*returnblks(winddata[:Σ],t,Nw)*ones(Nw)*α[g,t]] in RotatedSecondOrderCone())

    # generation constraints
    @constraint(m, gen_lim_max[g=1:Ng, t=1:settings[:T]], [networkdata[:gens][g].p̅ - p[g,t] ; Φ(settings,networkdata)*returnblks(winddata[:Σ],t,Nw)*ones(Nw)*α[g,t]] in SecondOrderCone())
    @constraint(m, gen_lim_min[g=1:Ng, t=1:settings[:T]], [p[g,t] - networkdata[:gens][g].p̲ ; Φ(settings,networkdata)*returnblks(winddata[:Σ],t,Nw)*ones(Nw)*α[g,t]] in SecondOrderCone())

    #reserve limit constraints
    @constraint(m, res_lim_max[g=1:Ng, t=1:settings[:T]], [networkdata[:gens][g].r̅; Φ(settings,networkdata)*returnblks(winddata[:Σ],t,Nw)*ones(Nw)*α[g,t]] in SecondOrderCone())
    @constraint(m, res_lim_min[g=1:Ng, t=1:settings[:T]], [networkdata[:gens][g].r̅; Φ(settings,networkdata)*returnblks(winddata[:Σ],t,Nw)*ones(Nw)*α[g,t]] in SecondOrderCone())

    if settings[:ramplimits] == true
        #ramping limit constraints
        @constraint(m, ram_lim_max[g=1:Ng, t=2:settings[:T]], [networkdata[:gens][g].Ra̅ + p[g,t-1] - p[g,t]; Φ(settings,networkdata)*returnblksramping(winddata[:Σ],t,Nw)*[-α[g,t-1]*ones(Nw) ;  α[g,t]*ones(Nw)]] in SecondOrderCone())
        @constraint(m, ram_lim_min[g=1:Ng, t=2:settings[:T]], [networkdata[:gens][g].Ra̅ - p[g,t-1] + p[g,t]; Φ(settings,networkdata)*returnblksramping(winddata[:Σ],t,Nw)*[ α[g,t-1]*ones(Nw) ; -α[g,t]*ones(Nw)]] in SecondOrderCone())
    end

    #balance constraints
    @constraint(m, λ_e[t=1:settings[:T]], sum(p[g,t] for g=1:Ng) + sum(b[s,t] for s=1:Ns) + sum(winddata[:ŵ][k,t] for k=1:Nw) - loaddata[:D][t] == 0)

    settings[:det] == false ? @constraint(m, λ_r[t=1:settings[:T]], sum(α[g,t] for g=1:Ng) + sum(γ[s,t] for s=1:Ns)  - 1 == 0) : NaN
    if settings[:det] == true
        @constraint(m, alpha_set_zero[g=1:Ng, t=1:settings[:T]], α[g,t]==0)
        @constraint(m, gamma_set_zero[s=1:Ns, t=1:settings[:T]], γ[s,t]==0)
        λ_r = zeros(settings[:T])
    end

    #ESR charging and discharging max constraints
    @constraint(m, esr_chg_lim[s=1:Ns, t=1:settings[:T]], [networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d - b[s,t]     ; Φ(settings,networkdata)*returnblks(winddata[:Σ],t,Nw)*ones(Nw)*γ[s,t]] in SecondOrderCone())
    @constraint(m, esr_dis_lim[s=1:Ns, t=1:settings[:T]], [b[s,t] + (1/networkdata[:esrs][s].η_c)*networkdata[:esrs][s].P̅c ; Φ(settings,networkdata)*returnblks(winddata[:Σ],t,Nw)*ones(Nw)*γ[s,t]] in SecondOrderCone())

    #ESR SOC limits minimum and maximum constraints
    @constraint(m, esr_spill_up[s=1:Ns, t=1:settings[:T]],  [networkdata[:esrs][s].s̅ - networkdata[:esrs][s].E𝟶 + sum(b[s,t̂] for t̂ in 1:t); Φ(settings,networkdata)*returnblksstorage(winddata[:Σ],t,Nw)*repeat(γ[s,1:t],inner=Nw)] in SecondOrderCone())
    @constraint(m, esr_spill_dn[s=1:Ns, t=1:settings[:T]],  [networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].s̲ - sum(b[s,t̂] for t̂ in 1:t); Φ(settings,networkdata)*returnblksstorage(winddata[:Σ],t,Nw)*repeat(γ[s,1:t],inner=Nw)] in SecondOrderCone())

    # #ESR terminal hour energy level constraints
    @constraint(m, esr_eod_max[s=1:Ns, t=settings[:T]], [networkdata[:esrs][s].B_s + sum(b[s,t̂] for t̂ in 1:t) ; Φ(settings,networkdata)*returnblksstorage(winddata[:Σ],settings[:T],Nw)*repeat(γ[s,1:settings[:T]], inner=Nw)] in SecondOrderCone())
    @constraint(m, esr_eod_min[s=1:Ns, t=settings[:T]], [networkdata[:esrs][s].B_s - sum(b[s,t̂] for t̂ in 1:t) ; Φ(settings,networkdata)*returnblksstorage(winddata[:Σ],settings[:T],Nw)*repeat(γ[s,1:settings[:T]], inner=Nw)] in SecondOrderCone())

    #Power flow in the lines - PTDF formulation
    ρ = Array{JuMP.ConstraintRef, 2}(undef, length(networkdata[:lines])*2, settings[:T])
    for t in 1:settings[:T]
        awinline = (networkdata[:Ψ]*networkdata[:Cgens]*p[:,t]
                 + networkdata[:Ψ]*networkdata[:Cesr]  *b[:,t]
                 + networkdata[:Ψ]*winddata[:Cwind]    *winddata[:ŵ][:,t]
                 - networkdata[:Ψ]*loaddata[:node_loads][:,t])
        bwinline = (networkdata[:Ψ]*networkdata[:Cgens]*α[:,t]*ones(1,Nw)
                  + networkdata[:Ψ]*networkdata[:Cesr] *γ[:,t]*ones(1,Nw)
                  - networkdata[:Ψ]*winddata[:Cwind])
        for l in 1:length(networkdata[:lines])
            ρ[l,t]                             = @constraint(m, [networkdata[:f̅][l] - awinline[l]; Φ(settings,networkdata) * returnblks(winddata[:Σ],t,Nw) * (bwinline[l,:])] in SecondOrderCone())
            ρ[length(networkdata[:lines])+l,t] = @constraint(m, [networkdata[:f̅][l] + awinline[l]; Φ(settings,networkdata) * returnblks(winddata[:Σ],t,Nw) * (bwinline[l,:])] in SecondOrderCone())
            set_name(ρ[l,t], "ρ_l$(l)_t$(t)")
        end
    end

    # solve
    @time optimize!(m)
    status = termination_status(m)
    @info("Primal Spatial SOCP status ---> $(status)")

    nomflows = zeros(length(networkdata[:lines]), settings[:T])
    for t in 1:settings[:T]
        nomflows[:,t] = networkdata[:Ψ]*(loaddata[:node_loads][:,t] - networkdata[:Cgens]*JuMP.value.(p[:,t]) - networkdata[:Cesr]*JuMP.value.(b[:,t]) - winddata[:Cwind]*winddata[:ŵ][:,t])
    end

    # return solution
    solution = Dict(
    :p          => JuMP.value.(p),
    :z_p        => JuMP.value.(z_p),
    :z_α        => JuMP.value.(z_α),
    :α          => JuMP.value.(α),
    :b          => JuMP.value.(b),
    :γ          => JuMP.value.(γ),
    :nomflows   => nomflows,
    :λ_E        => JuMP.dual.(λ_e),
    settings[:det] == true ? :λ_R => 0 : :λ_R        => JuMP.dual.(λ_r),
    :ρ          => JuMP.dual.(ρ),
    :cost       => JuMP.objective_value.(m),
    :cost_ener  => sum(JuMP.value.(z_p)) + sum(sum(networkdata[:gens][g].c * JuMP.value.(p[g,t]) for g=1:Ng) for t=1:settings[:T]) + sum(sum(networkdata[:esrs][s].c_l * JuMP.value.(b[s,t]) for s=1:Ns) for t=1:settings[:T]),
    :cost_flex  => sum(JuMP.value.(z_α)),
    :ŵ          => winddata[:ŵ],
    :D          => loaddata[:D],
    :solvetime  => JuMP.MOI.get(m, JuMP.MOI.SolveTimeSec()),
    :model      => m
    )
    return solution
end

function compute_ex_ante_violprob(settings,networkdata,winddata,sol_sto_SOCP)
    """ Ex-Ante constraint violation probability and expected cost - Without Reoptimization """
    p̃ = zeros(length(networkdata[:gens]), winddata[:S], settings[:T])
    b̃ = zeros(length(networkdata[:esrs]), winddata[:S], settings[:T])
    for t in 1:settings[:T]
        for scen in 1:winddata[:S]
            p̃[:,scen,t] .= sol_sto_SOCP[:p][:,t] .+ sum(winddata[:ξ][:,scen,t]) * sol_sto_SOCP[:α][:,t]
            b̃[:,scen,t] .= sol_sto_SOCP[:b][:,t] .+ sum(winddata[:ξ][:,scen,t]) * sol_sto_SOCP[:γ][:,t]
        end
    end
    exp_cost = zeros(winddata[:S])
    for scen in 1:winddata[:S]
        if settings[:ESRs] == true
            exp_cost[scen] = (sum(p̃[g,scen,:]'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*p̃[g,scen,:] + (ones(settings[:T])*networkdata[:gens][g].c)'p̃[g,scen,:] for g in 1:length(networkdata[:gens]))
                        + sum(b̃[s,scen,:]'*diagm(ones(settings[:T])*networkdata[:esrs][s].c_q)*b̃[s,scen,:]  + (ones(settings[:T])*networkdata[:esrs][s].c_l)'b̃[s,scen,:] for s in 1:length(networkdata[:esrs])))
        elseif settings[:ESRs] == false
            exp_cost[scen] = sum(p̃[g,scen,:]'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*p̃[g,scen,:] + (ones(settings[:T])*networkdata[:gens][g].c)'p̃[g,scen,:] for g in 1:length(networkdata[:gens]))
        else
            @warn("Something wrong with computing expected cost in violation probability calculation!")
        end
    end
    inf_flag = zeros(winddata[:S])
    num_tolerance = 0.00001
    for scen in 1:winddata[:S]
        for g in 1:length(networkdata[:gens])
            for t in 1:settings[:T]
                p̃[g,scen,t] >= networkdata[:gens][g].p̅ + num_tolerance ? inf_flag[scen] = 1 : NaN
                p̃[g,scen,t] <= networkdata[:gens][g].p̲ - num_tolerance ? inf_flag[scen] = 1 : NaN
                sum(winddata[:ξ][:,scen,t]) * sol_sto_SOCP[:α][g,t] >=  networkdata[:gens][g].r̅ + num_tolerance ? inf_flag[scen] = 1 : NaN
                sum(winddata[:ξ][:,scen,t]) * sol_sto_SOCP[:α][g,t] <= -networkdata[:gens][g].r̅ - num_tolerance ? inf_flag[scen] = 1 : NaN
            end
            for t in 2:settings[:T]
                p̃[g,scen,t] - p̃[g,scen,t-1] >=  networkdata[:gens][g].Ra̅ + num_tolerance ? inf_flag[scen] = 1 : NaN
                p̃[g,scen,t] - p̃[g,scen,t-1] <= -networkdata[:gens][g].Ra̅ - num_tolerance ? inf_flag[scen] = 1 : NaN
            end
        end
        for s in 1:length(networkdata[:esrs])
            for t in 1:settings[:T]
                b̃[s,scen,t] >=  networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d + num_tolerance ? inf_flag[scen] = 1 : NaN
                b̃[s,scen,t] <= -(1/networkdata[:esrs][s].η_c)*networkdata[:esrs][s].P̅c - num_tolerance ? inf_flag[scen] = 1 : NaN
                networkdata[:esrs][s].E𝟶 - sum(b̃[s,scen,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].s̅ + num_tolerance ? inf_flag[scen] = 1 : NaN
                networkdata[:esrs][s].E𝟶 - sum(b̃[s,scen,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].s̲ - num_tolerance ? inf_flag[scen] = 1 : NaN
            end
            networkdata[:esrs][s].E𝟶 - sum(b̃[s,scen,t̂] for t̂ in 1:settings[:T]) >= networkdata[:esrs][s].E𝟶 + networkdata[:esrs][s].B_s + num_tolerance ? inf_flag[scen] = 1 : NaN
            networkdata[:esrs][s].E𝟶 - sum(b̃[s,scen,t̂] for t̂ in 1:settings[:T]) <= networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].B_s - num_tolerance ? inf_flag[scen] = 1 : NaN
        end
        for l in 1:length(networkdata[:lines])
            for t in 1:settings[:T]
             (-(networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
              +(networkdata[:Ψ]networkdata[:Cgens]*(p̃[:,scen,t]))[l]
              +(networkdata[:Ψ]networkdata[:Cesr]*(b̃[:,scen,t]))[l]
              +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - winddata[:ξ][:,scen,t]))[l]) >=  networkdata[:f̅][l] + num_tolerance ? inf_flag[scen]=1 : NaN

             (-(networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
              +(networkdata[:Ψ]networkdata[:Cgens]*(p̃[:,scen,t]))[l]
              +(networkdata[:Ψ]networkdata[:Cesr]*(b̃[:,scen,t]))[l]
              +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - winddata[:ξ][:,scen,t]))[l]) <= -networkdata[:f̅][l] - num_tolerance ? inf_flag[scen]=1 : NaN
           end
        end
    end
    @info("empirical violations ---> $(sum(inf_flag))")
    @info("empirical violation probability ---> $(sum(inf_flag)/winddata[:S]*100)%")
    return Dict(:cost => mean(exp_cost),  :ε_stat => sum(inf_flag)/winddata[:S], :p̃ => p̃, :b̃ => b̃)
end

function compute_ex_post_violprob(settings,networkdata,winddata,sol_sto_SOCP,oosdata)
    """ Out-of-sample Simulations with re-optimization"""
    p̃ = zeros(length(networkdata[:gens]), oosdata[:S], settings[:T])
    b̃ = zeros(length(networkdata[:esrs]), oosdata[:S], settings[:T])
    exp_cost = zeros(winddata[:S])
    for t in 1:settings[:T]
        for scen in 1:oosdata[:S]
            p̃[:,scen,t] .= sol_sto_SOCP[:p][:,t] .+ sum(oosdata[:ξ_oos][:,scen,t]) * sol_sto_SOCP[:α][:,t]
            b̃[:,scen,t] .= sol_sto_SOCP[:b][:,t] .+ sum(oosdata[:ξ_oos][:,scen,t]) * sol_sto_SOCP[:γ][:,t]
        end
    end
    exp_cost = zeros(oosdata[:S])
    for scen in 1:oosdata[:S]
        exp_cost[scen] = (sum(p̃[g,scen,:]'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*p̃[g,scen,:] + (ones(settings[:T])*networkdata[:gens][g].c)'p̃[g,scen,:] for g in 1:length(networkdata[:gens]))
                        + sum(b̃[s,scen,:]'*diagm(ones(settings[:T])*networkdata[:esrs][s].c_q)*b̃[s,scen,:]  + (ones(settings[:T])*networkdata[:esrs][s].c_l)'b̃[s,scen,:] for s in 1:length(networkdata[:esrs])))
    end

    inf_flag = zeros(oosdata[:S])
    num_tolerance = 0.00001
    for scen in 1:oosdata[:S]
        for g in 1:length(networkdata[:gens])
            for t in 1:settings[:T]
                p̃[g,scen,t] >= networkdata[:gens][g].p̅ + num_tolerance ? inf_flag[scen] = 1 : NaN
                p̃[g,scen,t] <= networkdata[:gens][g].p̲ - num_tolerance ? inf_flag[scen] = 1 : NaN
                sum(oosdata[:ξ_oos][:,scen,t]) * sol_sto_SOCP[:α][g,t] >=  networkdata[:gens][g].r̅ + num_tolerance ? inf_flag[scen] = 1 : NaN
                sum(oosdata[:ξ_oos][:,scen,t]) * sol_sto_SOCP[:α][g,t] <= -networkdata[:gens][g].r̅ - num_tolerance ? inf_flag[scen] = 1 : NaN
            end
            for t in 2:settings[:T]
                p̃[g,scen,t] - p̃[g,scen,t-1] >=  networkdata[:gens][g].Ra̅ + num_tolerance ? inf_flag[scen] = 1 : NaN
                p̃[g,scen,t] - p̃[g,scen,t-1] <= -networkdata[:gens][g].Ra̅ - num_tolerance ? inf_flag[scen] = 1 : NaN
            end
        end
        for s in 1:length(networkdata[:esrs])
            for t in 1:settings[:T]
                b̃[s,scen,t] >=  networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d + num_tolerance ? inf_flag[scen] = 1 : NaN
                b̃[s,scen,t] <= -(1/networkdata[:esrs][s].η_c)*networkdata[:esrs][s].P̅c - num_tolerance ? inf_flag[scen] = 1 : NaN
                networkdata[:esrs][s].E𝟶 - sum(b̃[s,scen,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].s̅ + num_tolerance ? inf_flag[scen] = 1 : NaN
                networkdata[:esrs][s].E𝟶 - sum(b̃[s,scen,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].s̲ - num_tolerance ? inf_flag[scen] = 1 : NaN
            end
            networkdata[:esrs][s].E𝟶 - sum(b̃[s,scen,t̂] for t̂ in 1:settings[:T]) >= networkdata[:esrs][s].E𝟶 + networkdata[:esrs][s].B_s + num_tolerance ? inf_flag[scen] = 1 : NaN
            networkdata[:esrs][s].E𝟶 - sum(b̃[s,scen,t̂] for t̂ in 1:settings[:T]) <= networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].B_s - num_tolerance ? inf_flag[scen] = 1 : NaN
        end
        for l in 1:length(networkdata[:lines])
            for t in 1:settings[:T]
             (-(networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
              +(networkdata[:Ψ]networkdata[:Cgens]*(p̃[:,scen,t]))[l]
              +(networkdata[:Ψ]networkdata[:Cesr]*(b̃[:,scen,t]))[l]
              +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - oosdata[:ξ_oos][:,scen,t]))[l]) >=  networkdata[:f̅][l] + num_tolerance ? inf_flag[scen]=1 : NaN

             (-(networkdata[:Ψ]loaddata[:node_loads][:,t])[l]
              +(networkdata[:Ψ]networkdata[:Cgens]*(p̃[:,scen,t]))[l]
              +(networkdata[:Ψ]networkdata[:Cesr]*(b̃[:,scen,t]))[l]
              +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - oosdata[:ξ_oos][:,scen,t]))[l]) <= -networkdata[:f̅][l] - num_tolerance ? inf_flag[scen]=1 : NaN
           end
        end
    end
    @info("empirical violations ---> $(sum(inf_flag))")
    @info("empirical violation probability ---> $(sum(inf_flag)/oosdata[:S]*100)%")
    ex_post_res = Dict(:cost => mean(exp_cost),  :ε_stat => sum(inf_flag)/oosdata[:S], :p̃ => p̃, :b̃ => b̃)

    return ex_post_res
end


function run_sto_SOCP_oos_reopt(settings,networkdata,winddata,sol_sto_SOCP,oosdata,scen)
    """ Real-time reoptimization """
    #Agents
    Ng = length(networkdata[:gens])
    Ns = length(networkdata[:esrs])
    Nw = size(winddata[:ŵ],1)
    C_shed  = 500
    C_spill = 500
    Δ = sum(oosdata[:ξ_oos][w,scen,:] for w in 1:Nw)

    m = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0, "QUIET" => true))
    # ---- Define Primal VARIABLES ------ #
    #For power and participation factor
    @variable(m, r_p[1:Ng, 1:settings[:T]])    #real-time adjustment
    @variable(m, r_b[1:Ns, 1:settings[:T]])    #real-time adjustment
    @variable(m, p_shed[1:length(networkdata[:buses]), 1:settings[:T]]) #load_shedding to ensure feasibility
    @variable(m, w_spill[1:Nw,1:settings[:T]])      #windspillage to ensure feasibility

    @objective(m, Min,  sum(((sol_sto_SOCP[:p][g,:] + r_p[g,:])'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*(sol_sto_SOCP[:p][g,:] + r_p[g,:]) + (ones(settings[:T])*networkdata[:gens][g].c)' * (sol_sto_SOCP[:p][g,:] + r_p[g,:])  for g in 1:length(networkdata[:gens])))
                       + sum((ones(settings[:T])*networkdata[:esrs][s].c_l)' * (sol_sto_SOCP[:b][s,:] + r_b[s,:]) for s in 1:length(networkdata[:esrs]))
                       + sum((C_shed*ones(settings[:T]))'p_shed[d,:] for d in 1:length(networkdata[:buses]))
                       + sum((C_spill*ones(settings[:T]))'w_spill[w,:] for w in 1:Nw))

    @constraint(m, p_res_act[g=1:Ng,    t=1:settings[:T]],  r_p[g,t] == Δ[t] * sol_sto_SOCP[:α][g,t])
    @constraint(m, p_res_lim_mx[g=1:Ng, t=1:settings[:T]],  sol_sto_SOCP[:p][g,t] + r_p[g,t] <= networkdata[:gens][g].p̅)
    @constraint(m, p_res_lim_mn[g=1:Ng, t=1:settings[:T]],  sol_sto_SOCP[:p][g,t] + r_p[g,t] >= networkdata[:gens][g].p̲)
    @constraint(m, rese_lim_max[g=1:Ng, t=1:settings[:T]], -networkdata[:gens][g].r̅ <= r_p[g,t] <= networkdata[:gens][g].r̅)
    @constraint(m, p_res_ramp_a[g=1:Ng, t=2:settings[:T]], -networkdata[:gens][g].Ra̅ <= (sol_sto_SOCP[:p][g,t] + r_p[g,t]) - (sol_sto_SOCP[:p][g,t-1] + r_p[g,t-1]) <=  networkdata[:gens][g].Ra̅)

    #storage constraints
    @constraint(m, b_res_act[s=1:Ns,    t=1:settings[:T]], r_b[s,t] == Δ[t] * sol_sto_SOCP[:γ][s,t])
    @constraint(m, esr_chg_lim[s=1:Ns,  t=1:settings[:T]], sol_sto_SOCP[:b][s,t] + r_b[s,t] <= networkdata[:esrs][s].η_d*networkdata[:esrs][s].P̅d)
    @constraint(m, esr_dis_lim[s=1:Ns,  t=1:settings[:T]], sol_sto_SOCP[:b][s,t] + r_b[s,t] >= - 1/networkdata[:esrs][s].η_c*networkdata[:esrs][s].P̅c)
    @constraint(m, esr_spill_up[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(sol_sto_SOCP[:b][s,t̂] + r_b[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].s̅)
    @constraint(m, esr_spill_dn[s=1:Ns, t=1:settings[:T]], networkdata[:esrs][s].E𝟶 - sum(sol_sto_SOCP[:b][s,t̂] + r_b[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].s̲)
    @constraint(m, esr_eod_max[s=1:Ns,  t=settings[:T]],   networkdata[:esrs][s].E𝟶 - sum(sol_sto_SOCP[:b][s,t̂] + r_b[s,t̂] for t̂ in 1:t) <= networkdata[:esrs][s].E𝟶 + networkdata[:esrs][s].B_s)
    @constraint(m, esr_eod_min[s=1:Ns,  t=settings[:T]],   networkdata[:esrs][s].E𝟶 - sum(sol_sto_SOCP[:b][s,t̂] + r_b[s,t̂] for t̂ in 1:t) >= networkdata[:esrs][s].E𝟶 - networkdata[:esrs][s].B_s)

    #load shedding and wind spillage limits
    @constraint(m, shed_lims[d=1:length(networkdata[:buses]), t=1:settings[:T]], 0 <=  p_shed[d,t]  <= loaddata[:node_loads][d,t])
    @constraint(m, spillage[w=1:Nw, t=1:settings[:T]],  0 <=  w_spill[w,t] <= winddata[:ŵ][w,t] - oosdata[:ξ_oos][:,scen,:][w,t])

    # #Power flow in the lines - PTDF formulation
    @constraint(m, plims_max[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]*(loaddata[:node_loads][:,t] - p_shed[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cgens]*(sol_sto_SOCP[:p][:,t] + r_p[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cesr] *(sol_sto_SOCP[:b][:,t] + r_b[:,t]))[l]
                                                                                 +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - oosdata[:ξ_oos][:,scen,t] - w_spill[:,t]))[l] <=  networkdata[:f̅][l])
    @constraint(m, plims_min[l=1:length(networkdata[:lines]), t=1:settings[:T]], -(networkdata[:Ψ]*(loaddata[:node_loads][:,t] - p_shed[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cgens]*(sol_sto_SOCP[:p][:,t] + r_p[:,t]))[l]
                                                                                 +(networkdata[:Ψ]networkdata[:Cesr] *(sol_sto_SOCP[:b][:,t] + r_b[:,t]))[l]
                                                                                 +(networkdata[:Ψ]winddata[:Cwind]*(winddata[:ŵ][:,t] - oosdata[:ξ_oos][:,scen,t] - w_spill[:,t]))[l] >=  -networkdata[:f̅][l])

    # # balance constraints
    @constraint(m, r_p_act[t=1:settings[:T]], sum(r_p[g,t] for g in 1:Ng) + sum(r_b[s,t] for s in 1:Ns) + sum(p_shed[b,t] for b in 1:length(networkdata[:buses])) == Δ[t] + sum(w_spill[k,t] for k in 1:Nw))

    optimize!(m)

    solution = Dict(
    :cost    => JuMP.objective_value(m),
    :r_p        => JuMP.value.(r_p),
    :r_b        => JuMP.value.(r_b),
    :p_shed     => JuMP.value.(p_shed),
    :w_spill    => JuMP.value.(w_spill),
    :status     => termination_status(m),
    :solvetime  => JuMP.MOI.get(m, JuMP.MOI.SolveTimeSec()),
    :da_model   => "Mcc",
    )
    return solution
end
