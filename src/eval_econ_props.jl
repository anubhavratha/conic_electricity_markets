""" Compute LMPS """
function compute_lmps(settings,networkdata,sol)
    """computes LMPS based on duality"""
    if :ρ in keys(sol)
        #SOCP market clearing
        ρ̲ = zeros(length(networkdata[:lines]), settings[:T])
        ρ̅ = zeros(length(networkdata[:lines]), settings[:T])
        for t in 1:settings[:T]
            for l in 1:length(networkdata[:lines])
                ρ̅[l,t] = sol[:ρ][l,t][1]
            end
            for l in 1:length(networkdata[:lines])
                ρ̲[l,t] = sol[:ρ][l+length(networkdata[:lines]),t][1]
            end
        end
        Π_E = repeat(sol[:λ_E]', inner = (length(networkdata[:buses]),1)) .- networkdata[:Ψ]' * (ρ̅ .- ρ̲) * 1/Φ(settings,networkdata)
        Π_R = repeat(sol[:λ_R]', inner = (length(networkdata[:buses]),1)) .- networkdata[:Ψ]' * (ρ̅ .- ρ̲) * 1/Φ(settings,networkdata)
    elseif :ρ̅_da in keys(sol)
        #LP market clearing
        Π_E = repeat(sol[:λ_da_sys]', inner = (length(networkdata[:buses]),1)) .- networkdata[:Ψ]' * (sol[:ρ̅_da] .- sol[:ρ̲_da])
        Π_R = zeros(length(networkdata[:buses]), settings[:T], winddata[:S])
        for scen in 1:winddata[:S]
            Π_R[:,:,scen] = repeat(sol[:λ_rt_sys][:,scen]', inner = (length(networkdata[:buses]),1)) .- networkdata[:Ψ]' * (sol[:ρ̅_rt][:,:,scen] .- sol[:ρ̲_rt][:,:,scen])
        end
    else
        @warn("something went wrong while extracting LMPs")
    end
    Π_comms = Dict(:Π_E => Π_E, :Π_R => Π_R)
    return Π_comms
end

function analyze_gen_flex_payoffs(settings,networkdata,sol_det,sol_SOCP)
    """Determine the perceived price of flexibility payoffs for generators """
    Ng = length(networkdata[:gens])
    flex_payment  = zeros(Ng,settings[:T])
    flex_pay_rate = zeros(Ng,settings[:T])
    for t in 1:settings[:T]
        for g in 1:Ng
            if sol_SOCP[:α][g,t] > 1e-6
                flex_payment[g,t]   = sol_SOCP[:Π][:Π_R][networkdata[:gens][g].node,t] * sol_SOCP[:α][g,t]
                flex_pay_rate[g,t]  = flex_payment[g,t]/abs(sol_SOCP[:p][g,t] - sol_det[:p][g,t])
            else
                flex_payment[g,t]   = 0
                flex_pay_rate[g,t]  = 0
            end
        end
    end
    return Dict(:flex_pay => flex_payment, :flex_pay_rate => flex_pay_rate)
end

function analyze_esr_flex_payoffs(settings,networkdata,sol_det,sol_SOCP)
    """Determine the perceived price of flexibility payoffs for ESRs """
    Ns = length(networkdata[:esrs])
    flex_payment    = zeros(Ns,settings[:T])
    flex_pay_rate   = zeros(Ns,settings[:T])
    for t in 1:settings[:T]
        for s in 1:Ns
            if sol_SOCP[:γ][s,t] > 1e-6
                flex_payment[s,t]   = sol_SOCP[:Π][:Π_R][networkdata[:esrs][s].node,t] * sol_SOCP[:γ][s,t]
                flex_pay_rate[s,t]  = flex_payment[s,t]/abs(sol_SOCP[:b][s,t] - sol_det[:b][s,t])
            else
                flex_payment[s,t]   = 0
                flex_pay_rate[s,t]  = 0
            end
        end
    end
    return Dict(:flex_pay => flex_payment, :flex_pay_rate => flex_pay_rate)
end



""" COST RECOVERY """
#1. Generator Profits - Expected Cost Recovery
function process_gen_outcomes(settings,networkdata,sol)
    Ng = length(networkdata[:gens])
    Rev_gen = zeros(Ng)
    Cos_gen = zeros(Ng)
    Pro_gen = zeros(Ng)
    if :α in keys(sol)
        for g in 1:Ng
            Cos_gen[g] = round(sum(sol[:z_p][g,:]) + sum(sol[:z_α][g,:]) +  (ones(settings[:T])*networkdata[:gens][g].c)'sol[:p][g,:], digits=4)
            Rev_gen[g] = round(sol[:Π][:Π_E][networkdata[:gens][g].node,:]'sol[:p][g,:] + sol[:Π][:Π_R][networkdata[:gens][g].node,:]'sol[:α][g,:], digits=4)
            Pro_gen[g] = round(Rev_gen[g] - Cos_gen[g], digits=2)
        end
    elseif :p_da in keys(sol)
        for g in 1:Ng
            Cos_gen[g] = (sum(sol[:p_da][g,seg,:] for seg in 1:settings[:Nsteps])' * diagm(ones(settings[:T]) * networkdata[:gens][g].q) * sum(sol[:p_da][g,seg,:] for seg in 1:settings[:Nsteps]) + (ones(settings[:T]) * networkdata[:gens][g].c)' * sum(sol[:p_da][g,seg,:] for seg in 1:settings[:Nsteps])
                            + sum((1/winddata[:S])*((sum(sol[:p_rt][g,seg,:,scen] for seg in 1:settings[:Nsteps]))' * diagm(ones(settings[:T]) * networkdata[:gens][g].q)   * (sum(sol[:p_rt][g,seg,:,scen] for seg in 1:settings[:Nsteps])) + (ones(settings[:T])*networkdata[:gens][g].c)'* ((sum(sol[:p_rt][g,seg,:,scen] for seg in 1:settings[:Nsteps])))) for scen in 1:winddata[:S]))
            Rev_gen[g] = sol[:Π][:Π_E][networkdata[:gens][g].node,:]'*(sum(sol[:p_da][g,seg,:] for seg in 1:settings[:Nsteps])) + sum(sol[:Π][:Π_R][networkdata[:gens][g].node,:,scen]'*(sum(sol[:p_rt][g,seg,:,scen] for seg in 1:settings[:Nsteps])) for scen in 1:winddata[:S])
            Pro_gen[g] = round.(Rev_gen[g] - Cos_gen[g], digits=4)
        end
    else
        @warn("something is wrong with verifying cost recovery")
    end
    gen_outcomes = Dict(:Pro_gen => Pro_gen, :Cos_gen => Cos_gen, :Rev_gen => Rev_gen)
    return gen_outcomes
end

#Cost recovery by scenario
function check_cost_rec_scen(settings,networkdata,sol_da_mc,sol_rt_opt)
    """Checks cost recovery by scenario"""
    Ng = length(networkdata[:gens])
    Rev_gen = zeros(Ng)
    Cos_gen = zeros(Ng)
    Pro_gen = zeros(Ng)
    if sol_rt_opt[:da_model] == "Mcc"
        for g in 1:Ng
            Cos_gen[g] = (sol_da_mc[:p][g,:] + sol_rt_opt[:r_p][g,:])'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*(sol_da_mc[:p][g,:] + sol_rt_opt[:r_p][g,:]) + (ones(settings[:T])*networkdata[:gens][g].c)' * (sol_da_mc[:p][g,:] + sol_rt_opt[:r_p][g,:])
            Rev_gen[g] = sol_da_mc[:Π][:Π_E][networkdata[:gens][g].node,:]'*(sol_da_mc[:p][g,:] + sol_rt_opt[:r_p][g,:]) + sol_da_mc[:Π][:Π_R][networkdata[:gens][g].node,:]'sol_da_mc[:α][g,:]
            Pro_gen[g] = round.(Rev_gen[g] - Cos_gen[g], digits=4)
        end
    elseif sol_rt_opt[:da_model] == "R2"
        for g in 1:Ng
            Cos_gen[g] = sum(((sum(sol_da_mc[:p_da][g,seg,:] for seg in 1:settings[:Nsteps]) + sol_rt_opt[:r_p][g,:])'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*(sum(sol_da_mc[:p_da][g,seg,:] for seg in 1:settings[:Nsteps]) + sol_rt_opt[:r_p][g,:]) + (ones(settings[:T])*networkdata[:gens][g].c)' * (sum(sol_da_mc[:p_da][g,seg,:] for seg in 1:settings[:Nsteps]) + sol_rt_opt[:r_p][g,:])))
            Rev_gen[g] = sol_da_mc[:Π][:Π_E][networkdata[:gens][g].node,:]'*(sum(sol_da_mc[:p_da][g,seg,:] for seg in 1:settings[:Nsteps]) + sol_rt_opt[:r_p][g,:])
            Pro_gen[g] = round.(Rev_gen[g] - Cos_gen[g], digits=4)
        end
    elseif sol_rt_opt[:da_model] == "R1"
        for g in 1:Ng
            Cos_gen[g] = sum(((sum(sol_da_mc[:p][g,seg,:] for seg in 1:settings[:Nsteps]) + sol_rt_opt[:r_p][g,:])'*diagm(ones(settings[:T])*networkdata[:gens][g].q)*(sum(sol_da_mc[:p][g,seg,:] for seg in 1:settings[:Nsteps]) + sol_rt_opt[:r_p][g,:]) + (ones(settings[:T])*networkdata[:gens][g].c)' * (sum(sol_da_mc[:p][g,seg,:] for seg in 1:settings[:Nsteps]) + sol_rt_opt[:r_p][g,:])))
            Rev_gen[g] = sol_da_mc[:Π][:Π_E][networkdata[:gens][g].node,:]'*(sum(sol_da_mc[:p][g,seg,:] for seg in 1:settings[:Nsteps]) + sol_rt_opt[:r_p][g,:])
            Pro_gen[g] = round.(Rev_gen[g] - Cos_gen[g], digits=4)
        end
    else
        @warn("something is wrong with verifying cost recovery")
    end
    gen_outcomes = Dict(:Pro_gen => Pro_gen, :Cos_gen => Cos_gen, :Rev_gen => Rev_gen)
    return gen_outcomes
end

#2. ESR Profits
function process_esr_outcomes(settings,networkdata,sol)
    Ns = length(networkdata[:esrs])
    Rev_esr = zeros(Ns)
    Cos_esr = zeros(Ns)
    Pro_esr = zeros(Ns)
    if :γ in keys(sol)
        for s in 1:Ns
            Cos_esr[s] = (ones(settings[:T])*networkdata[:esrs][s].c_l)'sol[:b][s,:]
            Rev_esr[s] = sol[:Π][:Π_E][networkdata[:esrs][s].node,:]'sol[:b][s,:] + sol[:Π][:Π_R][networkdata[:esrs][s].node,:]'sol[:γ][s,:]
            Pro_esr[s] = round.(Rev_esr[s] - Cos_esr[s], digits=4)
        end
    elseif :b_da in keys(sol)
        for s in 1:Ns
            Cos_esr[s] = sol[:b_da][s,:]' * diagm(ones(settings[:T]) * networkdata[:esrs][s].c_q) * sol[:b_da][s,:] + (ones(settings[:T]) * networkdata[:esrs][s].c_l)'sol[:b_da][s,:] + sum((1/winddata[:S])*(sol[:b_rt][s,:,scen]' * diagm(ones(settings[:T]) * networkdata[:esrs][s].c_q)   * sol[:b_rt][s,:,scen] + (ones(settings[:T])*networkdata[:esrs][s].c_l)'sol[:b_rt][s,:,scen]) for scen in 1:winddata[:S])
            Rev_esr[s] = sol[:Π][:Π_E][networkdata[:esrs][s].node,:]'sol[:b_da][s,:] + sum(sol[:Π][:Π_R][networkdata[:esrs][s].node,:,scen]'sol[:b_rt][s,:,scen] for scen in 1:winddata[:S])
            Pro_esr[s] = round.(Rev_esr[s] - Cos_esr[s], digits=4)
        end
    else
        @warn("something is wrong with verifying cost recovery")
    end
    esr_outcomes = Dict(:Pro_esr => Pro_esr, :Cos_esr => Cos_esr, :Rev_esr => Rev_esr)
    return esr_outcomes
end

# p1_esr = groupedbar(["S1","S1","S2","S2","S3","S3"],
#             [esr_outcomes_SOCP[:Rev_esr] esr_outcomes_LP[:Rev_esr]],
#             label=["SOCP" "LP"],
#             xlabel = "ESR #",
#             ylabel="Revenues",
#             title="EOD Revenues: ESRs",
#             legend=:outerright)
#
# p2_esr = groupedbar(["S1","S1","S2","S2","S3","S3"],
#         [esr_outcomes_SOCP[:Cos_esr] esr_outcomes_LP[:Cos_esr]],
#         label=["SOCP" "LP"],
#         xlabel = "ESR #",
#         ylabel="Costs",
#         title="EOD Costs: ESRs",
#         legend=:outerright)
#
# p3_esr = groupedbar(["S1","S1","S2","S2","S3","S3"],
#             [esr_outcomes_SOCP[:Pro_esr] esr_outcomes_LP[:Pro_esr]],
#             label=["SOCP" "LP"],
#             xlabel = "ESR #",
#             ylabel="Profits",
#             title="EOD Profits: ESRs",
#             legend=:outerright)
#
# plot(p1_esr,p2_esr,p3_esr, layout = (3,1))
# savefig("results/plots/cost_recovery_LP_SOCP.pdf")

""" REVENUE ADEQUACY """
function check_revenue_adequacy(settings,networkdata,winddata,loaddata,sol)
    """ Checking revenue adequacy of the market operator """
    # Yet to implement -> network operator's congestion rent formulation
    Ng = length(networkdata[:gens])
    Ns = length(networkdata[:esrs])
    Nw = size(winddata[:ŵ],1)

    payment_gens = zeros(Ng)
    payment_esrs = zeros(Ns)
    ener_in = 0
    for g in 1:Ng
        payment_gens[g] = sol[:Π][:Π_E][networkdata[:gens][g].node,:]'sol[:p][g,:] + sol[:Π][:Π_R][networkdata[:gens][g].node,:]'sol[:α][g,:]
    end
    for s in 1:Ns
        payment_esrs[s] = sol[:Π][:Π_E][networkdata[:esrs][s].node,:]'sol[:b][s,:] + sol[:Π][:Π_R][networkdata[:esrs][s].node,:]'sol[:γ][s,:]
    end
    payment_wpps    = tr(sol[:Π][:Π_E]'winddata[:Cwind]*winddata[:ŵ])
    #payment_loads   = tr(sol_sto_SOCP[:Π][:Π_E]'loaddata[:node_loads]) + tr(sol_sto_SOCP[:Π][:Π_R]'ones(length(networkdata[:buses]), settings[:T]))
    payment_loads   = tr(sol[:Π][:Π_E]'loaddata[:node_loads]) + sum(sol[:λ_R])
    surplus         = payment_loads - payment_wpps - sum(payment_gens) - sum(payment_esrs)

    mo_account = Dict(:agents => Dict(:gens => payment_gens,
                                      :esrs => payment_esrs,
                                      :wpps => payment_wpps,
                                      :loads=> payment_loads),
                      :commodity => Dict(:ener_in   => ener_in,
                                         :ener_out  => ener_out,
                                         :flex_in   => flex_in,
                                         :flex_out  => flex_out))
    return mo_account
end
