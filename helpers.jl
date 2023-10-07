#auxiliary functions to return Cholesky decomposition matrices of various dimensions for CC reformulation

## Cholesky approach
returnblks(Σ,t,Nw)          = Matrix(cholesky(Σ[(t-1)*Nw+1 : t*Nw, (t-1)*Nw+1 : t*Nw]).L)
returnblksramping(Σ,t,Nw)   = Matrix(cholesky(Σ[(t-2)*Nw+1 : t*Nw, (t-2)*Nw+1 : t*Nw]).L)
returnblksstorage(Σ,t,Nw)   = Matrix(cholesky(Σ[1:t*Nw, 1:t*Nw]).L)

## Square root approach
# returnblks(Σ,t,Nw)          = sqrt((Σ[(t-1)*Nw+1 : t*Nw, (t-1)*Nw+1 : t*Nw]))
# returnblksramping(Σ,t,Nw)   = sqrt(Σ[(t-2)*Nw+1 : t*Nw, (t-2)*Nw+1 : t*Nw])
# returnblksstorage(Σ,t,Nw)   = sqrt(Σ[1:t*Nw, 1:t*Nw])

#auxiliary functions to give sending and receiving nodes for a given line
ns(l) = networkdata[:lines][l].b_f
nr(l) = networkdata[:lines][l].b_t

function Φ(settings,networkdata)
    """Risk parameter r_ε"""
    settings[:det] == true ? res = 0 : NaN
    (settings[:det] == false && settings[:drcc] == false) ? res = quantile(Normal(0,1), 1 - settings[:ε]/networkdata[:ncon]) : NaN
    (settings[:det] == false && settings[:drcc] == true)  ? res = sqrt((1 - settings[:ε]/networkdata[:ncon])/(settings[:ε]/networkdata[:ncon])) : NaN
    return res
end

function compute_num_scenarios(settings,networkdata)
    """Number of scenarios for probabilistic guarantee by the scenario program"""
    if settings[:S] == 0
        β = 0.99
        S = Int(ceil((settings[:ε])^(-1)(MathConstants.e/(MathConstants.e - 1))*(networkdata[:nvar] - 1 + log(1/β))))
    else
        S = settings[:S]
    end
    return S
end

function prep_insample_data(settings,ŵ)
    """Wind forecast data for CC reformulation and "in-sample" scenarios"""
    T = settings[:T]          #time periods
    S = settings[:S]                   #number of samples
    Nw = size(ŵ,1)    #number of uncertainty sources
    Σ = zeros(settings[:T]*Nw, settings[:T]*Nw) #large spatio temporal covariance matrix
    ξ = zeros(Nw,S,settings[:T])            #scenarios data
    for t in 1:T
        σ_t = zeros(Nw)             #vector of standard deviations
        c_t  = 0.1                 #standard deviation vector
        C_t  = zeros(Nw,Nw)         #correlation matrix
        for k in 1:Nw,j in 1:Nw
            σ_t[k] = settings[:σ]*ŵ[k,1]
            k != j ? C_t[k,j] = c_t : NaN
            k == j ? C_t[k,j] = 1 : NaN
        end
        Σ_t   = cor2cov(C_t, σ_t)
        Σ½_t  = cholesky(Σ_t).L
        ξ[1:Nw,:,t] = rand(MvNormal(zeros(Nw),Σ_t),S)
        Σ[(t-1)*Nw+1 : t*Nw, (t-1)*Nw+1 : t*Nw] = Σ_t
    end
    forecastdata = Dict(:ξ => ξ, :Σ => Σ, :S => S)
    return forecastdata
end

function prep_wind_load_data(settings,networkdata)
    """ Wind forecast and load data"""
    Wcap        = settings[:Wcap]     #MW
    wind_buses  = [3,5,7,15,21,23]
    ŵ           = (Wcap/100) .* Matrix(DataFrame(CSV.File("data/Point_Forecast_Values.csv", header=false)))[:,1:settings[:T]]                 #Point forecast
    if settings[:histdata] == true
        Σ       = Matrix(DataFrame(CSV.File("data/Covariance_Matrix_Data.csv", header=false))) #Large Spatial Temporal Covariance Matrix
        ξ       = 0
        S       = forecastdata[:S]
    else
        forecastdata = prep_insample_data(settings,ŵ)
        Σ       = forecastdata[:Σ]
        ξ       = forecastdata[:ξ]
        S       = forecastdata[:S]
    end
    Cwind = zeros(length(networkdata[:buses]),length(wind_buses))       #Wind PP Matrix
    for i=1:length(wind_buses)
        ThisWindPPBusNum = wind_buses[i]
        Cwind[ThisWindPPBusNum,i] = 1
    end
    winddata = Dict(:Σ => Σ, :ŵ => ŵ, :Wcap => Wcap, :Cwind => Cwind, :ξ => ξ, :S => S)

    load_buses = []
    for n in 1:length(networkdata[:buses])
        if networkdata[:buses][n].ls_pc > 0
            push!(load_buses,n)
        end
    end
    Cload = zeros(length(networkdata[:buses]), length(load_buses))   #load matrix
    for n in 1:length(load_buses)
        ThisLoadBusNum = load_buses[n]
        Cload[ThisLoadBusNum, n] = 1
    end
    load_profile = [1775.835,   1669.815,   1590.300,   1563.795,
                    1563.795,   1590.300,   1961.370,   2279.430,
                    2517.975,   2544.480,   2544.480,   2517.975,
                    2517.975,   2517.975,   2464.965,   2464.965,
                    2623.995,   2650.500,   2650.500,   2544.480,
                    2411.955,   2199.915,   1934.865,   1669.815]
    #demand_scaling_factor
    node_loads = zeros(length(networkdata[:buses]),settings[:T])
    for b in load_buses
        node_loads[b,:] = 0.01*networkdata[:buses][b].ls_pc * load_profile[1:settings[:T]]'
    end
    loaddata = Dict(:D => load_profile, :Cload => Cload, :node_loads => node_loads)
    return winddata,loaddata
end

function prep_outofsample_data(settings,winddata)
    T       = settings[:T]             #time periods
    ST      = settings[:TDsize]       #number of samples for test data
    Nw      = size(winddata[:ŵ],1)    #number of uncertainty sources
    σ_oos   = settings[:σ_oos]            #bias of covariance matrix of real-distribution
    ξ_oos   = zeros(Nw,ST,settings[:T])            #scenarios data
    for t in 1:T
        σ_t = zeros(Nw)             #vector of standard deviations
        c_t  = 0.1                 #standard deviation vector
        C_t  = zeros(Nw,Nw)         #correlation matrix
        for k in 1:Nw,j in 1:Nw
            σ_t[k] = σ_oos*winddata[:ŵ][k,1]
            k != j ? C_t[k,j] = c_t : NaN
            k == j ? C_t[k,j] = 1 : NaN
        end
        Σ_t   = cor2cov(C_t, σ_t)
        ξ_oos[1:Nw,:,t] = rand(MvNormal(zeros(Nw),Σ_t),ST)
        # for t in 1:T
        #     for k in 1:Nw
        #         ξ_oos[k,:,t] = 0.25*winddata[:ŵ][k,t] .* (-0.5 .+ rand(Beta(2,5),(1,ST)))
        #     end
        # end
    end
    oosdata = Dict(:ξ_oos => ξ_oos, :S => ST)
    return oosdata
end


function pwl_cost_approx(settings,networkdata,g,s)
    """Piecewise linear approximation of quadratic costs"""
    cost_coeffs = Dict()
    if s == 0
        pvals  = range(networkdata[:gens][g].p̲, networkdata[:gens][g].p̅, length=settings[:Nsteps]+1)
        cvals  = (networkdata[:gens][g].c) .* pvals + networkdata[:gens][g].q .* (pvals .^ 2)
        cpmin  = (networkdata[:gens][g].c) .* networkdata[:gens][g].p̲ + networkdata[:gens][g].q .* (networkdata[:gens][g].p̲ .^ 2)
        slopes = []
        for i in 1:length(cvals)-1
            slopes = push!(slopes, (cvals[i+1] - cvals[i]) / (pvals[i+1] - pvals[i]))
        end
        cost_coeffs[:slopes] = slopes
        cost_coeffs[:pvals]  = pvals
        cost_coeffs[:cvals]  = cvals
        cost_coeffs[:cpmin]  = cpmin
    elseif g == 0
        pvals = range(-networkdata[:esrs][s].P̅c, networkdata[:esrs][s].P̅d, length=settings[:Nsteps]+1)
        cvals  = (networkdata[:esrs][s].c_l) .* pvals + networkdata[:esrs][s].c_q .* (pvals .^ 2)
        cpmin  = (networkdata[:esrs][s].c_l) .* (-networkdata[:esrs][s].P̅c) + networkdata[:esrs][s].c_q .* ((-networkdata[:esrs][s].P̅c) .^ 2)
        slopes = []
        for i in 1:length(cvals)-1
            slopes = push!(slopes, (cvals[i+1] - cvals[i]) / (pvals[i+1] - pvals[i]))
        end
        cost_coeffs[:slopes] = slopes
        cost_coeffs[:pvals]  = pvals
        cost_coeffs[:cvals]  = cvals
        cost_coeffs[:cpmin]  = cpmin
    else
        @warn("something is wrong with piecewise linearization of costs")
    end
    return cost_coeffs
end


function kernel_unbiased_2d(xvar,yvar,xgrid,ygrid)
    xobs = vec([xvar yvar])
    xpred = vec([xgrid ygrid])
    dim_ker = length(xobs)
    σ_k = 0.1
    H_matrix = diagm(σ_k .* ones(dim_ker))
    return (2π)^(-dim_ker/2) * exp(-0.5 * (xobs - xpred)'*inv(H_matrix)*(xobs - xpred))
end

function kernel_biased_2d(xvar,yvar,ntype,xgrid,ygrid)
    xobs = vec([xvar yvar])
    xpred = vec([xgrid ygrid])
    dim_ker = length(xobs)
    σ_k = 2
    if     ntype == 1 #nodes along a horizontal pipeline
        H_matrix = diagm([σ_k, 0.5*σ_k])
    elseif ntype == 2 #nodes along a vertical pipeline
        H_matrix = diagm([0.5*σ_k, σ_k])
    elseif ntype == 3 #nodes along a acute angle w.r.t horizontal
        H_matrix = [0.4 σ_k; σ_k 0.4]
    elseif ntype == 4 # nodes along a obtuse angle w.r.t horizontal
        H_matrix = [σ_k 0.05; 0.05 2*σ_k]
    else
        H_matrix = diagm([σ_k, σ_k])
    end
    return (2π)^(-dim_ker/2) * exp(-0.5 * (xobs - xpred)'*inv(H_matrix)*(xobs - xpred))
end



function prepare_kde_data(kde_settings)
    rawdata=DataFrame(CSV.File("results/R_files/out_data_24node.csv"))

    xvar = rawdata[!, :x]
    yvar = rawdata[!, :y]
    ntype = rawdata[!, :ntype]
    fval = rawdata[!, :el_price]

    gridvals=kde_settings[:gridsize]
    xgrid = Array(0:gridvals:18)
    ygrid = Array(0:gridvals:30)

    kerData = DataFrame(xval=Any[],yval=Any[],fgval=Any[])
    fgval=zeros(length(ygrid))
    for gx in 1:length(xgrid)
        for gy in 1:length(ygrid)
            fgval_temp = zeros(length(xvar))
            ker=zeros(length(xvar))
            for d in 1:length(xvar)
                if kde_settings[:bias] == 1
                    ker[d] = kernel_biased_2d(xvar[d],yvar[d],ntype[d],xgrid[gx],ygrid[gy])
                elseif kde_settings[:bias] == 0
                    ker[d] = kernel_unbiased_2d(xvar[d],yvar[d],xgrid[gx],ygrid[gy])
                end
                fgval_temp[d] = fval[d] * ker[d]
            end
            fgval[gy] = sum(fgval_temp)
            push!(kerData,[xgrid[gx],ygrid[gy],fgval[gy]])
        end
    end
    return kerData
end
