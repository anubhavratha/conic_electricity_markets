using CSV, DataFrames
using DataStructures: SortedDict

mutable struct Generator
   ind::Int
   c::Any
   p̅::Any
   p̲::Any
   r̅::Any
   α̅::Any
   α̲::Any
   Ra̅::Any
   c_r::Any
   node::Int
   zone::Any
   q::Any
   function Generator(ind,c,p̅,p̲,r̅,α̅,α̲,Ra̅,c_r,node,zone,q)
      i = new()
      i.ind  = ind
      i.c = c
      i.p̅ = p̅
      i.p̲ = p̲
      i.r̅ = r̅
      i.α̅ = α̅
      i.α̲ = α̲
      i.Ra̅ = Ra̅
      i.c_r = c_r
      i.node = node
      i.zone = zone
      i.q = q
      return i
   end
end
mutable struct Bus
   ind::Int
   ls_pc::Float64
   w_cap::Float64
   σ2::Float64
   zone::Any
   function Bus(ind,ls_pc,w_cap,σ2,zone)
      i = new()
      i.ind  = ind
      i.ls_pc = ls_pc
      i.w_cap = w_cap
      i.σ2 = σ2
      i.zone = zone
      return i
   end
end
mutable struct Line
   ind::Int
   b_f::Int
   b_t::Int
   β::Float64
   f̅::Float64
   function Line(ind,b_f,b_t,β,f̅)
      i = new()
      i.ind  = ind
      i.b_f = b_f
      i.b_t = b_t
      i.β = β
      i.f̅ = f̅
      return i
   end
end
mutable struct EnergyStorageResource
    ind::Int
    node::Int
    E̅::Float64
    E𝟶::Float64
    c_q::Float64
    c_l::Float64
    c_r::Float64
    P̅c::Float64
    P̅d::Float64
    s̅::Float64
    s̲::Float64
    η_c::Float64
    η_d::Float64
    B_s::Float64
    γ̅::Int
    γ̲::Int
    function EnergyStorageResource(ind,node,E̅,E𝟶,c_q,c_l,c_r,P̅c,P̅d,s̅,s̲,η_c,η_d,B_s,γ̅,γ̲)
        i=new()
        i.ind = ind
        i.node = node
        i.E̅ = E̅
        i.c_q = c_q
        i.c_l = c_l
        i.c_r = c_r
        i.P̅c = P̅c
        i.P̅d = P̅d
        i.s̅ = s̅
        i.s̲ = s̲
        i.η_c = η_c
        i.η_d = η_d
        i.E𝟶 = E𝟶
        i.B_s = B_s
        i.γ̅ = γ̅
        i.γ̲ = γ̲
        return i
    end
end
function load_data(settings)
    gen_data    = DataFrame(CSV.File("data/generators.csv"))
    bus_data    = DataFrame(CSV.File("data/buses.csv"))
    if !settings[:bottleneck]
        line_data   = DataFrame(CSV.File("data/lines.csv"))
    else
        line_data   = DataFrame(CSV.File("data/lines_congested.csv"))
    end
    esr_data    = DataFrame(CSV.File("data/esr_data.csv"))

    gen = Dict()
    for i in 1:nrow(gen_data)
        ind     = gen_data[i, :index]
        c       = gen_data[i, :C]
        p̅       = gen_data[i, :P_max]
        p̲       = gen_data[i, :P_min]
        r̅       = gen_data[i, :R_max]
        α̅       = 1
        α̲       = -1
        Ra̅      = gen_data[i, :Ramp_max]
        c_r     = gen_data[i, :C_r]
        node    = gen_data[i, :node]
        zone    = gen_data[i, :zone]
        q       = gen_data[i, :Q]
        add_generator = Generator(ind,c,p̅,p̲,r̅,α̅,α̲,Ra̅,c_r,node,zone,q)
        gen[add_generator.ind] = add_generator
    end

    esr=Dict()
    if settings[:ESRs] != false
        for i in 1:nrow(esr_data)
            ind     = esr_data[i, :index]
            node    = esr_data[i, :Bus]
            E̅       = esr_data[i, :E_max]
            c_q     = esr_data[i, :Qcost]
            c_l     = esr_data[i, :Lcost]
            c_r     = esr_data[i, :C_r]
            P̅c      = esr_data[i, :PmaxC]
            P̅d      = esr_data[i, :PmaxD]
            s̅       = esr_data[i, :SOCMax]*E̅
            s̲       = esr_data[i, :SOCMin]*E̅
            η_c     = esr_data[i, :EffC]
            η_d     = esr_data[i, :EffD]
            E𝟶      = esr_data[i, :SOCinit]*E̅
            B_s     = esr_data[i, :SOCBound]*E̅
            γ̅       = 1
            γ̲       = -1
            #ESR Capacity Override by Wind Penetration, 6 wind farms and 3 ESRs
            if settings[:ESRoverride]
                P̅c  = settings[:Wcap]*0.25
                P̅d  = settings[:Wcap]*0.25
                E̅   = P̅c * 10
                s̅   = esr_data[i, :SOCMax]*E̅
                s̲   = esr_data[i, :SOCMin]*E̅
                E𝟶  = esr_data[i, :SOCinit]*E̅
                B_s = esr_data[i, :SOCBound]*E̅
            end
            add_storage = EnergyStorageResource(ind,node,E̅,E𝟶,c_q,c_l,c_r,P̅c,P̅d,s̅,s̲,η_c,η_d,B_s,γ̅,γ̲)
            esr[add_storage.ind] = add_storage
        end
    end

    bus = Dict()
    for i in 1:nrow(bus_data)
        ind = bus_data[i, :node]
        ls_pc = bus_data[i, :loadshare_pc]
        w_cap = bus_data[i, :w_cap]
        σ2 = 0
        zone = bus_data[i, :zone]
        add_bus = Bus(ind,ls_pc,w_cap,σ2,zone)
        bus[add_bus.ind] = add_bus
    end

    line = Dict()
    for i in 1:nrow(line_data)
        ind = line_data[i, :index]
        b_f = line_data[i, :b_f]
        b_t = line_data[i, :b_t]
        β   = 100/line_data[i, :x_ij]
        f̅   = line_data[i, :f_max]
        add_line = Line(ind,b_f,b_t,β,f̅)
        line[add_line.ind] = add_line
    end

    #sort dictionaries
    bus     = SortedDict(bus)
    gen     = SortedDict(gen)
    line    = SortedDict(line)
    esr     = SortedDict(esr)
    Nb      = length(bus)
    Ng      = length(gen)
    Ns      = length(esr)
    line_set= collect(keys(line))
    B = zeros(Nb,Nb)
    for i in 1:Nb
        for l in line_set
            if line[l].b_f == i || line[l].b_t  == i
                B[i,i] += line[l].β
            end
        end
    end
    for l in line_set
        B[line[l].b_f,line[l].b_t] = -line[l].β
        B[line[l].b_t,line[l].b_f] = -line[l].β
    end
    #Max flow limits - if a line exists between two buses, it takes the value of the flow limit, else 0
    f̅ = zeros(Nb,Nb)
    for l in line_set
        f̅[line[l].b_f,line[l].b_t] = line[l].f̅
        f̅[line[l].b_t,line[l].b_f] = line[l].f̅
    end

    β = zeros(Nb,Nb)
    for l in line_set
        β[line[l].b_f,line[l].b_t] = line[l].β
        β[line[l].b_t,line[l].b_f] = line[l].β
    end

    function remove_col_and_row(B,refbus)
        @assert size(B,1) == size(B,2)
        n = size(B,1)
        return B[1:n .!= refbus, 1:n .!= refbus]
    end

    function build_B̆(B̂inv,refbus)
        Nb = size(B̂inv,1)+1
        B̆ = zeros(Nb,Nb)
        for i in 1:Nb, j in 1:Nb
            if i < refbus && j < refbus
                B̆[i,j] = B̂inv[i,j]
            end
            if i > refbus && j > refbus
                B̆[i,j] = B̂inv[i-1,j-1]
            end
            if i > refbus && j < refbus
                B̆[i,j] = B̂inv[i-1,j]
            end
            if i < refbus && j > refbus
                B̆[i,j] = B̂inv[i,j-1]
            end
        end
        return B̆
    end

    refbus = 24
    B̂=remove_col_and_row(B,refbus)
    B̂inv = inv(B̂)
    𝛑=build_B̆(B̂inv,refbus)

    #Make Flow Susceptance Matrix
    numElLines = length(line_set)
    Bflow = zeros(numElLines,Nb)
    for l in 1:numElLines
        Bflow[l, line[l].b_f] =  line[l].β
        Bflow[l, line[l].b_t] = -line[l].β
    end

    #Forming the PTDF Matrix
    PTDF = Bflow*𝛑

    #==== Network Data Pre-processing for PTDF Formulation =====#
    Cgens = zeros(Nb,Ng)       #Generator matrix
    for i=1:Ng
        ThisGenBusNum = gen[i].node
        Cgens[ThisGenBusNum,i] = 1
    end
    Cesrs = zeros(Nb, Ns)
    for i=1:Ns
        ThisESRBusNum = esr[i].node
        Cesrs[ThisESRBusNum, i] = 1
    end
    PL = zeros(length(line),1)      #Vector of Line Flow Limits
    for l in 1:length(line)
        PL[l,1] = line_data[l, :f_max]
    end

    # Quantifying problem dimensions
    num_constraints  = (settings[:T]*length(gen)*4 +          #Flex gens reserve limits and production Limits
                            (settings[:T]-1)*length(gen)*2 +     #Flex gens ramping limits
                                settings[:T]*length(line)*2 +    #Line flow limits
                                    settings[:T]*length(esr)*4 + #Storage operation power limits
                                        2*length(esr))           #End of day storage bounds
    num_variables    = settings[:T]*(length(gen) + length(esr))*2

    networkdata = Dict(
    :gens   => gen,
    :esrs   => esr,
    :buses  => bus,
    :lines  => line,
    :B      => B,
    :f̅      => PL,
    :𝛑      => 𝛑,
    :ref    => refbus,
    :Ψ      => PTDF,
    :Cgens  => Cgens,
    :Cesr   => Cesrs,
    :ncon   => num_constraints,
    :nvar   => num_variables
    )
    return networkdata
end
