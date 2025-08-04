function GenerateDataEntropy(pre_quench, post_quench, style, L_cells, N_cells)
    # pre_quench: parameters before the quench
    # post_quench: parameters after the quench
    # style: "discrete" or "continuos"
    #L_cells = Number of cells in the partition
    #N_cells = Number of cells in the system
    
    time = 0:0.4:199.6 
    entropy = zeros(length(time))
    
    for (i, t) in enumerate(time)
        Val1 = EigenvaluesDensity(pre_quench, post_quench, style, L_cells, N_cells, "1", t) .+ 0.5
        V_s = filter(x -> (x > 0 && x < 1.0), Val1)
        entropy[i] = -sum((1 .- V_s) .* log.(1 .- V_s) .+ V_s .*log.(V_s)) 
    end
    return entropy  
end

function GenerateDataEntropyEquilibrium(pre_quench, style, N_cells; max_LA = nothing)
    # pre_quench: parameters before the quench
    # L_cells: maximum number of cells in the partition

    if max_LA === nothing
        max_LA = div(N_cells, 2)
    end

    LA_range = 1:max_LA
    entropy_equilibrium = zeros(length(LA_range))
    
    Threads.@threads for i in eachindex(LA_range)
        LA = LA_range[i]
        Val1 = EigenvaluesDensity(pre_quench, [0, 1, 1, 0], style, LA, N_cells, "1", 0.0) .+ 0.5
        V_s = filter(x -> (x > 0 && x < 1.0), Val1)
        entropy_equilibrium[i] = -sum((1 .- V_s) .* log.(1 .- V_s) .+ V_s .*log.(V_s)) 
    end
    return entropy_equilibrium
end
