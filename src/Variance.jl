function SingleVariance(pre_quench, post_quench, style, L_cells, N_cells, cut, t)
    ρ = DensityMatrix(pre_quench, post_quench, style, L_cells, N_cells, cut, t)
    ρ = ρ + 0.5.*I  
    norm_squared = sum(abs2, ρ)  # sum of squared absolute values of all elements
    return L_cells - norm_squared
end

function GenerateDataVarianceEquilibrium(pre_quench, style, N_cells; max_LA = nothing)
    if max_LA === nothing
        max_LA = div(N_cells, 2)
    end

    LA_range = 1:max_LA
    variances_equilibrium = zeros(length(LA_range))

    Threads.@threads for i in eachindex(LA_range)
        LA = LA_range[i]
        val = SingleVariance(pre_quench, [0,1,1,0], style, LA, N_cells, "1", 0.0)
        variances_equilibrium[i] = val
    end
    return variances_equilibrium
end

function GenerateDataVariance(pre_quench, post_quench, style, L_cells, N_cells)
    # Define time array
    time = 0:0.4:199.6
    nsteps = length(time)

    # Preallocate arrays
    variance_quench = zeros(nsteps)

    # Main loop
    @showprogress 1 for (i, t) in enumerate(time)
        ρ = DensityMatrix(pre_quench, post_quench, style, L_cells, N_cells, "1", t)
        ρ = ρ + 0.5.*I 
        variance_quench[i] = L_cells - sum(abs2, ρ)
    end
    return variance_quench
end
