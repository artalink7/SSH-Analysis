function SingleVariance(pre_quench, post_quench, style, L_cells, N_cells, t; filling_fraction=1.0)
    eigs = EigenvaluesDensity_sub(pre_quench, post_quench, style, L_cells, N_cells, t; filling_fraction=filling_fraction) 
    eigs = eigs[(eigs .> 1e-14) .& (eigs .< 1 - 1e-14)] 
    return sum(eigs .*(1 .- eigs))  # variance formula
end

function GenerateDataVarianceEquilibrium(pre_quench, style, N_cells; max_LA = nothing)
    if max_LA === nothing
        max_LA = div(N_cells, 2)
    end

    LA_range = 1:max_LA
    variances_equilibrium = zeros(length(LA_range))

    Threads.@threads for i in eachindex(LA_range)
        LA = LA_range[i]
        variances_equilibrium[i] = SingleVariance(pre_quench, [0,1,1,0], style, LA, N_cells, 0.0)
    end
    return variances_equilibrium
end

function GenerateDataVariance(pre_quench, post_quench, style, L_cells, N_cells)
    time = 0:0.4:199.6
    nsteps = length(time)
    variance_quench = zeros(nsteps)

    @showprogress 1 for (i, t) in enumerate(time)
        variance_quench[i] = SingleVariance(pre_quench, post_quench, style, L_cells, N_cells, t)
    end
    return variance_quench
end
