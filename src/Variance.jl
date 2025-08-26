function SingleVariance(pre_quench, post_quench, style, L_cells, N_cells, cut, t)
    ρ = DensityMatrix(pre_quench, post_quench, style, L_cells, N_cells, cut, t)
    C_A = ρ + 0.5 .* I  # recover correlation matrix
    return real(tr(C_A * (I - C_A)))  # variance formula
end

function GenerateDataVarianceEquilibrium(pre_quench, style, N_cells; max_LA = nothing)
    if max_LA === nothing
        max_LA = div(N_cells, 2)
    end

    LA_range = 1:max_LA
    variances_equilibrium = zeros(length(LA_range))

    Threads.@threads for i in eachindex(LA_range)
        LA = LA_range[i]
        variances_equilibrium[i] = SingleVariance(pre_quench, [0,1,1,0], style, LA, N_cells, "1", 0.0)
    end
    return variances_equilibrium
end

function GenerateDataVariance(pre_quench, post_quench, style, L_cells, N_cells)
    time = 0:0.4:199.6
    nsteps = length(time)
    variance_quench = zeros(nsteps)

    @showprogress 1 for (i, t) in enumerate(time)
        variance_quench[i] = SingleVariance(pre_quench, post_quench, style, L_cells, N_cells, "1", t)
    end
    return variance_quench
end
