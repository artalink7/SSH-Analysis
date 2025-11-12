# --- Quenched 2x2 density matrix ---
function QuenchedDensityMatrix(pre_quench, post_quench, k, t)
    # pre_quench: parameters before the quench
    # post_quench: parameters after the quench
    # k: wavevector
    # t: Time
    σ_x = [0 1; 1 0]
    σ_y = [0 -im; im 0]
    σ_z = [1 0; 0 -1]

    n_x, n_y, n_z = BlochEvolved(k, pre_quench, post_quench, t)

    ρ = 0.5 * (I(2) - (n_x * σ_x + n_y * σ_y + n_z * σ_z))

    return ρ
end

# --- Correlation function in momentum space ---
function CorrelationMomentumSpace(pre_quench, post_quench, corr, r, k, t; occ_k= 1.0)
    # pre_quench: parameters before the quench
    # post_quench: parameters after the quench
    # corr: matrix entrance and type (real or imag)
    # r: distance
    # k: wavevector
    # t: time
    correlation, tipo = corr 
    ρ_k = exp(-im * k * r) * QuenchedDensityMatrix(pre_quench, post_quench, k, t) * occ_k
    ρ_flat = vec(ρ_k)
    ρ_00, ρ_01, ρ_10, ρ_11 = ρ_flat

    if correlation == "00"
        return tipo == "Real" ? real(ρ_00) : imag(ρ_00)
    elseif correlation == "01"
        return tipo == "Real" ? real(ρ_01) : imag(ρ_01)
    elseif correlation == "10"
        return tipo == "Real" ? real(ρ_10) : imag(ρ_10)
    elseif correlation == "11"
        return tipo == "Real" ? real(ρ_11) : imag(ρ_11)
    else
        error("Invalid correlation type. Use '00', '01', '10', or '11'.")
    end

end

function OccupiedKs(pre_quench, N_cells, filling_fraction)
    # pre_quench: parameters before the quench
    # N_cells: number of cells
    # filling_fraction: filling fraction (0 to 1)
    k_vals = 2π * (0:N_cells-1) ./ N_cells
    energies = Float64[]
    for k in k_vals
        h = BlochPostNonNorm(k, pre_quench[1], pre_quench[2], pre_quench[3], pre_quench[4])
        push!(energies, -norm(h))
    end
    idx_sorted = sortperm(energies)
    M = clamp(round(Int, filling_fraction * N_cells), 0, N_cells)
    occ = zeros(Float64, N_cells)
    for j in 1:M
        occ[idx_sorted[j]] = 1.0
    end
    return occ  
end

# --- Real space 2x2 density matrix ---
function RealSpaceDensityM(pre_quench, post_quench, r, style, N_cells, t; filling_fraction=1.0)
    # pre_quench: parameters before the quench
    # post_quench: parameters after the quench
    # r: distance
    # style: "discrete" or "continuos"
    # N_cells: number of cells
    # t: time
    if style == "continuos"
        function integrand(corr)
            k -> CorrelationMomentumSpace(pre_quench, post_quench, corr, r, k, t)
        end

        if r!= 0 
            ρ_00 = quadgk(integrand(("00", "Real")), -π, π)[1] / (2π) + im * quadgk(integrand(("00", "Imag")), -π, π)[1] / (2π)
            ρ_01 = quadgk(integrand(("01", "Real")), -π, π)[1] / (2π) + im * quadgk(integrand(("01", "Imag")), -π, π)[1] / (2π)
            ρ_10 = quadgk(integrand(("10", "Real")), -π, π)[1] / (2π) + im * quadgk(integrand(("10", "Imag")), -π, π)[1] / (2π)
        else
            ρ_00 = quadgk(integrand(("00", "Real")), -π, π)[1] / (2π) 
            real_01 = quadgk(integrand(("01", "Real")), -π, π)[1] / (2π)
            imag_01 = quadgk(integrand(("01", "Imag")), -π, π)[1] / (2π)
            ρ_01 = real_01 + im * imag_01
            ρ_10 = real_01 - im * imag_01
        end
    elseif style == "discrete"
        k_vals = 2π * (0:N_cells-1) ./ N_cells
        occ_vec = OccupiedKs(pre_quench, N_cells, filling_fraction)
        #println("N_cells: $N_cells, OccupiedKs: ", sum(occ_vec))

        function sumcorr(corr)
            sum(CorrelationMomentumSpace(pre_quench, post_quench, corr, r, k_vals[i], t; occ_k = occ_vec[i]) for i in eachindex(k_vals)) / N_cells
            
        end

        if r!=0
            ρ_00 = sumcorr(("00", "Real")) + im * sumcorr(("00", "Imag"))
            ρ_01 = sumcorr(("01", "Real")) + im * sumcorr(("01", "Imag"))
            ρ_10 = sumcorr(("10", "Real")) + im * sumcorr(("10", "Imag"))
            ρ_11 = sumcorr(("11", "Real")) + im * sumcorr(("11", "Imag"))
        else
            real_01 = sumcorr(("01", "Real"))
            imag_01 = sumcorr(("01", "Imag"))
            ρ_00 = sumcorr(("00", "Real"))
            ρ_01 = real_01 + im * imag_01
            ρ_10 = real_01 - im * imag_01
            ρ_11 = sumcorr(("11", "Real"))
        end
    end
    return [ρ_00 ρ_01; ρ_10 ρ_11]

end

# Build correlation matrix only for a subsystem of length L_sub (number of cells in the subsystem)
function DensityMatrix_subsystem(pre_quench, post_quench, style, L_sub, N_cells, t; filling_fraction=1.0)
    # pre_quench, post_quench: same as before
    # style: "discrete" / "continuos"
    # L_sub: number of cells in the subsystem (smaller than full system)
    # N_cells: number of k-points for discrete sums (kept for RealSpaceDensityM)
    # t: time
    # returns: 2L_sub × 2L_sub correlation matrix (ComplexF64)
    twoL = 2 * L_sub
    ρ_sub = zeros(ComplexF64, twoL, twoL)

    # r = 0..L_sub-1 correlations fill blocks along diagonals
    # block indices: for cell index a (1-based), rows 2a-1:2a and columns likewise
    # Fill all required r (distance) blocks
    for r in 0:(L_sub - 1)
        ρ_bloch = RealSpaceDensityM(pre_quench, post_quench, r, style, N_cells, t; filling_fraction=filling_fraction)
        # place ρ_bloch at positions (a, a+r) and Hermitian conjugate at (a+r, a)
        max_a = L_sub - r
        for a in 1:max_a
            row = 2*(a - 1) + 1
            col = 2*(a + r - 1) + 1
            # assign 2x2 block without kron
            @inbounds ρ_sub[row:row+1, col:col+1] .= ρ_bloch
            if r != 0
                # Hermitian counterpart (conjugate transpose)
                @inbounds ρ_sub[col:col+1, row:row+1] .= ρ_bloch'
            end
        end
    end

    return ρ_sub
end

# eigenvalues of subsystem correlation matrix, memory-efficient
function EigenvaluesDensity_sub(pre_quench, post_quench, style, L_sub, N_cells, t; filling_fraction=1.0)
    ρ_sub = DensityMatrix_subsystem(pre_quench, post_quench, style, L_sub, N_cells, t; filling_fraction=filling_fraction)
    # Ensure Hermitian wrapper to use specialized LAPACK
    return eigvals(Hermitian(ρ_sub))
end

# and EntropyandVariance for subsystem
function EntropyandVariance_sub(pre_quench, post_quench, style, L_sub, N_cells, t; filling_fraction=1.0)
    eigs = EigenvaluesDensity_sub(pre_quench, post_quench, style, L_sub, N_cells, t; filling_fraction=filling_fraction)
    eigs = eigs[(eigs .> 1e-14) .& (eigs .< 1 - 1e-14)]
    S = -sum(eigs .* log.(eigs) .+ (1 .- eigs) .* log.(1 .- eigs))
    V = sum(eigs .* (1 .- eigs))
    return S, V
end
