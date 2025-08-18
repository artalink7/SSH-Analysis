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

    ρ = -0.5 * (n_x * σ_x + n_y * σ_y + n_z * σ_z)

    return ρ
end

# --- Correlation function in momentum space ---
function CorrelationMomentumSpace(pre_quench, post_quench, corr, r, k, t)
    # pre_quench: parameters before the quench
    # post_quench: parameters after the quench
    # corr: matrix entrance and type (real or imag)
    # r: distance
    # k: wavevector
    # t: time
    correlation, tipo = corr 
    ρ_k = exp(-im * k * r) * QuenchedDensityMatrix(pre_quench, post_quench, k, t)
    ρ_flat = vec(ρ_k)
    ρ_00, ρ_01, ρ_10, ρ_11 = ρ_flat

    if correlation == "00"
        return tipo == "Real" ? real(ρ_00) : imag(ρ_00)
    elseif correlation == "01"
        return tipo == "Real" ? real(ρ_01) : imag(ρ_01)
    elseif correlation == "10"
        return tipo == "Real" ? real(ρ_10) : imag(ρ_10)
    end

end

# --- Real space 2x2 density matrix ---
function RealSpaceDensityM(pre_quench, post_quench, r, style, N_cells, t)
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

        function sumcorr(corr)
            sum(CorrelationMomentumSpace(pre_quench, post_quench, corr, r, k, t) for k in k_vals) / N_cells
            
        end

        if r!=0
            ρ_00 = sumcorr(("00", "Real")) + im * sumcorr(("00", "Imag"))
            ρ_01 = sumcorr(("01", "Real")) + im * sumcorr(("01", "Imag"))
            ρ_10 = sumcorr(("10", "Real")) + im * sumcorr(("10", "Imag"))
        else
            real_01 = sumcorr(("01", "Real"))
            imag_01 = sumcorr(("01", "Imag"))
            ρ_00 = sumcorr(("00", "Real"))
            ρ_01 = real_01 + im * imag_01
            ρ_10 = real_01 - im * imag_01
        end
    end

    ρ_11 = -ρ_00
    return [ρ_00 ρ_01; ρ_10 ρ_11]

end

# --- Full many-body density matrix ---
function DensityMatrix(pre_quench, post_quench, style, L_cells, N_cells, cut, t)
    # pre_quench: parameters before the quench
    # post_quench: parameters after the quench
    # style: "discrete" or "continuos"
    # cut: cut position
    # t: time
    
    if style == "continuos"
        N_cells = L_cells
    end

    R = 1:L_cells - 1
    ρ_zero = zeros(ComplexF64, 2L_cells, 2L_cells)
    ρ_base = zeros(ComplexF64, 2L_cells, 2L_cells)

    ρ_diag = RealSpaceDensityM(pre_quench, post_quench, 0, style, N_cells, t)
    ρ_zero .= kron(I(L_cells), ρ_diag)
    
    for r in R 
        ρ_off = zeros(ComplexF64, 2L_cells, 2L_cells)
        ρ_bloch = RealSpaceDensityM(pre_quench, post_quench, r, style, N_cells, t)
        ρ_id = I(L_cells - r)
        ρ_off[1:end-2r, 2r+1:end] .= kron(ρ_id, ρ_bloch)
        ρ_base .+= ρ_off
    end

    ρ_total = ρ_base + ρ_base' + ρ_zero

    if cut == "1"
        return ρ_total
    elseif cut == "2"
        return ρ_total[1:end-1, 1:end-1]
    elseif cut == "3"
        return ρ_total[2:end-1, 2:end-1]
    end        

end

# --- Eigenvalues of the density matrix ---
function EigenvaluesDensity(pre_quench, post_quench, style, L_cells, N_cells, cut, t)
    # pre_quench: parameters before the quench
    # post_quench: parameters after the quench
    # style: "discrete" or "continuos"  

    ρ = DensityMatrix(pre_quench, post_quench, style, L_cells, N_cells, cut, t)
    return eigvals(ρ)
end
