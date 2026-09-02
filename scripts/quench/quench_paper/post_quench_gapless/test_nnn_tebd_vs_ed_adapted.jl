push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis
using ITensors
using LinearAlgebra

# ============================================================
# Minimal test setup
#
# This is deliberately kept single-threaded: N=3..5 is tiny, and
# we want to eliminate threading as a possible source of confusion.
# ============================================================

using ITensorMPS
using Printf

ITensors.disable_threaded_blocksparse()
BLAS.set_num_threads(1)

println("Julia threads: $(Threads.nthreads()), BLAS threads: $(BLAS.get_num_threads())")
println("ITensors version: 0.6.23 expected")
println("ITensorMPS version: 0.2.6 expected")


# ============================================================
# Your model
# ============================================================

struct ExtendedTBParams
    t0::Float64
    t2::Float64
    V::Float64
end

ExtendedTBParams(; t0=1.0, t2=0.0, V=0.0) =
    ExtendedTBParams(t0, t2, V)


function product_state_Nf(N::Int, Nf::Int; pattern::Symbol=:alternating)
    @assert 0 <= Nf <= N

    st = fill("0", N)

    if pattern == :alternating
        order = vcat(collect(1:2:N), collect(2:2:N))
    else
        error("Unknown pattern = $pattern")
    end

    for k in 1:Nf
        st[order[k]] = "1"
    end

    return st
end


# ============================================================
# EXACT FERMIONIC OPERATORS
#
# Basis convention:
#
#   |n1 n2 ... nN>
#
# where site 1 is the least-significant bit.
#
# This is the standard Jordan-Wigner convention:
#
# c†_i |n> = (-1)^(sum_{j<i} n_j)
#            |... 1_i ...>
# ============================================================

function exact_cdag(N::Int, i::Int)
    dim = 2^N
    Cdag = zeros(ComplexF64, dim, dim)

    for state in 0:(dim - 1)

        ni = (state >> (i - 1)) & 1

        # c† cannot act on an occupied site
        if ni == 1
            continue
        end

        # Jordan-Wigner parity to the left of site i
        parity = 0
        for j in 1:(i - 1)
            parity += (state >> (j - 1)) & 1
        end

        sign = parity % 2 == 0 ? 1.0 : -1.0

        newstate = state | (1 << (i - 1))

        Cdag[newstate + 1, state + 1] = sign
    end

    return Cdag
end


function exact_c(N::Int, i::Int)
    return exact_cdag(N, i)'
end


function exact_n(N::Int, i::Int)
    cdag = exact_cdag(N, i)
    c = cdag'
    return cdag * c
end


# ============================================================
# EXACT HAMILTONIAN
# ============================================================

function exact_extended_H(N::Int, p::ExtendedTBParams)

    dim = 2^N
    H = zeros(ComplexF64, dim, dim)

    # --------------------------------------------------------
    # NN hopping + NN interaction
    # --------------------------------------------------------

    for i in 1:(N - 1)

        ci_dag = exact_cdag(N, i)
        ci     = exact_c(N, i)

        cj_dag = exact_cdag(N, i + 1)
        cj     = exact_c(N, i + 1)

        # NN hopping
        H += -p.t0 * (ci_dag * cj + cj_dag * ci)

        # NN interaction:
        #
        # V n_i n_{i+1}
        # - V/2 n_i
        # - V/2 n_{i+1}
        #
        if p.V != 0.0

            ni = exact_n(N, i)
            nj = exact_n(N, i + 1)

            H += p.V * ni * nj
            H += -(p.V / 2) * ni
            H += -(p.V / 2) * nj
        end
    end

    # --------------------------------------------------------
    # NNN hopping
    # --------------------------------------------------------

    for i in 1:(N - 2)

        ci_dag = exact_cdag(N, i)
        ci     = exact_c(N, i)

        cj_dag = exact_cdag(N, i + 2)
        cj     = exact_c(N, i + 2)

        H += -p.t2 * (ci_dag * cj + cj_dag * ci)
    end

    return H
end



# ============================================================
# PURE NNN TEBD GATES
#
# This test deliberately contains ONLY
#
#     H_NNN = -t2 * sum_i (c†_i c_{i+2} + h.c.)
#
# so that we isolate the Jordan-Wigner / AutoFermion issue.
#
# include_F=false:
#     no explicit F operator
#
# include_F=true:
#     the original explicit-F construction
# ============================================================

function build_nnn_test_gates(
    sites;
    t2=1.0,
    dt=1e-4,
    include_F=false,
)

    N = length(sites)
    gates = ITensor[]

    for i in 1:(N - 2)

        s1 = sites[i]
        s2 = sites[i + 1]
        s3 = sites[i + 2]

        Id2 = op("Id", s2)
        Id3 = op("Id", s3)

        # IMPORTANT:
        # Do NOT initialize h as ITensor(0).
        # That is a rank-0 scalar and cannot be added to a
        # three-site/rank-6 tensor.
        h = zero(op("Id", s1) * Id2 * Id3)

        if include_F

            # ORIGINAL construction
            h += -t2 *
                op("Cdag", s1) *
                op("F", s2) *
                op("C", s3)

            h += -t2 *
                op("Cdag", s3) *
                op("F", s2) *
                op("C", s1)

        else

            # TESTED CANDIDATE: no explicit Jordan-Wigner F
            h += -t2 *
                op("Cdag", s1) *
                Id2 *
                op("C", s3)

            h += -t2 *
                op("Cdag", s3) *
                Id2 *
                op("C", s1)
        end

        push!(
            gates,
            exp(-1im * dt / 2 * h),
        )
    end

    # Symmetric second-order sweep
    append!(gates, reverse(gates))

    return gates
end


# ============================================================
# MPS -> DENSE VECTOR
#
# The exact ED code below uses |n1 n2 ... nN> with site 1 as the
# least-significant bit. vec(Array(...)) uses the first array
# index as the fastest-moving index, so this convention matches.
# ============================================================

function mps_to_dense(psi, sites)

    T = psi[1]

    for j in 2:length(psi)
        T *= psi[j]
    end

    T = permute(T, sites...)

    A = Array(T, sites...)

    return ComplexF64.(vec(A))
end


# ============================================================
# EXACT PRODUCT STATE
# ============================================================

function exact_product_vector(states::Vector{String})

    N = length(states)
    dim = 2^N

    v = zeros(ComplexF64, dim)

    state = 0

    for i in 1:N
        if states[i] == "1"
            state |= 1 << (i - 1)
        elseif states[i] == "0"
            nothing
        else
            error("Unknown state $(states[i])")
        end
    end

    v[state + 1] = 1.0

    return v
end


# ============================================================
# ONE PURE-NNN ED VS TEBD COMPARISON
# ============================================================

function compare_nnn_tebd_ed(
    N::Int;
    dt=1e-4,
    t2=1.0,
    Nf=N ÷ 2,
    include_F=false,
)

    p = ExtendedTBParams(
        t0=0.0,
        t2=t2,
        V=0.0,
    )

    states = product_state_Nf(N, Nf)

    # -------------------------
    # Exact evolution
    # -------------------------

    psi0_exact = exact_product_vector(states)

    H = exact_extended_H(N, p)

    psi_exact =
        exp(-1im * dt * H) *
        psi0_exact

    # -------------------------
    # TEBD evolution
    # -------------------------

    sites = siteinds(
        "Fermion",
        N;
        conserve_qns=true,
    )

    psi0 = MPS(
        sites,
        states,
    )

    gates = build_nnn_test_gates(
        sites;
        t2=t2,
        dt=dt,
        include_F=include_F,
    )

    psi_tebd = apply(
        gates,
        psi0;
        cutoff=1e-14,
        maxdim=1000,
    )

    psi_tebd_dense =
        mps_to_dense(
            psi_tebd,
            sites,
        )

    # -------------------------
    # Diagnostics
    # -------------------------

    norm_exact = norm(psi_exact)
    norm_tebd = norm(psi_tebd_dense)

    overlap =
        dot(
            psi_exact,
            psi_tebd_dense,
        )

    fidelity = abs(overlap)^2

    # Remove arbitrary global phase.
    phase =
        abs(overlap) > 0 ?
        overlap / abs(overlap) :
        1.0

    phase_aligned_error =
        norm(
            psi_tebd_dense -
            conj(phase) * psi_exact,
        )

    return (
        N=N,
        Nf=Nf,
        dt=dt,
        include_F=include_F,
        states=states,
        norm_exact=norm_exact,
        norm_tebd=norm_tebd,
        phase_aligned_error=phase_aligned_error,
        fidelity=fidelity,
        psi_exact=psi_exact,
        psi_tebd=psi_tebd_dense,
        H=H,
    )
end


# ============================================================
# DIRECT CUSTOM-STATE TEST
#
# |110> is useful because the middle site is occupied, so the
# Jordan-Wigner parity factor is actually probed.
# ============================================================

function compare_nnn_custom_state(
    states::Vector{String};
    dt=1e-4,
    t2=1.0,
    include_F=false,
)

    N = length(states)

    p = ExtendedTBParams(
        t0=0.0,
        t2=t2,
        V=0.0,
    )

    # Exact
    psi0_exact = exact_product_vector(states)

    H = exact_extended_H(N, p)

    psi_exact =
        exp(-1im * dt * H) *
        psi0_exact

    # TEBD
    sites = siteinds(
        "Fermion",
        N;
        conserve_qns=true,
    )

    psi0 = MPS(
        sites,
        states,
    )

    gates = build_nnn_test_gates(
        sites;
        t2=t2,
        dt=dt,
        include_F=include_F,
    )

    psi_tebd = apply(
        gates,
        psi0;
        cutoff=1e-14,
        maxdim=1000,
    )

    psi_tebd_dense =
        mps_to_dense(
            psi_tebd,
            sites,
        )

    overlap =
        dot(
            psi_exact,
            psi_tebd_dense,
        )

    fidelity = abs(overlap)^2

    phase =
        abs(overlap) > 0 ?
        overlap / abs(overlap) :
        1.0

    error =
        norm(
            psi_tebd_dense -
            conj(phase) * psi_exact,
        )

    return error, fidelity, norm(psi_tebd_dense)
end


# ============================================================
# MAIN TEST 1: N=3,4,5
# ============================================================

println()
println("="^78)
println("PURE NNN FERMIONIC TEBD vs EXACT DIAGONALIZATION")
println("="^78)
println()
println("Hamiltonian: t0=0, t2=1, V=0")
println("dt = 1e-4")
println("cutoff = 1e-14")
println("maxdim = 1000")
println("Fermion sites: conserve_qns=true")
println()

dt = 1e-4

for N in 3:5

    println("-"^78)
    println("N = $N")
    println("-"^78)

    r_noF = compare_nnn_tebd_ed(
        N;
        dt=dt,
        t2=1.0,
        Nf=N ÷ 2,
        include_F=false,
    )

    @printf(
        "WITHOUT F: states=%s  norm=%.12f  error=%.8e  fidelity=%.12f\n",
        join(r_noF.states, ""),
        r_noF.norm_tebd,
        r_noF.phase_aligned_error,
        r_noF.fidelity,
    )

    r_withF = compare_nnn_tebd_ed(
        N;
        dt=dt,
        t2=1.0,
        Nf=N ÷ 2,
        include_F=true,
    )

    @printf(
        "WITH F:    states=%s  norm=%.12f  error=%.8e  fidelity=%.12f\n",
        join(r_withF.states, ""),
        r_withF.norm_tebd,
        r_withF.phase_aligned_error,
        r_withF.fidelity,
    )

    println()
end


# ============================================================
# MAIN TEST 2: DIRECT |110> TEST
#
# N=3, |110> is deliberately chosen instead of the alternating
# product_state_Nf(3,2)=|101>, because it probes the middle-site
# parity factor directly.
# ============================================================

println()
println("="^78)
println("CRITICAL DIRECT |110> NNN TEST")
println("="^78)

states_110 = ["1", "1", "0"]

err_noF, fid_noF, norm_noF =
    compare_nnn_custom_state(
        states_110;
        dt=1e-4,
        t2=1.0,
        include_F=false,
    )

err_withF, fid_withF, norm_withF =
    compare_nnn_custom_state(
        states_110;
        dt=1e-4,
        t2=1.0,
        include_F=true,
    )

@printf(
    "WITHOUT F: norm=%.12f  error=%.8e  fidelity=%.12f\n",
    norm_noF,
    err_noF,
    fid_noF,
)

@printf(
    "WITH F:    norm=%.12f  error=%.8e  fidelity=%.12f\n",
    norm_withF,
    err_withF,
    fid_withF,
)


# ============================================================
# MAIN TEST 3: dt CONVERGENCE FOR N=5
#
# This checks that the residual error of the correct construction
# decreases as the Trotter step is reduced.
# ============================================================

println()
println("="^78)
println("N=5 dt-CONVERGENCE TEST (WITHOUT F)")
println("="^78)

for dt_test in [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-4]

    r = compare_nnn_tebd_ed(
        5;
        dt=dt_test,
        t2=1.0,
        Nf=2,
        include_F=false,
    )

    @printf(
        "dt=% .1e  error=% .8e  fidelity=%.12f  norm=%.12f\n",
        dt_test,
        r.phase_aligned_error,
        r.fidelity,
        r.norm_tebd,
    )
end

println()
println("="^78)
println("TEST COMPLETE")
println("="^78)
