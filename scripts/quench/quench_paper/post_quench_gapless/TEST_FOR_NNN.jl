push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis
using ITensors
using LinearAlgebra

# ---------------------------------------------------------------------
# Threading setup — do this BEFORE building any sites/tensors.
#
# You have 24 logical cores essentially idle (load avg ~1.2). This won't
# change *when* the entanglement wall hits, but it changes how much wall
# clock each χ costs you, which is what lets you push the wall further
# out in practice.
#
#   * Launch this script with:   julia -t 8 run_quench_paper_patched.jl
#     (or `julia -t auto`, or set JULIA_NUM_THREADS=8 in the environment)
#
#   * enable_threaded_blocksparse is specifically for QN-conserving
#     tensors (which you already have via conserve_qns=true) — it lets
#     ITensor parallelize contractions across the QN blocks.
#
#   * IMPORTANT: block-sparse threading parallelizes ACROSS your QN blocks
#     using Julia threads. It is coarse-grained parallelism, one block per
#     thread. If BLAS/Strided.jl then ALSO try to multithread the (small)
#     linear algebra work inside each individual block, you get nested
#     parallelism fighting over the same cores (8 Julia threads x 3 BLAS
#     threads = 24 threads contending, with only 24 physical cores and no
#     room left for anything else, plus constant thread-handoff overhead).
#     ITensors' own advice — confirmed by the warnings you saw — is to set
#     BOTH BLAS and Strided.jl to single-threaded and let the Julia threads
#     do all the parallel work. Do NOT split 24 cores proportionally here;
#     that was wrong in an earlier version of this script.
# ---------------------------------------------------------------------
using ITensors.NDTensors: Strided

n_cores = Threads.nthreads()
if n_cores == 1
    @warn "Julia started with only 1 thread. Relaunch with `julia -t 12`!"
end

# 1. Turn OFF block-sparse threading globally to banish the [13] crash
ITensors.disable_threaded_blocksparse()

# 2. Turn ON BLAS threading so dense matrix math uses all 12 cores!
BLAS.set_num_threads(n_cores)

# 3. Turn ON Strided threading to help with tensor permutations
Strided.set_num_threads(n_cores)

println("Julia threads: $n_cores, BLAS threads: $(BLAS.get_num_threads()), " *
        "Strided threads: $(Strided.get_num_threads())")

using LinearAlgebra
using ITensors
using ITensorMPS
using Printf

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
# TEBD GATES
#
# include_F = false:
#   proposed corrected implementation
#
# include_F = true:
#   your original implementation
# ============================================================

function build_test_tebd_gates(
    sites,
    p::ExtendedTBParams;
    dt=1e-4,
    include_F=false,
)

    N = length(sites)
    gates = ITensor[]

    # --------------------------------------------------------
    # Forward half step
    # --------------------------------------------------------

    for i in 1:(N - 2)

        s1 = sites[i]
        s2 = sites[i + 1]
        s3 = sites[i + 2]

        Id2 = op("Id", s2)
        Id3 = op("Id", s3)

        h = zero(op("Id", s1) * Id2 * Id3)

        # NN hopping
        if p.t0 != 0.0
            h = -p.t0 * op("Cdag", s1) * op("C", s2) * Id3
            h += -p.t0 * op("Cdag", s2) * op("C", s1) * Id3
        end

        # NN interaction
        if p.V != 0.0

            h += p.V *
                op("N", s1) *
                op("N", s2) *
                Id3

            h += -(p.V / 2) *
                op("N", s1) *
                Id2 *
                Id3

            h += -(p.V / 2) *
                op("Id", s1) *
                op("N", s2) *
                Id3
        end

        # ----------------------------------------------------
        # NNN hopping
        # ----------------------------------------------------

        if p.t2 != 0.0

            if include_F

                # ORIGINAL VERSION
                h += -p.t2 *
                    op("Cdag", s1) *
                    op("F", s2) *
                    op("C", s3)

                h += -p.t2 *
                    op("Cdag", s3) *
                    op("F", s2) *
                    op("C", s1)

            else

                # CORRECTED VERSION
                #
                # No explicit Jordan-Wigner F.
                h += -p.t2 *
                    op("Cdag", s1) *
                    Id2 *
                    op("C", s3)

                h += -p.t2 *
                    op("Cdag", s3) *
                    Id2 *
                    op("C", s1)
            end
        end

        push!(gates, exp(-1im * dt / 2 * h))
    end

    # --------------------------------------------------------
    # Final NN bond
    # --------------------------------------------------------

    s1 = sites[N - 1]
    s2 = sites[N]

    h_end = zero(op("Id", s1) * op("Id", s2))

    if p.t0 != 0.0

        h_end =
            -p.t0 *
            op("Cdag", s1) *
            op("C", s2)

        h_end +=
            -p.t0 *
            op("Cdag", s2) *
            op("C", s1)
    end

    if p.V != 0.0

        h_end +=
            p.V *
            op("N", s1) *
            op("N", s2)

        h_end +=
            -(p.V / 2) *
            op("N", s1) *
            op("Id", s2)

        h_end +=
            -(p.V / 2) *
            op("Id", s1) *
            op("N", s2)
    end

    push!(gates, exp(-1im * dt / 2 * h_end))

    # --------------------------------------------------------
    # Reverse sweep
    # --------------------------------------------------------

    append!(gates, reverse(gates))

    return gates
end


# ============================================================
# CONVERT MPS -> DENSE VECTOR
#
# Site ordering is explicitly restored to sites[1],...,sites[N].
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
# EXACT INITIAL BASIS VECTOR
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
            # nothing
        else
            error("Unknown state $(states[i])")
        end
    end

    v[state + 1] = 1.0

    return v
end


# ============================================================
# ONE TEBD STEP
# ============================================================

function tebd_one_step(
    N::Int,
    p::ExtendedTBParams;
    dt=1e-4,
    Nf=N ÷ 2,
    include_F=false,
)

    sites = siteinds("Fermion", N; conserve_qns=true)

    states = product_state_Nf(N, Nf)

    psi = MPS(sites, states)

    gates = build_test_tebd_gates(
        sites,
        p;
        dt=dt,
        include_F=include_F,
    )

    psi = apply(
        gates,
        psi;
        cutoff=1e-14,
        maxdim=1000,
    )

    normalize!(psi)

    return sites, states, psi
end


# ============================================================
# RUN ONE ED VS TEBD COMPARISON
# ============================================================

function compare_tebd_ed(
    N::Int;
    dt=1e-4,
    t0=0.0,
    t2=1.0,
    V=0.0,
    Nf=N ÷ 2,
    include_F=false,
)

    p = ExtendedTBParams(
        t0=t0,
        t2=t2,
        V=V,
    )

    # --------------------------------------------------------
    # Exact initial state
    # --------------------------------------------------------

    states = product_state_Nf(N, Nf)

    psi0_exact = exact_product_vector(states)

    # --------------------------------------------------------
    # Exact evolution
    # --------------------------------------------------------

    H = exact_extended_H(N, p)

    U = exp(-1im * dt * H)

    psi_exact = U * psi0_exact

    # --------------------------------------------------------
    # TEBD evolution
    # --------------------------------------------------------

    sites = siteinds("Fermion", N; conserve_qns=true)

    psi0 = MPS(sites, states)

    gates = build_test_tebd_gates(
        sites,
        p;
        dt=dt,
        include_F=include_F,
    )

    psi_tebd = apply(
        gates,
        psi0;
        cutoff=1e-14,
        maxdim=1000,
    )

    normalize!(psi_tebd)

    psi_tebd_dense = mps_to_dense(psi_tebd, sites)

    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------

    norm_exact = norm(psi_exact)
    norm_tebd  = norm(psi_tebd_dense)

    # Direct vector error
    direct_error = norm(psi_tebd_dense - psi_exact)

    # Phase-independent error
    #
    # This matters because two vectors that differ only by a
    # global phase are physically equivalent.
    overlap = dot(psi_exact, psi_tebd_dense)

    phase = abs(overlap) > 0 ? overlap / abs(overlap) : 1.0

    phase_aligned_error =
        norm(psi_tebd_dense - conj(phase) * psi_exact)

    fidelity = abs(overlap)^2

    return (
        N=N,
        Nf=Nf,
        dt=dt,
        include_F=include_F,
        norm_exact=norm_exact,
        norm_tebd=norm_tebd,
        direct_error=direct_error,
        phase_aligned_error=phase_aligned_error,
        fidelity=fidelity,
        states=states,
        psi_exact=psi_exact,
        psi_tebd=psi_tebd_dense,
        H=H,
    )
end


# ============================================================
# PRINT RESULT
# ============================================================

function print_result(r)

    @printf(
        "N=%d  Nf=%d  F=%s  direct=% .4e  phase-aligned=% .4e  fidelity=%.12f\n",
        r.N,
        r.Nf,
        r.include_F ? "YES" : "NO ",
        r.direct_error,
        r.phase_aligned_error,
        r.fidelity,
    )
end


# ============================================================
# MAIN TEST
# ============================================================

println()
println("="^78)
println("NNN FERMIONIC TEBD vs EXACT DIAGONALIZATION")
println("="^78)

println()
println("Parameters:")
println("  t0 = 0")
println("  t2 = 1")
println("  V  = 0")
println("  dt = 1e-4")
println("  cutoff = 1e-14")
println()

dt = 1e-4

for N in 3:5

    println("-"^78)
    println("N = $N")
    println("-"^78)

    # --------------------------------------------------------
    # WITHOUT F
    # --------------------------------------------------------

    r_noF = compare_tebd_ed(
        N;
        dt=dt,
        t0=0.0,
        t2=1.0,
        V=0.0,
        Nf=N ÷ 2,
        include_F=false,
    )

    print_result(r_noF)

    # --------------------------------------------------------
    # WITH F
    # --------------------------------------------------------

    r_withF = compare_tebd_ed(
        N;
        dt=dt,
        t0=0.0,
        t2=1.0,
        V=0.0,
        Nf=N ÷ 2,
        include_F=true,
    )

    print_result(r_withF)

    println()
end

println("="^78)
println("TEST COMPLETE")
println("="^78)

function compare_custom_state(
    states::Vector{String};
    dt=1e-4,
    t0=0.0,
    t2=1.0,
    V=0.0,
    include_F=false,
)

    N = length(states)

    p = ExtendedTBParams(
        t0=t0,
        t2=t2,
        V=V,
    )

    # Exact
    psi0_exact = exact_product_vector(states)

    H = exact_extended_H(N, p)

    psi_exact =
        exp(-1im * dt * H) * psi0_exact

    # TEBD
    sites = siteinds(
        "Fermion",
        N;
        conserve_qns=true,
    )

    psi0 = MPS(sites, states)

    gates = build_test_tebd_gates(
        sites,
        p;
        dt=dt,
        include_F=include_F,
    )

    psi_tebd = apply(
        gates,
        psi0;
        cutoff=1e-14,
        maxdim=1000,
    )

    normalize!(psi_tebd)

    psi_tebd_dense =
        mps_to_dense(psi_tebd, sites)

    overlap = dot(
        psi_exact,
        psi_tebd_dense,
    )

    phase =
        abs(overlap) > 0 ?
        overlap / abs(overlap) :
        1.0

    err =
        norm(
            psi_tebd_dense -
            conj(phase) * psi_exact,
        )

    fidelity = abs(overlap)^2

    return err, fidelity
end

println()
println("="^78)
println("DIRECT |110> TEST")
println("="^78)

states = ["1", "1", "0"]

err_noF, fid_noF =
    compare_custom_state(
        states;
        dt=1e-4,
        t0=0.0,
        t2=1.0,
        V=0.0,
        include_F=false,
    )

err_withF, fid_withF =
    compare_custom_state(
        states;
        dt=1e-4,
        t0=0.0,
        t2=1.0,
        V=0.0,
        include_F=true,
    )

@printf(
    "Without F: error = %.8e, fidelity = %.12f\n",
    err_noF,
    fid_noF,
)

@printf(
    "With F:    error = %.8e, fidelity = %.12f\n",
    err_withF,
    fid_withF,
)