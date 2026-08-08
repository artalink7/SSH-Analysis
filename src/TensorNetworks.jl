using LinearAlgebra
using ITensors.NDTensors: Strided
struct SSHParams
    v::Float64
    w::Float64
    Δ::Float64
    V::Float64
end
SSHParams(; v=1.0, w=1.0, Δ=0.0, V=0.0) = SSHParams(v, w, Δ, V)
struct ExtendedTBParams
    t0::Float64
    t2::Float64
    V::Float64
end
ExtendedTBParams(; t0=1.0, t2 = 0.0, V=0.0) = ExtendedTBParams(t0, t2, V)

function build_extended_tb_MPO_OBC(sites; p::ExtendedTBParams)
    N = length(sites)
    os = OpSum()

    # Nearest-neighbor tight-binding hopping
    for i in 1:N-1
        os += -p.t0, "Cdag", i, "C", i+1
        os += -p.t0, "Cdag", i+1, "C", i
    end

    # Next-nearest-neighbor hopping
    if p.t2 != 0.0
        for i in 1:N-2 
            os += -p.t2, "Cdag", i, "C", i+2
            os += -p.t2, "Cdag", i+2, "C", i 
        end
    end

    # Nearest-neighbor density-density interaction: V (n_i -1/2)(n_{i+1} - 1/2)
    if p.V !=0 
        for i in 1:N-1
            os += p.V, "N", i, "N", i+1
            os += -p.V/2.0, "N", i 
            os += -p.V/2.0, "N", i+1
        end
    end

    return MPO(os, sites)
end

function build_SSH_MPO_OBC(sites; p::SSHParams)
    N = length(sites)
    os = OpSum()

    # Hopping terms
    for i in 1:2:N-1
        os += -p.v, "Cdag", i, "C", i+1; os += -p.v, "Cdag", i+1, "C", i
    end
    for i in 2:2:N-1
        os += -p.w, "Cdag", i, "C", i+1; os += -p.w, "Cdag", i+1, "C", i
    end

    # Staggered Potential
    if p.Δ != 0.0
        for i in 1:N; os += p.Δ * (-1)^i, "N", i; end
    end

    # Nearest-neighbor interaction
    if p.V != 0.0
        for i in 1:N-1
            os += p.V, "N", i, "N", i+1
            os += -p.V/2.0, "N", i
            os += -p.V/2.0, "N", i+1
        end
    end

    return MPO(os, sites)
end

# ---------------------------------------------------------------------
# Bulk-interval entanglement entropy for an interval [a,b] that may or
# may not touch a physical edge (a=1 or b=N).
#
#   - If A touches an edge, it's a genuine single-cut bipartition and
#     we just do the cheap SVD at the correct bond (same math as your
#     original code, just generalized to whichever end is the "cut").
#   - If A touches NEITHER edge, rho_A is mixed and needs the two-sided
#     construction: contract ket-bra pairs site by site over sites a..b,
#     summing the physical index at each step and carrying only the two
#     boundary link indices (unprimed=ket, primed=bra) forward. This
#     never forms a 2^(b-a+1)-dimensional object -- cost stays
#     O(width * chi^3 * d), same order as a normal MPS overlap.
#     Its spectrum equals the entanglement spectrum of A exactly.
# ---------------------------------------------------------------------

function single_cut_entropy(psi::MPS, cut::Int)
    orthogonalize!(psi, cut)
    row_inds = uniqueinds(psi[cut], psi[cut + 1])
    _, S_vals, _ = svd(psi[cut], row_inds)
    S = 0.0
    for n in 1:dim(S_vals, 1)
        p = S_vals[n, n]^2
        if p > 1e-12
            S -= p * log(p)
        end
    end
    return real(S)
end

function bulk_entanglement_entropy(psi::MPS, a::Int, b::Int)
    orthogonalize!(psi, a)

    # Running transfer tensor: after processing site j, E carries only
    # the (still-open) left boundary link (a-1|a) and whatever the
    # current right boundary link is, each doubled into ket/bra copies.
    ket = psi[a]
    bra = prime(dag(psi[a]), "Link")
    E = ket * bra   # contracts the shared physical index at site a

    for j in (a + 1):b
        ket = psi[j]
        bra = prime(dag(psi[j]), "Link")
        E = E * ket   # extend ket-side link
        E = E * bra   # extend bra-side link, contract physical index s_j
    end

    # E now has exactly 4 open indices: the ket/bra copies of the left
    # boundary link (a-1|a) and the ket/bra copies of the right boundary
    # link (b|b+1). Diagonalize it as a (link x link) matrix.
    row_inds = filter(i -> plev(i) == 0, inds(E))
    col_inds = filter(i -> plev(i) == 1, inds(E))
    D, _ = eigen(E, row_inds, col_inds)

    S = 0.0
    for n in 1:dim(D, 1)
        λ = real(D[n, n])
        if λ > 1e-12
            S -= λ * log(λ)
        end
    end
    return S
end

function measure_observables(psi::MPS, N::Int, V::Float64; a::Int, b::Int)
    subsystem_inds = a:b

    S_EE = if a == 1
        single_cut_entropy(psi, b)          # A=[1,b]: cut at bond b|b+1
    elseif b == N
        single_cut_entropy(psi, a - 1)       # A=[a,N]: cut at bond a-1|a
    else
        bulk_entanglement_entropy(psi, a, b) # genuine bulk interval, two cuts
    end

    if iszero(V)
        CM = correlation_matrix(psi, "Cdag", "C", sites = subsystem_inds)
        F = real(tr(CM * (I - CM)))
    else
        NM = correlation_matrix(psi, "N", "N", sites = subsystem_inds)
        N_A = real(tr(NM))
        F = real(sum(NM) - N_A^2)
    end

    return S_EE, F
end

function measure_F(psi::MPS, N::Int, V::Float64; a::Int, b::Int)
    subsystem_inds = a:b
    if iszero(V)
        CM = correlation_matrix(psi, "Cdag", "C", sites = subsystem_inds)
        F = real(tr(CM * (I - CM)))
    else
        NM = correlation_matrix(psi, "N", "N", sites = subsystem_inds)
        N_A = real(tr(NM))
        F = real(sum(NM) - N_A^2)
    end
    return F
end

function measure_S_EE(psi::MPS, N::Int; a::Int, b::Int)
    if a==1
        S_EE = single_cut_entropy(psi, b)
    elseif b==N
        S_EE = single_cut_entropy(psi, a-1)
    else
        S_EE = bulk_entanglement_entropy(psi, a, b)
    end
    return S_EE
end

function product_state_Nf(N::Int, Nf::Int; pattern::Symbol=:alternating)
    @assert 0<= Nf <= N
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

function ground_energy(H::MPO, sites; init_state, nsweeps=12)
    # 1. Use randomMPS to break the strict 1,0,1,0 pinning
    psi0 = randomMPS(sites, init_state)
    
    # 2. Use the Sweeps object to apply enough noise to allow 
    # the edge states to tunnel across the chain and symmetrize
    sweeps = Sweeps(nsweeps)
    maxdim!(sweeps, 50, 100, 200, 400, 800)
    cutoff!(sweeps, 1e-10)
    noise!(sweeps, 1E-4, 1E-5, 1E-6, 1E-7, 0.0) # Slightly higher initial noise
    
    energy, psi = dmrg(H, psi0, sweeps; outputlevel=1)
    return energy, psi
end

function build_tebd_gates(sites, p::SSHParams; dt=0.05)
    N = length(sites)
    gates = ITensor[]

    # Forward sweep: dt/2
    for j in 1:(N - 1)
        s1 = sites[j]
        s2 = sites[j + 1]

        # Alternating hopping parameters for SSH
        hop = (j % 2 != 0) ? p.v : p.w

        # Non-interacting bond Hamiltonian
        hj = -hop * op("Cdag", s1) * op("C", s2) - hop * op("Cdag", s2) * op("C", s1)

        # Nearest-neighbor interaction
        if p.V != 0.0
            Id1 = op("Id", s1)
            Id2 = op("Id", s2)
            hj += p.V * op("N", s1) * op("N", s2)
            hj += -(p.V/2) * op("N", s1) * Id2
            hj += -(p.V/2) * Id1 * op("N", s2)
        end

        # Optional staggered potential 
        if p.Δ != 0.0
            Id1 = op("Id", s1)
            Id2 = op("Id", s2)
            hj += p.Δ * (-1)^j * op("N", s1) * Id2
            if j == N - 1
                # pick up the last site's on-site term on the final bond
                hj += p.Δ * (-1)^(j+1) * Id1 * op("N", s2)
            end
        end

        Gj = exp(-im * dt / 2 * hj)
        push!(gates, Gj)
    end

    # Reverse sweep: dt/2 to complete the 2nd order Trotter step
    append!(gates, reverse(gates))

    println("Constructed ", length(gates), " TEBD gates for quench evolution.")

    return gates
end

function build_extended_tb_tebd_gates(sites, p::ExtendedTBParams; dt=0.05)
    N = length(sites)
    gates = ITensor[]

    # Forward sweep: dt/2
    for i in 1:(N - 2)
        s1, s2, s3 = sites[i], sites[i+1], sites[i+2]
        Id2, Id3 = op("Id", s2), op("Id", s3)

        # NN hopping + NN interaction on bond (i, i+1)
        h  = -p.t0 * op("Cdag", s1) * op("C", s2) * Id3
        h += -p.t0 * op("Cdag", s2) * op("C", s1) * Id3

        if p.V != 0.0
            h += p.V * op("N", s1) * op("N", s2) * Id3
            h += -(p.V/2) * op("N", s1) * Id2 * Id3
            h += -(p.V/2) * op("Id", s1) * op("N", s2) * Id3
        end

        # NNN hopping on bond (i, i+2) — Jordan-Wigner string included
        if p.t2 != 0.0
            Fs2 = op("F", s2)
            h += -p.t2 * op("Cdag", s1) * Fs2 * op("C", s3)
            h += -p.t2 * op("Cdag", s3) * Fs2 * op("C", s1)
        end

        push!(gates, exp(-im * dt / 2 * h))
    end

    # --- FIX 1: Add the missing nearest-neighbor terms for the final (N-1, N) bond ---
    s_penult, s_last = sites[N-1], sites[N]
    h_end = -p.t0 * op("Cdag", s_penult) * op("C", s_last)
    h_end += -p.t0 * op("Cdag", s_last) * op("C", s_penult)
    
    if p.V != 0.0
        h_end += p.V * op("N", s_penult) * op("N", s_last)
        h_end += -(p.V/2) * op("N", s_penult) * op("Id", s_last)
        h_end += -(p.V/2) * op("Id", s_penult) * op("N", s_last)
    end
    push!(gates, exp(-im * dt / 2 * h_end))

    # --- FIX 2: Append the reverse sweep to complete the 2nd-order Trotter step ---
    append!(gates, reverse(gates))

    return gates
end

#
# PATCHED VERSION of the quench-related functions from your SSHAnalysis module.
#
# What changed vs. your original, and why:
#
#   1. simulate_quench() now logs `norm_before` (the norm of psi right after
#      `apply(...)` but before `normalize!`). You were already computing this
#      and discarding it — it's actually your truncation-error diagnostic.
#      Since the TEBD gates are unitary, exact evolution preserves the norm;
#      any drop below 1 comes purely from SVD truncation at maxdim/cutoff.
#      We turn it into a running "cumulative fidelity" estimate:
#           total_fidelity *= norm_before^2
#      so you get a concrete, per-step and cumulative number for how much
#      the χ=1700 cap is actually costing you past t≈3.7, instead of just
#      knowing "χ hit the ceiling".
#
#   2. Checkpointing: psi + accumulated observables are written to an HDF5
#      file every `checkpoint_every` steps. This means (a) a crash doesn't
#      cost you the full run, and (b) you can extend/rerun the *tail* of a
#      quench at a larger bond_dim without redoing DMRG + the cheap early
#      dynamics.
#
#   3. Resume support: pass `resume_from=<checkpoint path>` to continue an
#      evolution from a saved state, optionally with a different bond_dim
#      or cutoff than the original run used.
#
#   4. Threading hooks: see quench_run_script_patched.jl — the physics code
#      here is unchanged except for the two additions above, since threading
#      is a global/process-level setting, not something that belongs inside
#      this function.
#
# Everything else (build_SSH_MPO_OBC, build_tebd_gates, measure_observables,
# product_state_Nf, ground_energy, create_CDW_states) is untouched — paste
# this in place of your existing simulate_quench, and add the two small
# helper functions (save_checkpoint / load_checkpoint) alongside it.

# ---------------------------------------------------------------------
# Checkpoint helpers
# ---------------------------------------------------------------------

function save_checkpoint(path, psi, step, t, S_EE_vals, F_vals, total_fidelity)
    h5open(path, "w") do file
        write(file, "psi", psi)
        write(file, "step", step)
        write(file, "t", t)
        write(file, "S_EE_vals", S_EE_vals)
        write(file, "F_vals", F_vals)
        write(file, "total_fidelity", total_fidelity)
    end
end

function load_checkpoint(path)
    return h5open(path, "r") do file
        psi = read(file, "psi", MPS)
        step = read(file, "step")
        t = read(file, "t")
        S_EE_vals = read(file, "S_EE_vals")
        F_vals = read(file, "F_vals")
        total_fidelity = read(file, "total_fidelity")
        (psi=psi, step=step, t=t, S_EE_vals=S_EE_vals, F_vals=F_vals,
         total_fidelity=total_fidelity)
    end
end

# ---------------------------------------------------------------------
# simulate_quench, patched
# ---------------------------------------------------------------------

function simulate_quench(N_sites, T_max, dt, pre::SSHParams, post::SSHParams;
                          bond_dim=1000, cutoff=1e-8, logfile=nothing,
                          checkpoint_file=nothing, checkpoint_every=20,
                          resume_from=nothing, 
                          subsystem_a=1, subsystem_b=N_sites ÷ 2, 
                          entropy_every=1)

    log_io = logfile === nothing ? nothing : open(logfile, "w")
    function log(msg)
        println(msg)
        if log_io !== nothing
            println(log_io, msg)
            flush(log_io)
        end
    end

    log("="^60)
    log("Quench simulation started: $(now())")
    log("N_sites = $N_sites, T_max = $T_max, dt = $dt, bond_dim = $bond_dim, cutoff = $cutoff")
    log("Pre-quench parameters: v=$(pre.v), w=$(pre.w), Δ=$(pre.Δ), V=$(pre.V)")
    log("Post-quench parameters: v=$(post.v), w=$(post.w), Δ=$(post.Δ), V=$(post.V)")
    log("Subsystem for S_EE / F: sites $subsystem_a:$subsystem_b " *
        "(touches left edge: $(subsystem_a==1), touches right edge: $(subsystem_b==N_sites))")
    if resume_from !== nothing
        log("Resuming from checkpoint: $resume_from")
    end
    log("="^60)

    times = 0.0:dt:T_max
    use_interacting_formula = !iszero(post.V) || !iszero(pre.V)

    local sites, psi, gates
    S_EE_vals = Float64[]
    F_vals = Float64[]
    start_step = 1
    total_fidelity = 1.0

    if resume_from === nothing
        # --- Fresh start: sites, ground state, gates as before ---
        sites = siteinds("Fermion", N_sites; conserve_qns=true)

        log("Preparing interacting ground state (t < 0)...")
        H_init = build_SSH_MPO_OBC(sites; p=pre)
        init_state = product_state_Nf(N_sites, N_sites ÷ 2)
        t_dmrg_start = time()
        E0, psi = ground_energy(H_init, sites; init_state=init_state, nsweeps=12)
        log("Initial Ground State Energy: $E0 (DMRG took $(round(time() - t_dmrg_start, digits=2)) seconds)")
        log("Initial max bond dimension: $(maxlinkdim(psi))")

        log("Constructing TEBD gates for quench...")
        gates = build_tebd_gates(sites, post; dt=dt)
    else
        # --- Resume: load psi + history, rebuild gates against its sites ---
        ckpt = load_checkpoint(resume_from)
        psi = ckpt.psi
        sites = siteinds(psi)
        S_EE_vals = ckpt.S_EE_vals
        F_vals = ckpt.F_vals
        total_fidelity = ckpt.total_fidelity
        start_step = ckpt.step + 1
        log("Loaded checkpoint at step $(ckpt.step), t = $(ckpt.t), " *
            "cumulative fidelity so far = $(round(total_fidelity, digits=6))")

        log("Constructing TEBD gates for quench (post-resume bond_dim=$bond_dim, cutoff=$cutoff)...")
        gates = build_tebd_gates(sites, post; dt=dt)
    end

    log("Starting time evolution...")
    t_evolution_start = time()

    for (step, t) in enumerate(times)
        step < start_step && continue  # skip steps already done before the checkpoint

        if step >= start_step && length(F_vals) < step
            # only measure if we don't already have this step's values from a checkpoint
            F = measure_F(psi, N_sites, use_interacting_formula ? 1.0 : 0.0; 
                                            a=subsystem_a, b=subsystem_b)
            push!(F_vals, F)
            if (step -1) % entropy_every == 0
                S_EE = measure_S_EE(psi, N_sites; a =subsystem_a, b=subsystem_b)
            else
                S_EE = NaN
            end
            push!(S_EE_vals, S_EE)
        else
            S_EE, F = S_EE_vals[step], F_vals[step]
        end

        step_start = time()

        psi = apply(gates, psi; cutoff=cutoff, maxdim=bond_dim)
        norm_before = norm(psi)          # <-- truncation-error diagnostic
        normalize!(psi)
        step_time = round(time() - step_start, digits=4)

        # Cumulative fidelity estimate: since the gates are unitary, any norm
        # loss at this step is purely truncation. norm_before^2 approximates
        # the fraction of weight kept at this step; multiplying across steps
        # gives a running lower-bound-ish estimate of overlap with the exact
        # (untruncated) state.
        total_fidelity *= norm_before^2

        chi = maxlinkdim(psi)
        elapsed_time = round(time() - t_evolution_start, digits=4)
        eta = (elapsed_time / (step - start_step + 1)) * (length(times) - step)

        log("t = $(round(t, digits=3)), S_EE = $(round(S_EE, digits=4)), " *
            "F = $(round(F, digits=4)), χ = $chi / $bond_dim, " *
            "norm_before = $(round(norm_before, digits=8)), " *
            "cum_fidelity = $(round(total_fidelity, digits=6)), " *
            "step_time = $step_time s, elapsed_time = $elapsed_time s, " *
            "estimated_remaining_time = $(round(eta, digits=4)) s")

        if chi >= bond_dim
            log("Warning: Maximum bond dimension reached at t = $t. Consider increasing bond_dim.")
        end
        if total_fidelity < 0.99
            log("Warning: cumulative truncation error has exceeded 1% (cum_fidelity = " *
                "$(round(total_fidelity, digits=6))) at t = $t. Dynamics beyond this point " *
                "should be treated with caution.")
        end

        if checkpoint_file !== nothing && step % checkpoint_every == 0
            save_checkpoint(checkpoint_file, psi, step, t, S_EE_vals, F_vals, total_fidelity)
            log("Checkpoint saved at step $step (t = $t) -> $checkpoint_file")
        end
    end

    if checkpoint_file !== nothing
        save_checkpoint(checkpoint_file, psi, length(times), times[end], S_EE_vals, F_vals, total_fidelity)
        log("Final checkpoint saved -> $checkpoint_file")
    end

    log("="^60)
    log("Quench simulation completed: $(now())")
    log("Total time evolution duration: $(round((time() - t_evolution_start)/60, digits=2)) min")
    log("Final cumulative fidelity estimate: $(round(total_fidelity, digits=6))")
    log("="^60)

    if log_io !== nothing
        close(log_io)
    end

    return times, S_EE_vals, F_vals
end

function simulate_quench_only_fluctuations(N_sites, T_max, dt, pre::SSHParams, post::SSHParams;
                          bond_dim=1000, cutoff=1e-8, logfile=nothing,
                          checkpoint_file=nothing, checkpoint_every=20,
                          resume_from=nothing, 
                          subsystem_a=1, subsystem_b=N_sites ÷ 2)

    log_io = logfile === nothing ? nothing : open(logfile, "w")
    function log(msg)
        println(msg)
        if log_io !== nothing
            println(log_io, msg)
            flush(log_io)
        end
    end

    log("="^60)
    log("Quench simulation started: $(now())")
    log("N_sites = $N_sites, T_max = $T_max, dt = $dt, bond_dim = $bond_dim, cutoff = $cutoff")
    log("Pre-quench parameters: v=$(pre.v), w=$(pre.w), Δ=$(pre.Δ), V=$(pre.V)")
    log("Post-quench parameters: v=$(post.v), w=$(post.w), Δ=$(post.Δ), V=$(post.V)")
    log("Subsystem for F: sites $subsystem_a:$subsystem_b " *
        "(touches left edge: $(subsystem_a==1), touches right edge: $(subsystem_b==N_sites))")
    if resume_from !== nothing
        log("Resuming from checkpoint: $resume_from")
    end
    log("="^60)

    times = 0.0:dt:T_max
    use_interacting_formula = !iszero(post.V) || !iszero(pre.V)

    local sites, psi, gates
    S_EE_vals = Float64[]
    F_vals = Float64[]
    start_step = 1
    total_fidelity = 1.0

    if resume_from === nothing
        # --- Fresh start: sites, ground state, gates as before ---
        sites = siteinds("Fermion", N_sites; conserve_qns=true)

        log("Preparing interacting ground state (t < 0)...")
        H_init = build_SSH_MPO_OBC(sites; p=pre)
        init_state = product_state_Nf(N_sites, N_sites ÷ 2)
        t_dmrg_start = time()
        E0, psi = ground_energy(H_init, sites; init_state=init_state, nsweeps=12)
        log("Initial Ground State Energy: $E0 (DMRG took $(round(time() - t_dmrg_start, digits=2)) seconds)")
        log("Initial max bond dimension: $(maxlinkdim(psi))")

        log("Constructing TEBD gates for quench...")
        gates = build_tebd_gates(sites, post; dt=dt)
    else
        # --- Resume: load psi + history, rebuild gates against its sites ---
        ckpt = load_checkpoint(resume_from)
        psi = ckpt.psi
        sites = siteinds(psi)
        S_EE_vals = ckpt.S_EE_vals
        F_vals = ckpt.F_vals
        total_fidelity = ckpt.total_fidelity
        start_step = ckpt.step + 1
        log("Loaded checkpoint at step $(ckpt.step), t = $(ckpt.t), " *
            "cumulative fidelity so far = $(round(total_fidelity, digits=6))")

        log("Constructing TEBD gates for quench (post-resume bond_dim=$bond_dim, cutoff=$cutoff)...")
        gates = build_tebd_gates(sites, post; dt=dt)
    end

    log("Starting time evolution...")
    t_evolution_start = time()

    for (step, t) in enumerate(times)
        step < start_step && continue  # skip steps already done before the checkpoint

        if step >= start_step && length(F_vals) < step
            # only measure if we don't already have this step's values from a checkpoint
            F = measure_F(psi, N_sites, use_interacting_formula ? 1.0 : 0.0; 
                                            a=subsystem_a, b=subsystem_b)
            push!(F_vals, F)
        else
            F = F_vals[step]
        end

        step_start = time()

        psi = apply(gates, psi; cutoff=cutoff, maxdim=bond_dim)
        norm_before = norm(psi)          # <-- truncation-error diagnostic
        normalize!(psi)
        step_time = round(time() - step_start, digits=4)

        # Cumulative fidelity estimate: since the gates are unitary, any norm
        # loss at this step is purely truncation. norm_before^2 approximates
        # the fraction of weight kept at this step; multiplying across steps
        # gives a running lower-bound-ish estimate of overlap with the exact
        # (untruncated) state.
        total_fidelity *= norm_before^2

        chi = maxlinkdim(psi)
        elapsed_time = round(time() - t_evolution_start, digits=4)
        eta = (elapsed_time / (step - start_step + 1)) * (length(times) - step)

        log("t = $(round(t, digits=3)), " *
            "F = $(round(F, digits=4)), χ = $chi / $bond_dim, " *
            "norm_before = $(round(norm_before, digits=8)), " *
            "cum_fidelity = $(round(total_fidelity, digits=6)), " *
            "step_time = $step_time s, elapsed_time = $elapsed_time s, " *
            "estimated_remaining_time = $(round(eta, digits=4)) s")

        if chi >= bond_dim
            log("Warning: Maximum bond dimension reached at t = $t. Consider increasing bond_dim.")
        end
        if total_fidelity < 0.99
            log("Warning: cumulative truncation error has exceeded 1% (cum_fidelity = " *
                "$(round(total_fidelity, digits=6))) at t = $t. Dynamics beyond this point " *
                "should be treated with caution.")
        end

        if checkpoint_file !== nothing && step % checkpoint_every == 0
            save_checkpoint(checkpoint_file, psi, step, t, S_EE_vals, F_vals, total_fidelity)
            log("Checkpoint saved at step $step (t = $t) -> $checkpoint_file")
        end
    end

    if checkpoint_file !== nothing
        save_checkpoint(checkpoint_file, psi, length(times), times[end], S_EE_vals, F_vals, total_fidelity)
        log("Final checkpoint saved -> $checkpoint_file")
    end

    log("="^60)
    log("Quench simulation completed: $(now())")
    log("Total time evolution duration: $(round((time() - t_evolution_start)/60, digits=2)) min")
    log("Final cumulative fidelity estimate: $(round(total_fidelity, digits=6))")
    log("="^60)

    if log_io !== nothing
        close(log_io)
    end

    return times, F_vals
end

function simulate_quench_only_fluctuations_extended_tebd(N_sites, T_max, dt, pre::ExtendedTBParams, post::ExtendedTBParams;
                          bond_dim=1000, cutoff=1e-8, logfile=nothing,
                          checkpoint_file=nothing, checkpoint_every=20,
                          resume_from=nothing,
                          subsystem_a=1, subsystem_b=N_sites ÷ 2)

    log_io = logfile === nothing ? nothing : open(logfile, "w")
    function log(msg)
        println(msg)
        if log_io !== nothing
            println(log_io, msg)
            flush(log_io)
        end
    end

    log("="^60)
    log("Quench simulation started: $(now())")
    log("N_sites = $N_sites, T_max = $T_max, dt = $dt, bond_dim = $bond_dim, cutoff = $cutoff")
    log("Pre-quench parameters: t0=$(pre.t0), t2=$(pre.t2), V=$(pre.V)")
    log("Post-quench parameters: t0=$(post.t0), w=$(post.t2), V=$(post.V)")
    log("Subsystem for F: sites $subsystem_a:$subsystem_b " *
        "(touches left edge: $(subsystem_a==1), touches right edge: $(subsystem_b==N_sites))")
    if resume_from !== nothing
        log("Resuming from checkpoint: $resume_from")
    end
    log("="^60)

    times = 0.0:dt:T_max
    use_interacting_formula = !iszero(post.V) || !iszero(pre.V)

    local sites, psi, gates
    F_vals = Float64[]
    start_step = 1
    total_fidelity = 1.0

    if resume_from === nothing
        sites = siteinds("Fermion", N_sites; conserve_qns=true)

        # Threading is managed in the run script now!
        
        H_init = build_extended_tb_MPO_OBC(sites; p=pre)
        init_state = product_state_Nf(N_sites, N_sites ÷ 2)
        E0, psi = ground_energy(H_init, sites; init_state=init_state, nsweeps=12)
        log("Initial Ground State Energy: $E0")
    else
        ckpt = load_checkpoint(resume_from)
        psi = ckpt.psi
        sites = siteinds(psi)
        F_vals = ckpt.F_vals
        total_fidelity = ckpt.total_fidelity
        start_step = ckpt.step + 1
        log("Loaded checkpoint at step $(ckpt.step), t = $(ckpt.t)")
    end

    # Build the gates using your corrected 3-site gate function
    gates = build_extended_tb_tebd_gates(sites, post; dt=dt)

    for (step, t) in enumerate(times)
        step < start_step && continue
        if length(F_vals) < step
            F = measure_F(psi, N_sites, use_interacting_formula ? 1.0 : 0.0; a=subsystem_a, b=subsystem_b)
            push!(F_vals, F)
        end

        psi = apply(gates, psi; cutoff=cutoff, maxdim=bond_dim)
        norm_before = norm(psi)
        normalize!(psi)
        total_fidelity *= norm_before^2
        chi = maxlinkdim(psi)

        log("t = $(round(t, digits=3)), F = $(round(F_vals[step], digits=4)), χ = $chi / $bond_dim, cum_fidelity = $(round(total_fidelity, digits=6))")

        if checkpoint_file !== nothing && step % checkpoint_every == 0
            save_checkpoint(checkpoint_file, psi, step, t, Float64[], F_vals, total_fidelity)
            log("Checkpoint saved at step $step")
        end
    end

    if log_io !== nothing; close(log_io); end
    return times, F_vals
end

function extract_bipartite_fluctuations_centered(N_sites::Int, p::ExtendedTBParams; 
                                                 l_min::Int=5, l_max::Union{Int, Nothing}=nothing,
                                                 logfile::Union{String, Nothing}=nothing)
    
    # Setup logger
    log_io = logfile === nothing ? nothing : open(logfile, "w")
    function log_msg(msg)
        println(msg)
        if log_io !== nothing
            println(log_io, msg)
            flush(log_io)
        end
    end

    if l_max === nothing
        l_max = N_sites - 10 
    end

    log_msg("="^60)
    log_msg("Fluctuation measurement started: $(now())")
    log_msg("N_sites = $N_sites, Subsystem range: l = $l_min to $l_max")
    log_msg("Parameters: t0=$(p.t0), t2=$(p.t2), V=$(p.V)")
    log_msg("="^60)
    
    # 1. Setup the system
    sites = siteinds("Fermion", N_sites; conserve_qns=true)
    H = build_extended_tb_MPO_OBC(sites; p=p)
    init_state = product_state_Nf(N_sites, N_sites ÷ 2)
    
    # 2. Get the ground state
    log_msg("Running DMRG to obtain the ground state (nsweeps=12)...")
    t_dmrg_start = time()
    E0, psi = ground_energy(H, sites; init_state=init_state, nsweeps=12)
    
    log_msg("DMRG completed in $(round(time() - t_dmrg_start, digits=2)) seconds.")
    log_msg("Ground state energy: $E0")
    log_msg("Final ground state max bond dimension: $(maxlinkdim(psi))")
    log_msg("-"^60)
    
    # 3. Measure fluctuations
    log_msg("Measuring bipartite charge fluctuations for centered subsystems...")
    l_vals = Int[]
    F_vals = Float64[]
    
    t_measure_start = time()
    for l in l_min:l_max
        a = div(N_sites - l, 2) + 1
        b = a + l - 1
        
        # Pass p.V so measure_F knows whether to use the interacting or non-interacting observable
        F = measure_F(psi, N_sites, p.V; a=a, b=b)
        push!(l_vals, l)
        push!(F_vals, F)
        
        # Log progress per length
        log_msg("  -> subsystem length l = $l (sites $a to $b): F = $(round(F, digits=6))")
    end
    
    log_msg("-"^60)
    log_msg("Measurement completed in $(round(time() - t_measure_start, digits=2)) seconds.")
    log_msg("="^60)
    
    if log_io !== nothing
        close(log_io)
    end
    
    return l_vals, F_vals
end

function extract_bipartite_fluctuations_edge(N_sites::Int, p::ExtendedTBParams; 
                                             l_min::Int=5, l_max::Union{Int, Nothing}=nothing,
                                             logfile::Union{String, Nothing}=nothing)
    
    # Setup logger
    log_io = logfile === nothing ? nothing : open(logfile, "w")
    function log_msg(msg)
        println(msg)
        if log_io !== nothing
            println(log_io, msg)
            flush(log_io)
        end
    end

    # For an edge subsystem, we measure up to the middle of the chain
    if l_max === nothing
        l_max = N_sites ÷ 2 
    end

    log_msg("="^60)
    log_msg("Fluctuation measurement started (EDGE SUBSYSTEM): $(now())")
    log_msg("N_sites = $N_sites, Subsystem range: l = $l_min to $l_max")
    log_msg("Parameters: t0=$(p.t0), t2=$(p.t2), V=$(p.V)")
    log_msg("="^60)
    
    sites = siteinds("Fermion", N_sites; conserve_qns=true)
    H = build_extended_tb_MPO_OBC(sites; p=p)
    init_state = product_state_Nf(N_sites, N_sites ÷ 2)
    
    log_msg("Running DMRG to obtain the ground state (nsweeps=12)...")
    t_dmrg_start = time()
    E0, psi = ground_energy(H, sites; init_state=init_state, nsweeps=12)
    
    log_msg("DMRG completed in $(round(time() - t_dmrg_start, digits=2)) seconds.")
    log_msg("Ground state energy: $E0")
    log_msg("Final ground state max bond dimension: $(maxlinkdim(psi))")
    log_msg("-"^60)
    
    log_msg("Measuring bipartite charge fluctuations for edge subsystems (1 to l)...")
    l_vals = Int[]
    F_vals = Float64[]
    
    t_measure_start = time()
    for l in l_min:l_max
        a = 1
        b = l
        
        F = measure_F(psi, N_sites, p.V; a=a, b=b)
        push!(l_vals, l)
        push!(F_vals, F)
        
        log_msg("  -> subsystem length l = $l (sites $a to $b): F = $(round(F, digits=6))")
    end
    
    log_msg("-"^60)
    log_msg("Measurement completed in $(round(time() - t_measure_start, digits=2)) seconds.")
    log_msg("="^60)
    
    if log_io !== nothing
        close(log_io)
    end
    
    return l_vals, F_vals
end