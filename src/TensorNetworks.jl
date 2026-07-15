struct SSHParams
    v::Float64
    w::Float64
    Δ::Float64
    V::Float64
end
SSHParams(; v=1.0, w=1.0, Δ=0.0, V=0.0) = SSHParams(v, w, Δ, V)

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

function measure_observables(psi, N, V)
    La_width::Int = N ÷ 2
    subsystem_inds = 1:La_width

    # 1. Entanglement Entropy
    orthogonalize!(psi, La_width) 
    row_inds = uniqueinds(psi[La_width], psi[La_width+1])
    U, S_vals, V_mat = svd(psi[La_width], row_inds)

    S_EE = 0.0
    for n in 1:dim(S_vals, 1)
        λ = S_vals[n, n]
        p = λ^2
        if p > 1e-12 # Avoid log(0) issues
            S_EE -= p * log(p)
        end
    end
    S_EE = real(S_EE) 

    # 2. Particle Number Fluctuations
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

function create_CDW_states(N::Int)
    # 1. Ground state: perfectly alternating |1 0 1 0 1 0 ...>
    st0 = fill("0", N)
    for i in 1:2:N; st0[i] = "1"; end
    
    # 2. Ep State: Add particle to an even site as close to the center as possible
    stp = copy(st0)
    center_even = (N ÷ 2) % 2 == 0 ? (N ÷ 2) : (N ÷ 2) + 1
    stp[center_even] = "1"
    
    # 3. Em State: Remove particle from an odd site as close to the center as possible
    stm = copy(st0)
    center_odd = (N ÷ 2) % 2 != 0 ? (N ÷ 2) : (N ÷ 2) + 1
    stm[center_odd] = "0"
    
    return st0, stp, stm
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
                          resume_from=nothing)

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

        if step >= start_step && length(S_EE_vals) < step
            # only measure if we don't already have this step's values from a checkpoint
            S_EE, F = measure_observables(psi, N_sites, use_interacting_formula ? 1.0 : 0.0)
            push!(S_EE_vals, S_EE)
            push!(F_vals, F)
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

# function build_SSH_MPO_OBC_disordered(sites, p::SSHParams; disorder_strength=0.5)
#     N = length(sites)
#     os = OpSum()

#     # Hopping terms with disorder
#     for i in 1:2:N-1
#         v_disorder = p.v + disorder_strength * (rand() - 0.5) 
#         os += -v_disorder, "Cdag", i, "C", i+1; os += -v_disorder, "Cdag", i+1, "C", i
#     end
#     for i in 2:2:N-1
#         w_disorder = p.w + disorder_strength * (rand() - 0.5)
#         os += -w_disorder, "Cdag", i, "C", i+1; os += -w_disorder, "Cdag", i+1, "C", i
#     end

#     # Staggered Potential
#     if p.Δ != 0.0
#         for i in 1:N; os += p.Δ * (-1)^i, "N", i; end
#     end

#     # Nearest-neighbor interaction
#     if p.V != 0.0
#         for i in 1:N-1
#             os += p.V, "N", i, "N", i+1
#             os += -p.V/2.0, "N", i
#             os += -p.V/2.0, "N", i+1
#         end
#     end

#     return MPO(os, sites)

# end

# function simulate_quench_initial_disorder(N_sites, T_max, dt, pre::SSHParams, post::SSHParams; disorder_strength=0.5, bond_dim=1000)
#     # 1. Initialize sites with particle number conservation
#     sites = siteinds("Fermion", N_sites; conserve_qns=true)
    
#     # 2. Prepare the Initial State (t < 0)
#     println("Preparing interacting ground state (t < 0)...")
#     H_init = build_SSH_MPO_OBC_disordered(sites, pre; disorder_strength=disorder_strength)
    
#     # Using half-filling based on product state function
#     init_state = product_state_Nf(N_sites, N_sites ÷ 2) 
#     E0, psi = ground_energy(H_init, sites; init_state=init_state, nsweeps=12)
#     println("Initial Ground State Energy: ", E0)
    
#     # 3. Build TEBD gates for post-quench Hamiltonian (t > 0)
#     println("Constructing TEBD gates for quench...")
#     gates = build_tebd_gates(sites, post; dt=dt)
    
#     # 4. Time Evolution Loop
#     times = 0.0:dt:T_max
#     S_EE_vals = Float64[]
#     F_vals = Float64[]
#     use_interacting_formula = !iszero(post.V) || !iszero(pre.V)

    
#     println("Starting time evolution...")
#     for t in times
#         # Measure observables before applying the time step
#         S_EE, F = measure_observables(psi, N_sites, use_interacting_formula ? 1.0 : 0.0)
#         push!(S_EE_vals, S_EE)
#         push!(F_vals, F)
        
#         # Apply the Trotter gates to evolve the state by dt
#         # 'cutoff' and 'maxdim' are critical here to manage entanglement growth
#         psi = apply(gates, psi; cutoff=1e-10, maxdim=bond_dim)
#         normalize!(psi) # Normalize after each full Trotter step
#         println("Time: ", round(t, digits=3), " | S_EE: ", round(S_EE, digits=4), " | F: ", round(F, digits=4))
#     end
    
#     return times, S_EE_vals, F_vals
# end