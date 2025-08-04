# --- Components ---
function Components(k, φ, w, v, V)
    # k: wavevector
    # φ: phase
    # w: intercell hopping
    # v: intracell hopping
    # V: potential
    h_x = v*cos(φ) + w*(cos(φ)*cos(k) - sin(φ)*sin(k))
    h_y = -v*sin(φ) + w*(sin(φ)*cos(k) + cos(φ)*sin(k))
    h_z = V
    return h_x, h_y, h_z   
end

# --- Normalized Bloch vector before quench ---
function BlochPre(k, φ_i, w_i, v_i, V_i)
    # k: wavevector
    # φ: phase
    # w: intercell hopping
    # v: intracell hopping
    # V: potential
    h_x, h_y, h_z = Components(k, φ_i, w_i, v_i, V_i)
    E_plus = sqrt(h_x^2 + h_y^2 + h_z^2)
    return [h_x, h_y, h_z] ./ E_plus
end

# --- Normalized Bloch vector after quench ---
function BlochPost(k, φ_f, w_f, v_f, V_f)
    # k: wavevector
    # φ: phase
    # w: intercell hopping
    # v: intracell hopping
    # V: potential
    h_x, h_y, h_z = Components(k, φ_f, w_f, v_f, V_f)
    E_plus = sqrt(h_x^2 + h_y^2 + h_z^2)
    return [h_x, h_y, h_z] ./ E_plus
end

# --- Non-Normalized Bloch vector after quench ---
function BlochPostNonNorm(k, φ_f, w_f, v_f, V_f)
    # k: wavevector
    # φ: phase
    # w: intercell hopping
    # v: intracell hopping
    # V: potential
    h_x, h_y, h_z = Components(k, φ_f, w_f, v_f, V_f)
    return [h_x, h_y, h_z]
end

# --- Time-evolved Bloch vector ---
function BlochEvolved(k, pre_quench, post_quench, t)
    # pre_quench: parameters before the quench
    # post_quench: parameters after the quench
    # t: time
    φ_i, w_i, v_i, V_i = pre_quench
    φ_f, w_f, v_f, V_f = post_quench

    b_pre = BlochPre(k, φ_i, w_i, v_i, V_i)
    b_post = BlochPost(k, φ_f, w_f, v_f, V_f)
    bloch_post_non_norm = BlochPostNonNorm(k, φ_f, w_f, v_f, V_f)

    h1 = dot(b_pre, b_post)*b_post
    h2 = (b_pre - dot(b_pre, b_post) * b_post) * cos(2 *norm(bloch_post_non_norm) * t)
    h3 = cross(b_post, b_pre) * sin(2 * norm(bloch_post_non_norm) * t)

    return h1 + h2 + h3
end
