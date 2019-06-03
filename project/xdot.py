from numpy import cos, sin

def xdot(t, x, V, kp, ki, R, B0, wg, V_base, XL_base, w_base, P_ref, Q_ref):
    """
    returns change in state dynamics.
    x representst system state, where x = (theta_g, i_d, i_q, z_d, z_q, v_c_d, v_c_q)
    All quantities are in per unit except theta_g, which is in rad/sec
    
    Arguments:
        t {float} -- time
        x {tuple} -- tuple containing state
        V {float} -- nominal voltage
        kp {float} -- current controller proportional gain
        ki {float} -- current controller integral gain
        R {float} -- per unit resistance of RL branch
        B0 {float} -- per unit susceptance of LCL cap
        V_base {float} -- voltage base for per-unitization
        XL_base {float} -- reactance base for inductive impedance per-unitization
        w_base {float} -- frequency base for per-unitization
        P_ref {float} -- input real power reference
        Q_ref {float} -- input reactive power reference
    
    Returns:
        tuple -- change in state variablles
    """
    # x = [thetag id iq zd zq vc]

    (theta_g, i_d, i_q, z_d, z_q, v_c_d, v_c_q) = x

    # grid angle
    dtheta_g = wg

    # grid voltage
    v_d = V * cos(theta_g - theta_g) / V_base
    v_q = V * sin(theta_g - theta_g) / V_base

    # reference currents
    i_d_ref = (P_ref * v_d + Q_ref * v_q) / (v_d**2 + v_q**2)
    i_q_ref = (P_ref * v_q - Q_ref * v_d) / (v_d**2 + v_q**2)

    # current control error
    e_id = i_d_ref - i_d
    e_iq = i_q_ref - i_q

    # terminal voltage
    vt_d = z_d + (e_id * kp) - (XL_base * (dtheta_g/w_base) * i_q) + v_d
    vt_q = z_q + (e_iq * kp) + (XL_base * (dtheta_g/w_base) * i_d) + v_q

    # change in currents and current controller states
    di_d = (w_base / XL_base) * (vt_d - (i_d * R) - v_d + (XL_base * (dtheta_g/w_base) * i_q))
    di_q = (w_base / XL_base) * (vt_q - (i_q * R) - v_q - (XL_base * (dtheta_g/w_base) * i_d))
    dz_d = w_base * ki * e_id
    dz_q = w_base * ki * e_iq

    # change in LCL cap voltage
    # TODO: solve for grid current i_g_d and i_g_q
    # i_g_dq = i_dq - (C0 * dv_c_d)
    dv_c_d = 0.0 #(w_base / B0) * (i_d - i_g_d + dtheta_g * B0 * v_c_q)
    dv_c_q = 0.0 #(w_base / B0) * (i_q - i_g_q - dtheta_g * B0 * v_c_d)

    return (dtheta_g, di_d, di_q, dz_d, dz_q, dv_c_d, dv_c_q)