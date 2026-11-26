function [u_raw, state] = ctrl_lqri_bc(x_est, u_prev, ref, state, missile, cfg)
    if contains(cfg.controller_type, '5DOF')
        x_5    = x_est(1:5);
        y_curr = missile.C_5 * x_5 + missile.D_5 * u_prev;
        err    = ref(1:2) - y_curr;
        u_raw  = -state.K_x * x_5 - state.K_i * state.int_err;
    else
        y_curr = missile.C_track * x_est + missile.D_track * u_prev;
        err    = ref - y_curr;
        u_raw  = -state.K_x * x_est - state.K_i * state.int_err;
    end
    
    u_sat = max(min(u_raw, cfg.delta_max), -cfg.delta_max);
    state.int_err = state.int_err + (err + state.K_bc * (u_sat - u_raw)) * cfg.dt;
end
