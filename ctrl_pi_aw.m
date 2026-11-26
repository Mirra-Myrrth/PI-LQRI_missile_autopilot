function [u_raw, state] = ctrl_pi_aw(x_est, u_prev, ref, state, missile, cfg)
    % Текущие перегрузки: y_curr = [n_y; n_z; phi]
    y_curr = missile.C_track * x_est + missile.D_track * u_prev;
    err    = ref(1:2) - y_curr(1:2); % [err_ny; err_nz]
    
    % Оценки линейных и угловых скоростей из фильтра Калмана
    v = x_est(1); w = x_est(2);
    p = x_est(3); q = x_est(4); r = x_est(5);
    
    % Коэффициенты [K_vel, K_rate, Kp, Ki]
    kp = state.pi_pitch; 
    ky = state.pi_yaw;   
    
    % 1. Канал крена (демпфирование)
    u_da = -state.Kp_roll * p;
    
    % 2. Трёхконтурный канал тангажа: w (угол атаки) + q (ДУС) + ПИ по перегрузке n_y
    u_de = -kp(1)*w - kp(2)*q + kp(3)*err(1) + state.int_err(1);
    
    % 3. Трёхконтурный канал рыскания: v (скольжение) + r (ДУС) + ПИ по перегрузке n_z
    u_dr = -ky(1)*v - ky(2)*r + ky(3)*err(2) + state.int_err(2);
    
    u_raw = [u_da; u_de; u_dr];
    u_sat = max(min(u_raw, cfg.delta_max), -cfg.delta_max);
    
    % Анти-виндап (Back-Calculation) для интеграторов тангажа и рыскания
    sat_diff = [u_sat(2) - u_raw(2); u_sat(3) - u_raw(3)];
    Ki_vec   = [kp(4); ky(4)];
    
    state.int_err = state.int_err + (Ki_vec .* err + state.Kb * sat_diff) * cfg.dt;
end

