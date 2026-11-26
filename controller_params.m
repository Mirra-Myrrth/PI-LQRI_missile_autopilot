%% РАСЧЕТ МАТРИЦ И КОЭФФИЦИЕНТОВ РЕГУЛЯТОРОВ (Базовый синтез vs Оптимизация ITAE)
ctrl_state = struct();
opts_itae = optimset('Display', 'off', 'MaxIter', 400, 'TolX', 1e-4, 'TolFun', 1e-4);
t_opt = 0:cfg.dt:1.5;

switch cfg.controller_type
    %% 1. LQRI (7-DOF) — ВЕСА ПО ПРАВИЛУ БРАЙСОНА
    case 'LQRI_7DOF'
        A_ext = [missile.A,       zeros(7, 3);
                -missile.C_track, zeros(3, 3)];
        B_ext = [missile.B;
                -missile.D_track];

        Q_bryson = diag([1/50^2, 1/50^2, 1/5^2, 1/5^2, 1/5^2, 0.1, 0.1, 100, 100, 50]);
        R_bryson = diag([1/cfg.delta_max^2, 1/cfg.delta_max^2, 1/cfg.delta_max^2]);

        K_ext = lqr(A_ext, B_ext, Q_bryson, R_bryson);
        ctrl_state.K_x     = K_ext(:, 1:7);
        ctrl_state.K_i     = K_ext(:, 8:10);
        ctrl_state.int_err = zeros(3, 1);
        ctrl_state.K_bc    = 5.0 * pinv(ctrl_state.K_i);

    %% 2. LQRI (7-DOF) — ОПТИМИЗАЦИЯ ВЕСОВ ПО КРИТЕРИЮ ITAE
    case 'LQRI_7DOF_ITAE'
        fprintf('  -> Оптимизация весов Q и R для LQRI_7DOF_ITAE...\n');
        A_ext = [missile.A,       zeros(7, 3);
                -missile.C_track, zeros(3, 3)];
        B_ext = [missile.B;
                -missile.D_track];
            
        %x0_7 = [0.01, 0.01, 0.1, 0.1, 0.1, 0.5, 0.5, 100, 100, 50, 10, 10, 10];
        % Агрессивный старт: высокий вес интеграторов (2500) и низкий штраф на рули (0.5)
        x0_7 = [0.01, 0.01, 0.5, 0.3, 0.3, 0.2, 0.2, 1500, 1500, 200, 2.0, 0.6, 0.6];

        ref_opt = [cfg.cmd_ny; cfg.cmd_nz; cfg.cmd_phi];
        
        [x_opt, ~] = fminsearch(@(w) obj_lqri_itae(w, missile.A, missile.B, missile.C_track, ...
                                missile.D_track, A_ext, B_ext, t_opt, ref_opt, cfg.delta_max), ...
                                x0_7, opts_itae);
                            
        Q_opt = diag(abs(x_opt(1:10))) + eye(10)*1e-4;
        R_opt = diag(abs(x_opt(11:13))) + eye(3)*1e-2;
        
        K_ext = lqrd(A_ext, B_ext, Q_opt, R_opt, cfg.dt);
        ctrl_state.K_x     = K_ext(:, 1:7);
        ctrl_state.K_i     = K_ext(:, 8:10);
        ctrl_state.int_err = zeros(3, 1);
        ctrl_state.K_bc    = 5.0 * pinv(ctrl_state.K_i);

    %% 3. LQRI (5-DOF) — ВЕСА ПО ПРАВИЛУ БРАЙСОНА
    case 'LQRI_5DOF'
        A_ext = [missile.A_5,  zeros(5, 2);
                -missile.C_5,  zeros(2, 2)];
        B_ext = [missile.B_5;
                -missile.D_5];

        Q_bryson = diag([1/50^2, 1/50^2, 1/5^2, 1/5^2, 1/5^2, 100, 100]);
        R_bryson = diag([1/cfg.delta_max^2, 1/cfg.delta_max^2, 1/cfg.delta_max^2]);

        K_ext = lqr(A_ext, B_ext, Q_bryson, R_bryson);
        ctrl_state.K_x     = K_ext(:, 1:5);
        ctrl_state.K_i     = K_ext(:, 6:7);
        ctrl_state.int_err = zeros(2, 1);
        ctrl_state.K_bc    = 5.0 * pinv(ctrl_state.K_i);

    %% 4. LQRI (5-DOF) — ОПТИМИЗАЦИЯ ВЕСОВ ПО КРИТЕРИЮ ITAE
    case 'LQRI_5DOF_ITAE'
        fprintf('  -> Оптимизация весов Q и R для LQRI_5DOF_ITAE...\n');
        A_ext = [missile.A_5,  zeros(5, 2);
                -missile.C_5,  zeros(2, 2)];
        B_ext = [missile.B_5;
                -missile.D_5];
            
        x0_5 = [0.01, 0.01, 0.1, 0.1, 0.1, 100, 100, 1, 1, 1];
        ref_opt = [cfg.cmd_ny; cfg.cmd_nz];
        
        [x_opt, ~] = fminsearch(@(w) obj_lqri_itae(w, missile.A_5, missile.B_5, missile.C_5, ...
                                missile.D_5, A_ext, B_ext, t_opt, ref_opt, cfg.delta_max), ...
                                x0_5, opts_itae);
                            
        Q_opt = diag(abs(x_opt(1:7))) + eye(7)*1e-4;
        R_opt = diag(abs(x_opt(8:10))) + eye(3)*1e-2;
        
        K_ext = lqr(A_ext, B_ext, Q_opt, R_opt);
        ctrl_state.K_x     = K_ext(:, 1:5);
        ctrl_state.K_i     = K_ext(:, 6:7);
        ctrl_state.int_err = zeros(2, 1);
        ctrl_state.K_bc    = 5.0 * pinv(ctrl_state.K_i);

  
  
   %% 5. ТРЁХКОНТУРНЫЙ ПИ-РЕГУЛЯТОР (ЧАСТОТНЫЙ / МОДАЛЬНЫЙ СИНТЕЗ) + АНТИ-ВИНДАП
    case 'PI_AW'
        % --- Канал Тангажа: 3 контура [w (угол атаки), q (ДУС), int_ny] -> delta_e ---
        idx_p = [2 4];
        A_p = missile.A_5(idx_p, idx_p);
        B_p = missile.B_5(idx_p, 2);
        C_ny = missile.C_5(1, idx_p);
        
        % 1) Внутренние два контура [Kw, Kq] стабилизируют (w, q) с полюсами [-12 + 8i, -12 - 8i]
        K_in_p = place(A_p, B_p, [-12 + 8i, -12 - 8i]);
        Kw_p   = K_in_p(1);
        Kq_p   = K_in_p(2);
        
        % 2) Внешний ПИ-контур по перегрузке n_y
        A_cl_p   = A_p - B_p * K_in_p;
        k_stat_p = C_ny * (-A_cl_p \ B_p); % Статический коэффициент передачи (~ +16.4 g/рад)
        Kp_p     = 0.4 / k_stat_p;
        Ki_p     = 5.0 / k_stat_p;

        % --- Канал Рыскания: 3 контура [v (скольжение), r (ДУС), int_nz] -> delta_r ---
        idx_y = [1 5];
        A_y = missile.A_5(idx_y, idx_y);
        B_y = missile.B_5(idx_y, 3);
        C_nz = missile.C_5(2, idx_y);
        
        K_in_y = place(A_y, B_y, [-12 + 8i, -12 - 8i]);
        Kv_y   = K_in_y(1);
        Kr_y   = K_in_y(2);
        
        A_cl_y   = A_y - B_y * K_in_y;
        k_stat_y = C_nz * (-A_cl_y \ B_y);
        Kp_y     = 0.4 / k_stat_y;
        Ki_y     = 5.0 / k_stat_y;

        % Сохраняем 4 коэффициента для каждого канала: [K_vel, K_rate, Kp, Ki]
        ctrl_state.pi_pitch = [Kw_p, Kq_p, Kp_p, Ki_p];
        ctrl_state.pi_yaw   = [Kv_y, Kr_y, Kp_y, Ki_y];
        ctrl_state.Kp_roll  = 0.05;
        ctrl_state.int_err  = zeros(2, 1);
        ctrl_state.Kb       = 5.0;

    %% 6. ТРЁХКОНТУРНЫЙ ПИ-РЕГУЛЯТОР (ОПТИМИЗАЦИЯ ПО КРИТЕРИЮ ITAE) + АНТИ-ВИНДАП
    case 'PI_AW_ITAE'
        fprintf('  -> Оптимизация 3-контурных ПИ-регуляторов по критерию ITAE...\n');
        idx_p = [2 4];
        A_p = missile.A_5(idx_p, idx_p); B_p = missile.B_5(idx_p, 2); C_ny = missile.C_5(1, idx_p);
        sys_p = ss(A_p, B_p, C_ny, missile.D_5(1, 2));
        
        K_in_p = place(A_p, B_p, [-12 + 8i, -12 - 8i]);
        k_stat_p = C_ny * (-(A_p - B_p * K_in_p) \ B_p);
        x0_p = [K_in_p(1), K_in_p(2), 0.4/k_stat_p, 5.0/k_stat_p];
        
        [pi_p, ~] = fminsearch(@(k) obj_pi3_itae(k, sys_p, t_opt, cfg.cmd_ny, cfg.delta_max), x0_p, opts_itae);
                           
        idx_y = [1 5];
        A_y = missile.A_5(idx_y, idx_y); B_y = missile.B_5(idx_y, 3); C_nz = missile.C_5(2, idx_y);
        sys_y = ss(A_y, B_y, C_nz, missile.D_5(2, 3));
        
        K_in_y = place(A_y, B_y, [-12 + 8i, -12 - 8i]);
        k_stat_y = C_nz * (-(A_y - B_y * K_in_y) \ B_y);
        x0_y = [K_in_y(1), K_in_y(2), 0.4/k_stat_y, 5.0/k_stat_y];
        
        [pi_y, ~] = fminsearch(@(k) obj_pi3_itae(k, sys_y, t_opt, cfg.cmd_nz, cfg.delta_max), x0_y, opts_itae);
                           
        ctrl_state.pi_pitch = pi_p;
        ctrl_state.pi_yaw   = pi_y;
        ctrl_state.Kp_roll  = 0.05;
        ctrl_state.int_err  = zeros(2, 1);
        ctrl_state.Kb       = 5.0;


    otherwise
        error('Неизвестный тип регулятора: %s', cfg.controller_type);
end

%% --- ЦЕЛЕВАЯ ФУНКЦИЯ ITAE ДЛЯ LQRI ---
function J = obj_lqri_itae(w, A, B, C, D, A_ext, B_ext, t, ref, u_max)
    nx = size(A, 1); ny = size(C, 1); nu = size(B, 2);
    try
        N = length(t); dt = t(2) - t(1);
        Q = diag(abs(w(1 : nx+ny))) + eye(nx+ny)*1e-4;
        R = diag(abs(w(nx+ny+1 : nx+ny+nu))) + eye(nu)*1e-2;
        K_ext = lqrd(A_ext, B_ext, Q, R, dt);
        
        if isempty(K_ext) || any(isnan(K_ext(:))), J = 1e12; return; end

        if isempty(K_ext), J = 1e12; return; end
        
        K_x = K_ext(:, 1:nx);
        K_i = K_ext(:, nx+1 : nx+ny);
        K_bc = 5.0 * pinv(K_i);
        
        N = length(t); dt = t(2) - t(1);
        x_st = zeros(nx, 1); x_i = zeros(ny, 1);
        y_hist = zeros(ny, N); u_hist = zeros(nu, N);
        
        for k = 1:N
            u_prev = u_hist(:, max(1, k-1));
            y_hist(:, k) = C * x_st + D * u_prev;
            err = ref - y_hist(:, k);
            
            u_raw = -K_x * x_st - K_i * x_i;
            u_sat = max(-u_max, min(u_max, u_raw));
            u_hist(:, k) = u_sat;
            
            x_i  = x_i + (err + K_bc * (u_sat - u_raw)) * dt;
            x_st = x_st + (A * x_st + B * u_sat) * dt;
        end
        
        %err_all = ref - y_hist;            
        %u_dot = diff(u_hist, 1, 2) / dt;
        %J = sum(t .* sum(abs(err_all), 1)) * dt + 0.001 * sum(u_dot(:).^2) * dt;
        
        err_all = ref - y_hist;
        u_dot = diff(u_hist, 1, 2) / dt;
        J = sum(t .* sum(abs(err_all), 1)) * dt + 0.002 * sum(u_dot(:).^2) * dt;
        
        if any(isnan(y_hist(:))) || max(abs(y_hist(:))) > 50, J = 1e12; end
        if any(isnan(y_hist(:))) || max(abs(y_hist(:))) > 50, J = 1e12; end
    catch
        J = 1e12;
    end
end


%% --- ЦЕЛЕВАЯ ФУНКЦИЯ ITAE ДЛЯ ТРЁХКОНТУРНОГО ПИ-РЕГУЛЯТОРА ---
function J = obj_pi3_itae(k_vec, sys, t, r_cmd, u_max)
    % k_vec = [K_vel, K_rate, Kp, Ki]
    N = length(t); dt = t(2) - t(1);
    x_st = zeros(2, 1); v_int = 0; K_aw = 5.0;
    y_hist = zeros(1, N); u_hist = zeros(1, N);
    
    for k = 1:N
        u_prev = u_hist(max(1, k-1));
        y_hist(k) = sys.C * x_st + sys.D * u_prev;
        err = r_cmd - y_hist(k);
        
        % Трёхконтурный закон: -K_vel*x(1) - K_rate*x(2) + Kp*err + int_err
        u_raw = -k_vec(1)*x_st(1) - k_vec(2)*x_st(2) + k_vec(3)*err + v_int;
        u_sat = max(-u_max, min(u_max, u_raw));
        u_hist(k) = u_sat;
        
        v_int = v_int + (k_vec(4)*err + K_aw*(u_sat - u_raw)) * dt;
        x_st  = x_st + (sys.A * x_st + sys.B * u_sat) * dt;
    end
    
    err_all = r_cmd - y_hist;
    u_dot = diff(u_hist) / dt;
    J = sum(t .* abs(err_all)) * dt + 0.001 * sum(u_dot.^2) * dt;
    if any(isnan(y_hist)) || max(abs(y_hist)) > 50, J = 1e12; end
end
 