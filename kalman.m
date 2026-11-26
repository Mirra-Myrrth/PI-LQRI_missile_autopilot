function kf = init_kalman(missile, cfg)
    kf.A_d = eye(6) + missile.A * cfg.dt;
    kf.B_d = missile.B * cfg.dt;
    kf.H   = eye(6);
    kf.Q   = diag([1e-4, 1e-3, 1e-4, 1e-4, 1e-3, 1e-4]); % Шум процесса
    kf.R   = diag([0.005, 0.002, 0.005, 0.005, 0.002, 0.005].^2); % Шум датчиков
    kf.P   = eye(6);
    kf.x   = zeros(6, 1);
end

function [x_est, kf] = kalman_filter(y, u, kf, ~, ~)
    % Этап предсказания (Predict)
    x_pred = kf.A_d * kf.x + kf.B_d * u;
    P_pred = kf.A_d * kf.P * kf.A_d' + kf.Q;
    
    % Этап коррекции (Update)
    K_gain = P_pred * kf.H' / (kf.H * P_pred * kf.H' + kf.R);
    kf.x = x_pred + K_gain * (y - kf.H * x_pred);
    kf.P = (eye(6) - K_gain * kf.H) * P_pred;
    
    x_est = kf.x;
end