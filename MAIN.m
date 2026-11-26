clear; clc; close all;

%% 1. Подключение параметров ракеты (матрицы 7x7 и 7x3)
run('missile_params.m');

%% 2. Настройки моделирования (Конфигурация)
cfg.dt = 0.005;                    % Шаг интегрирования (с)
cfg.T_end = 2.0;                   % Время моделирования (с)

% Список регуляторов для последовательного запуска (можно оставить один или несколько)
cfg.controllers = {  'PI_AW', 'PI_AW_ITAE'};
%'LQRI_7DOF', 'LQRI_7DOF_ITAE', 'LQRI_5DOF', 'LQRI_5DOF_ITAE', 'PI_AW', 'PI_AW_ITAE'
cfg.use_kalman = true;             % Включить фильтр Калмана (true/false)
cfg.add_process_noise = true;      % Внешние возмущения (порывы ветра)
cfg.add_sensor_noise = true;       % Шумы датчиков

cfg.plot_individual = true;  % Строить отдельные окна для каждого регулятора (true/false)
cfg.plot_comparison = true;  % Строить сводные окна сравнения всех регуляторов (true/false)


% Заданные значения (уставки): [n_y (g); n_z (g); phi (рад)]
cfg.cmd_ny  = 3.5;                 
cfg.cmd_nz  = 1.5;                 
cfg.cmd_phi = 0.0;                 

% Ограничения отклонения рулей (~15 град) [рад]
cfg.delta_max = 0.26;

%% 3. Подготовка времени и общих профилей шумов
t = 0:cfg.dt:cfg.T_end;
N = length(t);

rng(42); 
% Вектор измерений датчиков (6 каналов): [n_y; n_z; p; q; r; phi]
% 1-2: акселерометры (g), 3-5: ДУС (рад/с), 6: датчик угла крена (рад)
sensor_noise_all = zeros(6, N);
if cfg.add_sensor_noise
    noise_std = [0.04; 0.04; 0.003; 0.003; 0.003; 0.002];
    sensor_noise_all = noise_std .* randn(6, N);
end

disturbance_all = zeros(7, N);
if cfg.add_process_noise
    disturbance_all(1, :) = 0.05 * (sin(3*t) + 0.2*randn(1, N));
    disturbance_all(2, :) = 0.05 * (cos(2*t) + 0.2*randn(1, N));
    disturbance_all(3:5, :) = 0.005 * randn(3, N);
end

%% 4. Цикл симуляции по всем выбранным регуляторам
for m = 1:length(cfg.controllers)
    cfg.controller_type = cfg.controllers{m};
    fprintf('Запуск симуляции для регулятора: %s...\n', cfg.controller_type);
    
    % Сброс массивов состояний перед каждым прогоном
    x      = zeros(7, N);     % Истинный вектор состояния [v; w; p; q; r; theta; phi]
    x_est  = zeros(7, N);     % Оценка состояния (наблюдатель Калмана)
    y_meas = zeros(6, N);     % Показания датчиков [n_y; n_z; p; q; r; phi]
    u_ctrl = zeros(3, N);     % Отклонения рулей [delta_a; delta_e; delta_r]
    y_out  = zeros(3, N);     % Истинные перегрузки и крен [n_y; n_z; phi]
    
    % Загрузка предрассчитанных коэффициентов и матриц текущего регулятора
    run('controller_params.m');
    
    % Инициализация наблюдателя Калмана (Акселерометры + ДУС + Датчик крена)
    kf_state = init_kf(missile, cfg);
    
    % --- Главный цикл интегрирования (РК4) ---
    for k = 1:N-1
        if k == 1
            u_prev = zeros(3, 1);
        else
            u_prev = u_ctrl(:, k-1);
        end
        
        % 4.1. Модель датчиков (6 каналов: [n_y; n_z; p; q; r; phi])
        y_true_sens  = kf_state.H * x(:, k) + kf_state.D * u_prev;
        y_meas(:, k) = y_true_sens + sensor_noise_all(:, k);
        
        % 4.2. Наблюдатель Калмана
        if cfg.use_kalman
            [x_est(:, k), kf_state] = step_kf(y_meas(:, k), u_prev, kf_state);
        else
            % Режим без фильтра Калмана: прямая подстановка зашумленных показаний ДУС и акселерометров
            x_est(:, k)    = x(:, k);
            x_est(1, k)    = (y_meas(2, k) - missile.D_track(2,3)*u_prev(3)) / missile.C_track(2,1);
            x_est(2, k)    = (y_meas(1, k) - missile.D_track(1,2)*u_prev(2)) / missile.C_track(1,2);
            x_est(3:5, k)  = y_meas(3:5, k); % Сырые показания ДУС [p; q; r]
            x_est(7, k)    = y_meas(6, k);   % Сырой угол крена phi
        end
        
        % 4.3. Расчет управляющего сигнала (Автопилот)
        ref = [cfg.cmd_ny; cfg.cmd_nz; cfg.cmd_phi];
        switch cfg.controller_type
            case {'PI_AW', 'PI_AW_ITAE'}
                [u_raw, ctrl_state] = ctrl_pi_aw(x_est(:, k), u_prev, ref, ctrl_state, missile, cfg);
            case {'LQRI_5DOF', 'LQRI_5DOF_ITAE', 'LQRI_7DOF', 'LQRI_7DOF_ITAE'}
                [u_raw, ctrl_state] = ctrl_lqri_bc(x_est(:, k), u_prev, ref, ctrl_state, missile, cfg);
            otherwise
                error('Неизвестный тип регулятора: %s', cfg.controller_type);
        end
        
        % Насыщение рулей (Saturation)
        u_ctrl(:, k) = max(min(u_raw, cfg.delta_max), -cfg.delta_max);
        
        % Сохранение истинных выходов для графиков
        y_out(:, k) = missile.C_track * x(:, k) + missile.D_track * u_ctrl(:, k);
        
        % 4.4. Внешние возмущения и интегрирование динамики ракеты (РК4)
        dist = disturbance_all(:, k);
        k1 = missile_dynamics(x(:, k), u_ctrl(:, k), dist, missile);
        k2 = missile_dynamics(x(:, k) + 0.5*cfg.dt*k1, u_ctrl(:, k), dist, missile);
        k3 = missile_dynamics(x(:, k) + 0.5*cfg.dt*k2, u_ctrl(:, k), dist, missile);
        k4 = missile_dynamics(x(:, k) + cfg.dt*k3, u_ctrl(:, k), dist, missile);
        
        x(:, k+1) = x(:, k) + (cfg.dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    
    y_meas(:, N) = kf_state.H * x(:, N) + kf_state.D * u_ctrl(:, N-1) + sensor_noise_all(:, N);
    x_est(:, N)  = x_est(:, N-1);
    u_ctrl(:, N) = u_ctrl(:, N-1);
    y_out(:, N)  = missile.C_track * x(:, N) + missile.D_track * u_ctrl(:, N);
    
    % Сохранение результатов текущего прогона для сводного сравнения
    results(m).type   = cfg.controller_type;
    results(m).x      = x;
    results(m).y_out  = y_out;
    results(m).u_ctrl = u_ctrl;

    % Построение индивидуальных графиков (если включено)
    if cfg.plot_individual
        plot_results(t, x, y_meas, x_est, u_ctrl, y_out, missile, cfg);
    end
end

%% 5. Построение сводных графиков сравнения регуляторов
    if cfg.plot_comparison && length(cfg.controllers) > 1
        plot_comparison_results(t, results, cfg);
    if cfg.plot_comparison && length(cfg.controllers) > 1
    plot_comparison_results(t, results, cfg);
end

print_quality_metrics(t, results, cfg);

fprintf('Все симуляции успешно завершены!\n');
    end

fprintf('Все симуляции успешно завершены!\n');

   



%% --- НАБЛЮДАТЕЛЬ КАЛМАНА (АКСЕЛЕРОМЕТРЫ + ДУС + ДАТЧИК КРЕНА) ---
function kf = init_kf(missile, cfg)
    % Формируем матрицы измерений для 6 датчиков: [n_y; n_z; p; q; r; phi]
    C_sens = [missile.C_track(1:2, :);    % Акселерометры: n_y, n_z
              0, 0, 1, 0, 0, 0, 0;       % ДУС крена: p
              0, 0, 0, 1, 0, 0, 0;       % ДУС тангажа: q
              0, 0, 0, 0, 1, 0, 0;       % ДУС рыскания: r
              0, 0, 0, 0, 0, 0, 1];      % Датчик угла крена: phi
          
    D_sens = [missile.D_track(1:2, :);
              zeros(4, 3)];

    % 1. Точная дискретизация полной 7-мерной модели методом ZOH
    sys_full_d = c2d(ss(missile.A, missile.B, C_sens, D_sens), cfg.dt);
    kf.Phi   = sys_full_d.A;
    kf.Gamma = sys_full_d.B;
    kf.H     = sys_full_d.C;
    kf.D     = sys_full_d.D;
    
    % 2. Синтез стационарного фильтра Калмана (dlqe) для 5-мерного динамического ядра
    % [v; w; p; q; r] (индексы 1 2 3 4 5) по 5 датчикам [n_y; n_z; p; q; r]
    idx_obs = [1 2 3 4 5];
    A_obs = missile.A(idx_obs, idx_obs);
    B_obs = missile.B(idx_obs, :);
    C_obs = C_sens(1:5, idx_obs);
    
    sys_obs_d = c2d(ss(A_obs, B_obs, C_obs, zeros(5, 3)), cfg.dt);
    
    % Ковариация шумов процесса (5x5) и шумов датчиков (2 акселерометра + 3 ДУС)
    Q_proc_obs = diag([1e-4, 1e-4, 1e-6, 1e-5, 1e-5]);
    R_sens_5   = diag([0.04^2, 0.04^2, 0.005^2, 0.005^2, 0.005^2]);

    
    L_obs = dlqe(sys_obs_d.A, eye(5), sys_obs_d.C, Q_proc_obs, R_sens_5);
    
    % 3. Расширяем матрицу усиления Калмана до 7x6
    kf.L = zeros(7, 6);
    kf.L(idx_obs, 1:5) = L_obs;
    kf.L(7, 6)         = 0.35; % Фильтрация угла крена phi по 6-му датчику
    
    kf.x = zeros(7, 1);
end

function [x_est, kf] = step_kf(y_meas, u_prev, kf)
    % Этап прогноза (Predict)
    x_pred = kf.Phi * kf.x + kf.Gamma * u_prev;
    
    % Прогнозируемые показания всех 6 датчиков
    y_pred = kf.H * x_pred + kf.D * u_prev;
    
    % Этап коррекции (Update)
    kf.x = x_pred + kf.L * (y_meas - y_pred);
    
    x_est = kf.x;
end

%% --- ПОСТРОЕНИЕ ГРАФИКОВ В ОТДЕЛЬНЫХ ОКНАХ (ПЕРЕГРУЗКИ + ДУС + РУЛИ + СПОСОБ НАСТРОЙКИ) ---
function plot_results(t, x, y_meas, x_est, u_ctrl, y_out, missile, cfg)
    % Вычисляем отфильтрованные перегрузки из восстановленного вектора x_est
    y_est_all = missile.C_track * x_est + missile.D_track * u_ctrl;
    
    % Формируем понятную подпись регулятора и способа его настройки
    full_label = get_tuning_label(cfg.controller_type);
    
    % Стиль для шумов (черный цвет и увеличенная толщина)
    noise_color = [1.0, 0.0, 0.0]; % Черный цвет
    noise_width = 1.2;

        %% ОКНО 1: КАНАЛ ТАНГАЖА (Зеленый сплошной для истинных значений)
    figure('Name', [full_label ' — Канал тангажа'], 'Color', 'w', 'NumberTitle', 'off');
    
    % 1. Нормальная перегрузка n_y
    subplot(3, 1, 1); hold on;
    if cfg.add_sensor_noise
        plot(t, y_meas(1,:), 'Color', noise_color, 'LineWidth', noise_width, 'DisplayName', 'Акселерометр n_y (шум)');
    end
    plot(t, y_out(1,:), 'g-', 'LineWidth', 1.6, 'DisplayName', 'Истинная n_y');
    if cfg.use_kalman
        plot(t, y_est_all(1,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Фильтр Калмана');
    end
    yline(cfg.cmd_ny, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Уставка n_{y,cmd}');
    title({full_label, 'Канал тангажа: Нормальная перегрузка n_y'}); 
    ylabel('Перегрузка (g)');
    legend('Location', 'best'); grid on;
    
    % 2. Угловая скорость тангажа q (показания ДУС и фильтрация)
    subplot(3, 1, 2); hold on;
    if cfg.add_sensor_noise
        plot(t, rad2deg(y_meas(4,:)), 'Color', noise_color, 'LineWidth', noise_width, 'DisplayName', 'ДУС q (шум)');
    end
    plot(t, rad2deg(x(4,:)), 'g-', 'LineWidth', 1.6, 'DisplayName', 'Истинная q');
    if cfg.use_kalman
        plot(t, rad2deg(x_est(4,:)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Фильтр Калмана');
    end
    title('Угловая скорость тангажа q (ДУС)');
    ylabel('Скорость (град/с)');
    legend('Location', 'best'); grid on;
    
    % 3. Отклонение руля высоты delta_e
    subplot(3, 1, 3); hold on;
    plot(t, rad2deg(u_ctrl(2,:)), 'm', 'LineWidth', 1.6, 'DisplayName', 'Руль высоты \delta_e');
    yline(rad2deg(cfg.delta_max), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Ограничение');
    yline(-rad2deg(cfg.delta_max), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    title('Отклонение руля высоты \delta_e'); 
    xlabel('Время (с)'); ylabel('Отклонение (град)');
    legend('Location', 'best'); grid on;

    %% ОКНО 2: КАНАЛ РЫСКАНИЯ (Синий сплошной для истинных значений)
    figure('Name', [full_label ' — Канал рыскания'], 'Color', 'w', 'NumberTitle', 'off');
    
    % 1. Боковая перегрузка n_z
    subplot(3, 1, 1); hold on;
    if cfg.add_sensor_noise
        plot(t, y_meas(2,:), 'Color', noise_color, 'LineWidth', noise_width, 'DisplayName', 'Акселерометр n_z (шум)');
    end
    plot(t, y_out(2,:), 'b-', 'LineWidth', 1.6, 'DisplayName', 'Истинная n_z');
    if cfg.use_kalman
        plot(t, y_est_all(2,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Фильтр Калмана');
    end
    yline(cfg.cmd_nz, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Уставка n_{z,cmd}');
    title({full_label, 'Канал рыскания: Боковая перегрузка n_z'}); 
    ylabel('Перегрузка (g)');
    legend('Location', 'best'); grid on;
    
    % 2. Угловая скорость рыскания r (показания ДУС и фильтрация)
    subplot(3, 1, 2); hold on;
    if cfg.add_sensor_noise
        plot(t, rad2deg(y_meas(5,:)), 'Color', noise_color, 'LineWidth', noise_width, 'DisplayName', 'ДУС r (шум)');
    end
    plot(t, rad2deg(x(5,:)), 'b-', 'LineWidth', 1.6, 'DisplayName', 'Истинная r');
    if cfg.use_kalman
        plot(t, rad2deg(x_est(5,:)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Фильтр Калмана');
    end
    title('Угловая скорость рыскания r (ДУС)');
    ylabel('Скорость (град/с)');
    legend('Location', 'best'); grid on;
    
    % 3. Отклонение руля направления delta_r
    subplot(3, 1, 3); hold on;
    plot(t, rad2deg(u_ctrl(3,:)), 'c', 'LineWidth', 1.6, 'DisplayName', 'Руль направления \delta_r');
    yline(rad2deg(cfg.delta_max), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Ограничение');
    yline(-rad2deg(cfg.delta_max), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    title('Отклонение руля направления \delta_r'); 
    xlabel('Время (с)'); ylabel('Отклонение (град)');
    legend('Location', 'best'); grid on;

end


%% --- РАСШИФРОВКА ТИПА РЕГУЛЯТОРА И СПОСОБА ЕГО НАСТРОЙКИ ---
function label = get_tuning_label(ctrl_type)
    switch ctrl_type
        case 'LQRI_7DOF'
            label = 'LQRI (7-DOF) | Настройка: Правило Брайсона';
        case 'LQRI_7DOF_ITAE'
            label = 'LQRI (7-DOF) | Настройка: Оптимизация ITAE (fminsearch)';
        case 'LQRI_5DOF'
            label = 'LQRI (5-DOF) | Настройка: Правило Брайсона';
        case 'LQRI_5DOF_ITAE'
            label = 'LQRI (5-DOF) | Настройка: Оптимизация ITAE (fminsearch)';
        case 'PI_AW'
            label = 'ПИ-регулятор + Anti-Windup | Настройка: Частотный синтез';
        case 'PI_AW_ITAE'
            label = 'ПИ-регулятор + Anti-Windup | Настройка: Оптимизация ITAE (fminsearch)';
        otherwise
            label = strrep(ctrl_type, '_', '\_');
    end
end

%% --- СВОДНЫЕ ГРАФИКИ СРАВНЕНИЯ РЕГУЛЯТОРОВ (ТАНГАЖ И РЫСКАНИЕ В ОТДЕЛЬНЫХ ОКНАХ) ---
function plot_comparison_results(t, results, cfg)
    num_ctrl = length(results);
    colors   = lines(num_ctrl); % Контрастная палитра для каждого регулятора
    styles   = {'-', '--', '-.', '-', '--', '-.'};

    %% ОКНО СРАВНЕНИЯ 1: КАНАЛ ТАНГАЖА
    figure('Name', 'Сравнение регуляторов — Канал тангажа', 'Color', 'w', 'NumberTitle', 'off');
    
    % 1. Нормальная перегрузка n_y
    subplot(3, 1, 1); hold on;
    for m = 1:num_ctrl
        lbl = get_tuning_label(results(m).type);
        st  = styles{mod(m-1, length(styles)) + 1};
        plot(t, results(m).y_out(1,:), 'Color', colors(m,:), 'LineStyle', st, ...
             'LineWidth', 1.6, 'DisplayName', lbl);
    end
    yline(cfg.cmd_ny, 'r:', 'LineWidth', 1.6, 'DisplayName', 'Уставка n_{y,cmd}');
    title('Сравнение регуляторов — Канал тангажа: Нормальная перегрузка n_y');
    ylabel('Перегрузка (g)');
    legend('Location', 'best'); grid on;

    % 2. Угловая скорость тангажа q
    subplot(3, 1, 2); hold on;
    for m = 1:num_ctrl
        lbl = get_tuning_label(results(m).type);
        st  = styles{mod(m-1, length(styles)) + 1};
        plot(t, rad2deg(results(m).x(4,:)), 'Color', colors(m,:), 'LineStyle', st, ...
             'LineWidth', 1.6, 'DisplayName', lbl);
    end
    title('Сравнение регуляторов — Угловая скорость тангажа q');
    ylabel('Скорость (град/с)');
    legend('Location', 'best'); grid on;

    % 3. Отклонение руля высоты delta_e
    subplot(3, 1, 3); hold on;
    for m = 1:num_ctrl
        lbl = get_tuning_label(results(m).type);
        st  = styles{mod(m-1, length(styles)) + 1};
        plot(t, rad2deg(results(m).u_ctrl(2,:)), 'Color', colors(m,:), 'LineStyle', st, ...
             'LineWidth', 1.6, 'DisplayName', lbl);
    end
    yline(rad2deg(cfg.delta_max), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Ограничение');
    yline(-rad2deg(cfg.delta_max), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    title('Сравнение регуляторов — Отклонение руля высоты \delta_e');
    xlabel('Время (с)'); ylabel('Отклонение (град)');
    legend('Location', 'best'); grid on;

    %% ОКНО СРАВНЕНИЯ 2: КАНАЛ РЫСКАНИЯ
    figure('Name', 'Сравнение регуляторов — Канал рыскания', 'Color', 'w', 'NumberTitle', 'off');
    
    % 1. Боковая перегрузка n_z
    subplot(3, 1, 1); hold on;
    for m = 1:num_ctrl
        lbl = get_tuning_label(results(m).type);
        st  = styles{mod(m-1, length(styles)) + 1};
        plot(t, results(m).y_out(2,:), 'Color', colors(m,:), 'LineStyle', st, ...
             'LineWidth', 1.6, 'DisplayName', lbl);
    end
    yline(cfg.cmd_nz, 'r:', 'LineWidth', 1.6, 'DisplayName', 'Уставка n_{z,cmd}');
    title('Сравнение регуляторов — Канал рыскания: Боковая перегрузка n_z');
    ylabel('Перегрузка (g)');
    legend('Location', 'best'); grid on;

    % 2. Угловая скорость рыскания r
    subplot(3, 1, 2); hold on;
    for m = 1:num_ctrl
        lbl = get_tuning_label(results(m).type);
        st  = styles{mod(m-1, length(styles)) + 1};
        plot(t, rad2deg(results(m).x(5,:)), 'Color', colors(m,:), 'LineStyle', st, ...
             'LineWidth', 1.6, 'DisplayName', lbl);
    end
    title('Сравнение регуляторов — Угловая скорость рыскания r');
    ylabel('Скорость (град/с)');
    legend('Location', 'best'); grid on;

    % 3. Отклонение руля направления delta_r
    subplot(3, 1, 3); hold on;
    for m = 1:num_ctrl
        lbl = get_tuning_label(results(m).type);
        st  = styles{mod(m-1, length(styles)) + 1};
        plot(t, rad2deg(results(m).u_ctrl(3,:)), 'Color', colors(m,:), 'LineStyle', st, ...
             'LineWidth', 1.6, 'DisplayName', lbl);
    end
    yline(rad2deg(cfg.delta_max), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Ограничение');
    yline(-rad2deg(cfg.delta_max), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    title('Сравнение регуляторов — Отклонение руля направления \delta_r');
    xlabel('Время (с)'); ylabel('Отклонение (град)');
    legend('Location', 'best'); grid on;
end
%% --- РАСЧЕТ И ВЫВОД ТАБЛИЦЫ ПОКАЗАТЕЛЕЙ КАЧЕСТВА (ВРЕМЯ РЕГУЛИРОВАНИЯ И ПЕРЕРЕГУЛИРОВАНИЕ) ---
function print_quality_metrics(t, results, cfg)
    fprintf('\n=========================================================================================================\n');
    fprintf(' СРАВНЕНИЕ КАЧЕСТВА ПЕРЕХОДНЫХ ПРОЦЕССОВ (Трубка точности: +-2%% от уставки)\n');
    fprintf('=========================================================================================================\n');
    fprintf(' %-56s | %10s | %10s | %10s | %10s\n', ...
            'Регулятор и способ настройки', 't_рег ny, с', 'Перер ny,%', 't_рег nz, с', 'Перер nz,%');
    fprintf('---------------------------------------------------------------------------------------------------------\n');

    for m = 1:length(results)
        lbl = get_tuning_label(results(m).type);
        
        % Метрики по каналу тангажа (n_y) и рыскания (n_z)
        [ts_ny, os_ny] = calc_step_metrics(t, results(m).y_out(1, :), cfg.cmd_ny);
        [ts_nz, os_nz] = calc_step_metrics(t, results(m).y_out(2, :), cfg.cmd_nz);
        
        % Форматирование строк (если процесс не устоялся за время T_end, выводим '> T_end')
        str_ts_ny = format_ts(ts_ny, cfg.T_end);
        str_ts_nz = format_ts(ts_nz, cfg.T_end);
        
        fprintf(' %-56s | %10s | %9.2f%% | %10s | %9.2f%%\n', ...
                lbl, str_ts_ny, os_ny, str_ts_nz, os_nz);
    end
    fprintf('=========================================================================================================\n\n');
end

function [t_settle, overshoot] = calc_step_metrics(t, y, cmd)
    if abs(cmd) < 1e-6
        t_settle = NaN;
        overshoot = 0;
        return;
    end
    
    % 1. Перерегулирование (%)
    y_aligned = y * sign(cmd);
    peak_val  = max(y_aligned);
    overshoot = max(0, (peak_val - abs(cmd)) / abs(cmd) * 100);
    
    % 2. Время регулирования по зоне +-2% от установившегося (заданного) значения
    tol_band = 0.02 * abs(cmd);
    idx_out  = find(abs(y - cmd) > tol_band, 1, 'last');
    
    if isempty(idx_out)
        t_settle = 0;
    elseif idx_out == length(t)
        t_settle = NaN; % Не вошел в трубку 2% к концу симуляции
    else
        t_settle = t(idx_out + 1);
    end
end

function str = format_ts(ts, t_end)
    if isnan(ts)
        str = sprintf('> %.2f', t_end);
    else
        str = sprintf('%.3f', ts);
    end
end
