%% ПАРАМЕТРЫ РАКЕТЫ ASTRA MK-1 (M=3.0, H=7000м)
missile.Mach = 3.0; 
missile.H    = 7000;    
T_curr   = 288.15 - 0.0065 * missile.H; 
P_curr   = 101325 * (1 - 0.0065 * missile.H / 288.15)^5.25588; 
rho_curr = P_curr / (287.05 * T_curr); 
a_sound  = sqrt(1.4 * 287.05 * T_curr); 

missile.V_M   = missile.Mach * a_sound; 
missile.q_inf = 0.5 * rho_curr * missile.V_M^2;

% Геометрия и масса Astra Mk-1
missile.m       = 154.0; 
missile.g       = 9.81; 
d_num          = 0.178; 
missile.l_total = 3.63;             
R              = d_num / 2; 
missile.S       = pi * R^2;              
missile.I_a     = 0.5 * missile.m * R^2;                          
missile.I_t     = (1/12) * missile.m * (3 * R^2 + missile.l_total^2); 

% Аэродинамические коэффициенты Astra Mk-1
C_M_de    =  1.85;  
C_L_da    =  0.35; 
C_Y_dr    = -0.45;  
C_N_dr    = -C_M_de; 
C_M_alpha = -1.95;  
C_N_beta  =  1.95; 
C_L_alpha =  3.80;  
C_Y_beta  = -3.80;   

% Вспомогательные аэродинамические комплексы
qS_mV   = (missile.q_inf * missile.S) / (missile.m * missile.V_M);
qS_m    = (missile.q_inf * missile.S) / missile.m;
qSl_Ia  = (missile.q_inf * missile.S * missile.l_total) / missile.I_a;
qSl_It  = (missile.q_inf * missile.S * missile.l_total) / missile.I_t;
qSl_ItV = qSl_It / missile.V_M;
k_g     = (missile.q_inf * missile.S) / (missile.m * missile.g);

%% 1. МАТРИЦЫ 5-DOF МОДЕЛИ: x_5 = [v; w; p; q; r], u = [delta_a; delta_e; delta_r]
% Связанная система координат (X вперед, Y вправо, Z вниз):
% При кабрировании (q > 0) проекция скорости w (ось Z вниз) уменьшается (-V_M * q),
% а нормальная перегрузка n_y (направленная вверх) положительна при w < 0 (-k_g * C_L_alpha * w / V_M).
missile.A_5 = [
    qS_mV * C_Y_beta,   0,                 0,  0,           -missile.V_M;
    0,                 -qS_mV * C_L_alpha, 0, -missile.V_M,   0;
    0,                  0,                 0,  0,            0;
    0,                  qSl_ItV * C_M_alpha, 0, 0,           0;
    qSl_ItV * C_N_beta, 0,                 0,  0,            0
];

missile.B_5 = [
    0,               0,              -qS_m * C_Y_dr;
    0,               0,               0;
    qSl_Ia * C_L_da, 0,               0;
    0,               qSl_It * C_M_de, 0;
    0,               0,               qSl_It * C_N_dr
];

% Выходы перегрузок: y_5 = [n_y (верт. вверх); n_z (бок. вправо)]
missile.C_5 = [
    0,                     -k_g * (C_L_alpha / missile.V_M), 0, 0, 0;
    k_g * (C_Y_beta / missile.V_M), 0,                       0, 0, 0
];

missile.D_5 = [
    0, 0, 0;
    0, 0, k_g * C_Y_dr
];

%% 2. МАТРИЦЫ 7-DOF МОДЕЛИ: x_7 = [v; w; p; q; r; theta; phi]
missile.theta_0 = (missile.m * missile.g) / (missile.q_inf * missile.S * C_L_alpha);

missile.A = [
    missile.A_5, [0, missile.g * cos(missile.theta_0);
                -missile.g * sin(missile.theta_0), 0;
                 zeros(3, 2)];
    0, 0, 0, 1, 0,                      0, 0;
    0, 0, 1, 0, tan(missile.theta_0),    0, 0
];

missile.B = [missile.B_5; zeros(2, 3)];

% Выходы 7-DOF модели: y_7 = [n_y; n_z; phi]
missile.C_track = [missile.C_5, zeros(2, 2);
                  0, 0, 0, 0, 0, 0, 1];

missile.D_track = [missile.D_5;
                  0, 0, 0];
