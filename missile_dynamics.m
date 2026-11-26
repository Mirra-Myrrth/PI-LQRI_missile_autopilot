function dx = missile_dynamics(x, u, dist, missile)
    % x = [v; w; p; q; r; theta; phi] (7x1)
    % u = [delta_a; delta_e; delta_r] (3x1)
    dx = missile.A * x + missile.B * u + dist;
end
