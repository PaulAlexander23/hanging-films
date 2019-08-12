function J = jbenney2d(domain, u, params)
    Dx = domain.getDiffMatrix([1; 0]) ;
    Dy = domain.getDiffMatrix([0; 1]);
    Lap = domain.getDiffMatrix([2; 0]) + domain.getDiffMatrix([0; 2]);
    Lap_Dx = Dx * Lap;
    Lap_Dy = Dy * Lap;
    
    U = domain.reshapeToVector(u);
    dUdx = domain.reshapeToVector(domain.diff(u, [1; 0]));
    dUdy = domain.reshapeToVector(domain.diff(u, [0; 1]));
    
    n = numel(u);
    
    q1_h = spdiags(2 * U.^2 .* (1 - dUdx * cot(params.theta) + ...
        (Lap_Dx * U) / (2 * params.C)), 0, n, n) + ...
        2 * U.^3/3 .* (-cot(params.theta) * Dx + Lap_Dx / (2 * params.C));
    
    q1_h = q1_h + ...
        params.Re * 8 * U.^6/15 .* Dx + ...
        spdiags(params.Re * 16 * U.^5 .* dUdx / 5, 0, n, n);
    
    q2_h = spdiags(2 * U.^2 .* (-dUdy * cot(params.theta) + ...
        (Lap_Dy * U) / (2 * params.C)), 0, n, n) + ...
        2 * U.^3/3 .* (-cot(params.theta) * Dy + Lap_Dy / (2 * params.C));
    
    J = -(Dx * q1_h + Dy * q2_h);
end