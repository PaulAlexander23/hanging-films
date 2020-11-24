function dYdt = fwibl1(domain, Y, params)
    Y = domain.reshapeToDomain(Y);

    h = Y(1:end/2, :, :);
    F1 = Y(end/2+1:end, :, :);
    
    P = 2 * h * cot(params.theta) ...
        - 1 / params.C * (domain.diff(h, [2; 0]) + domain.diff(h, [0; 2]));

    F2 = - domain.multiply(h, domain.diff(P, [0; 1]), [3, 1]) / 3;

    dydt = - domain.diff(F1, [1; 0]) - domain.diff(F2, [0; 1]);

    dF1dt = 5 / (3 * params.Re) * h ...
        - 5 / (6 * params.Re) * domain.multiply(h, domain.diff(P, [1; 0])) ...
        - 5 / (2 * params.Re) * domain.multiply(F1, h, [1, -2]) ...
        + 9 / 7 * domain.multiply(domain.diff(h, [1; 0]), ...
        domain.multiply(F1, h, [2, -2])) ...
        - 17 / 7 * domain.multiply(domain.diff(F1, [1; 0]), ...
        domain.multiply(F1, h, [1, -1]));
        
    dYdt = cat(1, dydt, dF1dt);

    dYdt = domain.reshapeToVector(dYdt);
end

