function dydt = fwibl1Explicit(domain, y, params)
    y = domain.reshapeToDomain(y);

    h = y(1:end/2, :, :);
    F1 = y(end/2+1:end, :, :);
    
    P = 2 * h * cot(params.theta);

    F2 = - domain.multiply(h, domain.diff(P, [0; 1]), [3, 1]) / 3;

    dhdt = - domain.diff(F2, [0; 1]) ...
        - 1 / params.C * domain.multiply( ...
        domain.multiply(h, domain.diff(h, [0; 1]), [2, 1]), ...
        domain.diff(domain.diff(h, [2; 0]) + domain.diff(h, [0; 2]), [0; 1]));

    dF1dt =  domain.multiply( ...
        9 * domain.multiply( ...
        domain.multiply(domain.diff(h, [1; 0]), F1, [1, 2]), ...
        h, [1, -1]) - 17 * domain.multiply(F1, domain.diff(F1, [1; 0]), [1, 1]), ...
        (7 * h), [1, -1]) + ...
        (10 * h - 15 * domain.multiply(F1, h, [1, -2]) - ...
        5 * domain.multiply(h, domain.diff(P, [1; 0]), [1, 1])) / (6 * params.Re);

    dydt = cat(1, dhdt, dF1dt);

    dydt = domain.reshapeToVector(dydt);
end
