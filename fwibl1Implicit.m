function dydt = fwibl1Implicit(domain, y, params)
    y = domain.reshapeToDomain(y);

    h = y(1:end/2, :, :);
    F1 = y(end/2+1:end, :, :);
    
    P = - 1 / params.C * (domain.diff(h, [2; 0]) + domain.diff(h, [0; 2]));

    Px = domain.diff(P, [1; 0]);

    Pyy = domain.diff(P, [0; 2]);

    dydt = - domain.diff(F1, [1; 0]) ...
        + domain.multiply(h, Pyy, [3, 1])/3;

    dF1dt = - 5 * domain.multiply(h, Px) / ...
        (6 * params.Re);

    dydt = cat(1, dydt, dF1dt);

    dydt = domain.reshapeToVector(dydt);
end
