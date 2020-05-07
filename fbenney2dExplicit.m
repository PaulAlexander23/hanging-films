function f = fbenney2dExplicit(domain, h, params)
    h = domain.reshapeToDomain(h);

    z = h.^4 / 4;

    dzdx = domain.diff(z, [1; 0]); 
    dzdy = domain.diff(z, [0; 1]); 

    P = 2 * cot(params.theta) * h;

    R = 3 * (dzdx.^2 + dzdy.^2) ...
        ./(4 * params.C * z);
    
    Q = - 1 / params.C * (domain.diff(z, [2; 0]) + domain.diff(z, [0; 2])) ...
        + R;
    
    F1 = 2 / 3 * domain.multiply(h, h, [2, 1]) - ...
        1 / 3 * domain.multiply(h, domain.diff(P, [1; 0]), [3, 1]) + ...
        8 * params.Re / 15 * domain.multiply(h, domain.diff(h, [1; 0]), [6, 1]) ...
        - 1 / 3 * domain.diff(R, [1; 0]) ...
        + 1 / 4 * dzdx ./ z .* Q;

    F2 = - 1 / 3 * domain.multiply(h, domain.diff(P, [0; 1]), [3, 1]) ...
        - 1 / 3 * domain.diff(R, [0; 1]) ...
        + 1 / 4 * dzdy ./ z .* Q;

    f = - domain.diff(F1, [1; 0]) - domain.diff(F2, [0; 1]);
    
    f = domain.reshapeToVector(f);
end
