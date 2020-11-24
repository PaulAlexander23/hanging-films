function f = fbenney2dImplicit(domain, h, params)
    h = domain.reshapeToDomain(h);

    P = - (domain.diff(h, [2; 0]) + domain.diff(h, [0; 2])) / params.C;

    f = domain.multiply(h, ...
        domain.diff(P, [2; 0]) + domain.diff(P, [0; 2]), ...
        [3, 1]) / 3;

    f = domain.reshapeToVector(f);
end
