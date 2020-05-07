function f = fbenney2dImplicit(domain, h, params)
    h = domain.reshapeToDomain(h);

    z = h.^4 / 4;

    Q = - (domain.diff(z, [2; 0]) + domain.diff(z, [0; 2])) / params.C;

    f = (domain.diff(Q, [2; 0]) + domain.diff(Q, [0; 2])) / 3;

    f = domain.reshapeToVector(f);
end
