function f = fbenney2dSuppressErrorAfterEval(domain, y, params)
    y = domain.reshapeToDomain(y);

    P = 2 * cot(params.theta) * y - ...
        (domain.diff(y, [2; 0]) + domain.diff(y, [0; 2])) / params.C;

    F1 = 2 / 3 * domain.multiply(y, y, [2, 1]) - ...
        1 / 3 * domain.multiply(y, domain.diff(P, [1; 0]), [3, 1]) + ...
        8 * params.Re / 15 * domain.multiply(y, domain.diff(y, [1; 0]), [6, 1]);

    F2 = - 1 / 3 * domain.multiply(y, domain.diff(P, [0; 1]), [3, 1]);

    f = - domain.diff(F1, [1; 0]) - ...
        domain.diff(F2, [0; 1]);

    f = suppressSmallModes(f, 1e-7);

    f = domain.reshapeToVector(f);
end
