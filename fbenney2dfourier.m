function f = fbenney2dfourier(domain, y, params)
    delta = params(1);
    theta = params(2);
    Re = params(3);
    C = params(4);

    P = 2 * cot(theta) * y - ...
        (domain.diff(y, [2; 0]) + domain.diff(y, [0; 2])) / C;

    F1 = 2 / 3 * domain.multiply(y, y, [2, 1]) + 8 * Re / 15 * domain.multiply(y, domain.diff(y, [1; 0]), [6, 1]) - 1 / 3 * domain.multiply(y, domain.diff(P, [1; 0]), [3, 1]);
    F2 = - 1 / 3 * domain.multiply(y, domain.diff(P, [0; 1]), [3, 1]);

    f = -domain.diff(F1, [1; 0]) - ...
        domain.diff(F2, [0; 1]);
end
