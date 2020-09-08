function f = fhybrid(domain, y, params)
    n = length(y);
    h = domain.reshapeToDomain(y(1:n/2));
    F1 = domain.reshapeToDomain(y(1+n/2:n));

    P = 2 * h * cot(params.theta) ...
        - 1 / params.C * (domain.diff(h, [2; 0]) + domain.diff(h, [0; 2]));

    F2 = - domain.multiply(h, domain.diff(P, [0; 1]), [3, 1]) / 3;

    fh = - domain.diff(F1, [1; 0]) - domain.diff(F2, [0; 1]);

    h3 = domain.multiply(h, h, [2, 1]);

    fF1 = F1 - 2/3*h3 ...
        - params.epsilon * ( ...
        8 * params.Re / 15 * domain.multiply(h, domain.diff(h, [1; 0]), [6, 1]) ...
        - 1 / 3 * domain.multiply(h, domain.diff(P, [1; 0]), [3, 1]) ...
        ) ...
        - params.delta * ( ...
        18*params.Re/35*domain.multiply(domain.diff(h,[1;0]), F1,[1,2]) ...
        - 34*params.Re/35*domain.multiply(domain.multiply(h,F1),domain.diff(F1,[1;0])) ...
        - 1 / 3 * domain.multiply(h, domain.diff(P, [1; 0]), [3, 1]) ...
        );

    f = [domain.reshapeToVector(fh); domain.reshapeToVector(fF1)];

end
