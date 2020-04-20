function f = fbenney2dTransformed(domain, z, params)
    r2 = sqrt(2);

    P = 2*r2 * cot(params.theta) * z.^(1/4) - ...
        r2/(4*params.C) * (...
        z.^(-3/4) .* (domain.diff(z, [2; 0]) + domain.diff(z, [0; 2])) - ...
        3/4 * z.^(-7/4) .* (domain.diff(z, [1; 0]).^2 + domain.diff(z, [0; 1]).^2) ...
        );

    F1 = 4*r2/3 * z.^(3/4) - ...
        1/3 * 2*r2 * z.^(3/4) .* domain.diff(P, [1; 0]) + ...
        16*r2/15 * params.Re * z.^(3/4) .* domain.diff(z, [1; 0]);

    F2 = - 1/3 * 2*r2 * z.^(3/4) .* domain.diff(P, [0; 1]);

    f = - 2*r2 * z.^(3/4) .* (domain.diff(F1, [1; 0]) + domain.diff(F2, [0; 1]));
end
