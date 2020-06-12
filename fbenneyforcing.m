function f = fbenneyforcing(t, domain, h0, params)
    % h0 = domain.reshapeToDomain(h0);

    % P = 2 * cot(params.theta) * h0 ...
    %     - (domain.diff(h0, [2;0]) + domain.diff(h0, [0;2]))/params.C;

    % F1 = 2/3 * h0.^3 ...
    %     - 1/3 * h0.^3 .* domain.diff(P, [1;0]) ...
    %     + 8 * params.Re / 15 * h0.^6 .* domain.diff(h0, [1;0]);

    % F2 = - 1/3 * h0.^3 .* domain.diff(P, [0;1]);

    % f = params.c * params.a * sin(domain.x{1} - params.c * t) + domain.diff(F1, [1;0]) + domain.diff(F2, [0;1]);

    x = domain.x{1};
    y = domain.x{2};

    h0 = 1 - params.a * cos(x) - params.b * cos(y);

    sx = params.a * sin(x - params.c*t) + 0 * y;
    cx = params.a * cos(x - params.c*t) + 0 * y;
    sy = params.b * sin(y) + 0 * x;
    cy = params.b * cos(y) + 0 * x;

    F1x = 2 * h0.^2 .* sx ...
        - h0.^2 .* sx.^2 * (2 * cot(params.theta) + 1 / params.C) ...
        - h0.^3 / 3 .* cx * (2 * cot(params.theta) + 1 / params.C) ...
        + 16 * params.Re / 5 * h0.^5 .* sx.^2 ...
        + 8 * params.Re / 15 * h0.^6 .* cx;

    F2y = - h0.^2 .* sy.^2 * (2 * cot(params.theta) + 1 / params.C) ...
        - h0.^3 / 3 .* cy * (2 * cot(params.theta) + 1 / params.C);

    f = - params.c * sx + F1x + F2y;

    f = domain.reshapeToVector(f);
end
