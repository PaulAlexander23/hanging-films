function J = jhybrid(domain, y, params)
    n = length(y);
    h = domain.reshapeToDomain(y(1:n/2));
    F1 = domain.reshapeToDomain(y(1+n/2:n));

    P = 2 * h * cot(params.theta) ...
        - 1 / params.C * (domain.diff(h, [2; 0]) + domain.diff(h, [0; 2]));

    h3 = domain.multiply(h, h, [2, 1]);

    N = n/2;

    sv = @(in) toSparseVector(domain, in, N);

    hV = sv(h);

    dPdh = 2*cot(params.theta)*speye(N) ...
        - (domain.diffMat([2, 0]') + domain.diffMat([0, 2]'))/params.C;
    dF2dh = - hV^2 * sv(domain.diff(P, [0, 1]')) ...
        - hV^3 * dPdh * domain.diffMat([0, 1]') / 3;

    dQ1dh = - domain.diffMat([0; 1]) * dF2dh;
    dQ1dF1 = - domain.diffMat([1, 0]');

    Jh = [dQ1dh, dQ1dF1];

    fF1 = domain.multiply( ...
        2/3*h3 - F1 ...
        + params.epsilon * ( ...
        8 * params.Re / 15 * domain.multiply(h, domain.diff(h, [1; 0]), [6, 1]) ...
        ) ...
        + params.delta * ( ...
        18*params.Re/35*domain.multiply(domain.diff(h,[1;0]), F1,[1,2]) ...
        - 34*params.Re/35*domain.multiply(domain.multiply(h,F1),domain.diff(F1,[1;0])) ...
        ) ...
        + (params.epsilon + params.delta) * ( ...
        - 1 / 3 * domain.multiply(h, domain.diff(P, [1; 0]), [3, 1]) ...
        ) ...
        , h, [1, -2]) / (params.delta * 2 / 5 * params.Re);

    dQ2dh = sv(h.^-2) * ( ...
        2 * sv(domain.multiply(h, h)) ...
        + params.epsilon * ( ...
        8 * params.Re / 15 * (6 * sv(domain.multiply(h, domain.diff(h, [1; 0]), [5, 1])) ...
        + domain.diffMat([1;0]) * sv(domain.multiply(h, h, [5, 1])) ) ...
        ) ...
        + params.delta * ( ...
        18*params.Re/35 * domain.diffMat([1;0]) * sv(domain.multiply(F1, F1)) ...
        - 34*params.Re/35 * sv(domain.multiply(F1,domain.diff(F1,[1;0]))) ...
        ) ...
        + (params.epsilon + params.delta) * ( ...
        - sv(domain.multiply(h, domain.diff(P, [1; 0]), [2, 1])) ...
        - sv(h3) * domain.diffMat([1;0]) * dPdh / 3 ...
        ) ...
        ) / (params.delta * 2 / 5 * params.Re) ...
        - 2 * sv(h.^-1) * sv(fF1);

    dQ2dF1 = sv(h.^-2) * ( ...
        - speye(N) ...
        + params.delta * ( ...
        36*params.Re/35 * sv(domain.multiply(domain.diff(h,[1;0]), F1)) ...
        - 34*params.Re/35 * sv(domain.multiply(h,domain.diff(F1,[1;0]))) ...
        - 34*params.Re/35 * sv(domain.multiply(h,F1)) * domain.diffMat([1;0]) ...
        ) ...
        ) / (params.delta * 2 / 5 * params.Re);

    JF1 = [dQ2dh, dQ2dF1];
    
    J = [Jh; JF1];

    function sparseVector = toSparseVector(domain, in, N)
        vector = domain.reshapeToVector(in);
        sparseVector = spdiags(vector, 0, N, N);
    end
end
