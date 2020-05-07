function [dydt, J] = fwibl1Implicit(domain, y, params)
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

    if nargout == 2
        sd = @(v) spdiags(v,0,length(v),length(v));
        Dx = domain.diffMat([1; 0]);
        Dy = domain.diffMat([0; 1]);
        Dyy = domain.diffMat([0; 2]);
        
        DPDh = - 1 / params.C * (domain.diffMat([2; 0]) + domain.diffMat([0; 2]));

        DhtDh = sd(domain.reshapeToVector(h).^2 .* domain.reshapeToVector(Pyy))...
            * Dy + ...
            sd(domain.reshapeToVector(h).^3) * DPDh/3 ...
            * Dyy;

        DhtDF1 = - Dx;

        DF1tDh = -5/(6*params.Re) * (sd(domain.reshapeToVector(Px)) + ...
            sd(domain.reshapeToVector(h)) * DPDh * Dx);
        
        DF1tDF1 = zeros(numel(h));

        J = [DhtDh, DhtDF1; ...
            DF1tDh, DF1tDF1];
    end
end
