function [f, J] = fbenney2dImplicit(domain, h, params)
    h = domain.reshapeToDomain(h);

    z = h.^4 / 4;

    Q = - (domain.diff(z, [2; 0]) + domain.diff(z, [0; 2])) / params.C;

    f = (domain.diff(Q, [2; 0]) + domain.diff(Q, [0; 2])) / 3;

    f = domain.reshapeToVector(f);

    if nargout == 2
        D2x = domain.diffMat([2; 0]);
        D2y = domain.diffMat([0; 2]);

        dQdz = - (D2x + D2y) .* domain.reshapeToVector(h).^3 / params.C;

        J = (D2x + D2y) * dQdz / 3;
    end
end
