function [f, J] = fbenney2dImplicit(domain, h, params)
    % P = - (domain.diff(h, [2; 0]) + domain.diff(h, [0; 2])) / params.C;

    % F1 = - 1 / 3 * domain.multiply(h, domain.diff(P, [1; 0]), [3, 1]);

    % F2 = - 1 / 3 * domain.multiply(h, domain.diff(P, [0; 1]), [3, 1]);

    % f = - domain.diff(F1, [1; 0]) - domain.diff(F2, [0; 1]);

    % z = h.^4 / 4;
    % r2 = sqrt(2);

    % dzdx = domain.diff(z, [1; 0]); 
    % dzdy = domain.diff(z, [0; 1]); 

    % P = - 1 / params.C * ( ...
    %     r2 / 4 * z.^(-3/4) .* (domain.diff(z, [2; 0]) + domain.diff(z, [0; 2])) - ...
    %     3 * r2 / 16 * z.^(-7/4) .* (dzdx.^2 + dzdy.^2));

    % F1 = - 2 * r2 / 3 * z.^(3/4) .* domain.diff(P, [1; 0]);

    % F2 = - 2 * r2 / 3 * z.^(3/4) .* domain.diff(P, [0; 1]);

    % f = - domain.diff(F1, [1; 0]) - domain.diff(F2, [0; 1]);

    % z = h.^4 / 4;
    % r2 = sqrt(2);

    % dzdx = domain.diff(z, [1; 0]); 
    % dzdy = domain.diff(z, [0; 1]); 

    % Q = - 1 / params.C * ( ...
    %     (domain.diff(z, [2; 0]) + domain.diff(z, [0; 2])) - ...
    %     3 / 4 * z.^(-1) .* (dzdx.^2 + dzdy.^2));

    % F1 = - 1 / 3 * domain.diff(Q, [1; 0]) + 1 / 4 * domain.diff(z, [1; 0]) ./ z .* Q;

    % F2 = - 1 / 3 * domain.diff(Q, [0; 1]) + 1 / 4 * domain.diff(z, [0; 1]) ./ z .* Q;

    % f = - domain.diff(F1, [1; 0]) - domain.diff(F2, [0; 1]);

    % z = h.^4 / 4;

    % Q = - (domain.diff(z, [2; 0]) + domain.diff(z, [0; 2])) / params.C;

    % F1 = - 1 / 3 * domain.diff(Q, [1; 0]);

    % F2 = - 1 / 3 * domain.diff(Q, [0; 1]);

    % f = - domain.diff(F1, [1; 0]) - domain.diff(F2, [0; 1]);

    z = h.^4 / 4;

    Q = - (domain.diff(z, [2; 0]) + domain.diff(z, [0; 2])) / params.C;

    f = (domain.diff(Q, [2; 0]) + domain.diff(Q, [0; 2])) / 3;

    if nargout == 2
        D2x = domain.diffMat([2; 0]);
        D2y = domain.diffMat([0; 2]);

        dQdz = - (D2x + D2y) .* domain.reshapeToVector(h).^3 / params.C;

        J = (D2x + D2y) * dQdz / 3;
    end
end
