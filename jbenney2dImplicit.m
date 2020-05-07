function J = jbenney2dImplicit(domain, h, params)
    D2x = domain.diffMat([2; 0]);
    D2y = domain.diffMat([0; 2]);

    dQdz = - (D2x + D2y) .* h.^3 / params.C;

    J = (D2x + D2y) * dQdz / 3;
end
