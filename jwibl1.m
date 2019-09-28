function jacobian = jwibl1(domain, y, params)
    h = y(1:end/2, :, :);
    F1 = y(end/2+1:end, :, :);
    
    N = numel(h(:,:,1));
    sv = @(in) toSparseVector(domain, in, N);
    
    hV = sv(h);
    F1V = sv(F1);
    
    P = 2*cot(params.theta)*h - (domain.diff(h, [2, 0]') + domain.diff(h, [0, 2]'))/params.C;
    dPdh = 2*cot(params.theta)*speye(N) - (domain.diffMat([2, 0]') + domain.diffMat([0, 2]'))/params.C;
    dF2dh = - hV^2 * sv(domain.diff(P, [0, 1]')) - hV^3 * dPdh * domain.diffMat([0, 1]') / 3;
    
    dQ1dh = - domain.diffMat([0, 1]') * dF2dh;
    dQ1dF1 = - domain.diffMat([1, 0]');
    
    R2 = - F1V + 2/3 * hV.^3 + 18/35 * params.Re * sv(domain.diff(h, [1, 0]')) * F1V.^2 - 34/35 * params.Re * hV * F1V * sv(domain.diff(F1, [1, 0]')) - 1/3 * hV.^3 * sv(domain.diff(P, [1, 0]'));
    dR2dh = 2 * hV.^2 + 18/35 * params.Re * F1V.^2  * domain.diffMat([1, 0]') - 34/35 * params.Re * F1V * sv(domain.diff(F1, [1, 0]')) - hV.^2 * sv(domain.diff(P, [1, 0]')) - 1/3 * hV.^3 * domain.diffMat([1, 0]') * dPdh;
    dR2dF1 = - speye(N) + 36/35 * params.Re * sv(domain.diff(h, [1, 0]')) * F1V - 34/35 * params.Re * hV * sv(domain.diff(F1, [1, 0]')) - 34/35 * params.Re * hV * F1V * domain.diffMat([1, 0]');
    
    dQ2dh = (5/(2 * params.Re)) * (hV^(-2) * dR2dh - 2 * hV^(-3) * R2);
    dQ2dF1 = (5/(2 * params.Re)) * hV^(-2) * dR2dF1;
    
    jacobian = [dQ1dh, dQ1dF1; dQ2dh, dQ2dF1];
    
    function sparseVector = toSparseVector(domain, in, N)
        vector = domain.reshapeToVector(in);
        sparseVector = spdiags(vector, 0, N, N);
    end
end