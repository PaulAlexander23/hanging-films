function eigF = eigenfunctionWIBL1STF(domain, hbar, F1bar, params, alpha)
    tolerance = 1e-3;

    A = linearisedWIBL1STF(domain, hbar, F1bar, params, alpha);

    [V, D] = eig(A);
    val = diag(D);

    volume = abs(real(sum(V(1:end/2,:))));
    flux1 = abs(real(sum(V(1+end/2:end,:))));
    zeroVolumeIndices = volume > tolerance;
    zeroFlux1Indices = flux1 > tolerance;

    val = val(logical(zeroVolumeIndices));
    V = V(:,logical(zeroVolumeIndices));

    [~,I] = sort(real(val), "descend");

    eigF = V(:,I(1));
end
