function eigF = eigenfunctionWIBL2(domain, hbar, F1bar, F2bar, params, alpha)
    tolerance = 1e-3;

    A = linearisedWIBL2(domain, hbar, F1bar, F2bar, params, alpha);

    [V, D] = eig(A);
    val = diag(D);

    volume = abs(real(sum(V(1:end/3,:))));
    flux1 = abs(real(sum(V(1+end/3:2*end/3,:))));
    flux2 = abs(real(sum(V(1+2*end/3:end,:))));

    zeroVolumeIndices = volume > tolerance;
    zeroFlux1Indices = flux1 > tolerance;
    zeroFlux2Indices = flux2 > tolerance;

    val = val(logical(zeroVolumeIndices));
    V = V(:,logical(zeroVolumeIndices));

    [~,I] = sort(real(val), "descend");

    eigF = V(:,I(1));
end
