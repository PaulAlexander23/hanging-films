function eigF = eigenfunctionBenney(domain, hbar, params, alpha)
    tolerance = 1e-3;

    A = linearisedBenney(domain, hbar, params, alpha);

    [V, D] = eig(A);

    val = diag(D);

    volume = abs(real(sum(V)));
    zeroVolumeIndices = volume > tolerance;

    val = val(zeroVolumeIndices);
    V = V(:,logical(zeroVolumeIndices));

    [~,I] = sort(real(val), "descend");

    eigF = V(:,I(1));
end
