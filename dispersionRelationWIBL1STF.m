function omega = dispersionRelationWIBL1STF(domain, hbar, F1bar, params, alpha, modes)
    tolerance = 1e-3;

    omega = nan(length(alpha), modes);

    for n = 1:length(alpha)
        A = linearisedWIBL1STF(domain, hbar, F1bar, params, alpha(n));

        [V, D] = eig(A);
        val = diag(D);

        volume = abs(real(sum(V(1:end/2,:))));
        flux1 = abs(real(sum(V(1+end/2:end,:))));

        zeroVolumeIndices = volume > tolerance;
        zeroFlux1Indices = flux1 > tolerance;

        %val = val(logical(zeroVolumeIndices));

        [~,I] = sort(real(val), "descend");

        m = min(length(I),modes);

        omega(n,1:m) = val(I(1:m));
    end
end
