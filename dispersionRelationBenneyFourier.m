function omega = dispersionRelationBenneyFourier(domain, hbar, params, alpha, modes)
    tolerance = 1e-3;

    omega = nan(length(alpha), modes);
    
    for n = 1:length(alpha)

        A = linearisedBenneyFourier(domain, hbar, params, alpha(n));

        [V, D] = eig(A);

        val = diag(D);

        volume = abs(real(sum(V)));
        zeroVolumeIndices = volume > tolerance;

        %val = val(zeroVolumeIndices);

        [~,I] = sort(real(val), "descend");

        m = min(length(I),modes);

        omega(n,1:m) = val(I(1:m));
    end
end

