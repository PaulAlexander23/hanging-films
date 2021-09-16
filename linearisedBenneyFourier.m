function A = linearisedBenneyFourier(domain, hbar, params, alpha)
    beta = 0;
    epsilon = 1;

    N = prod(domain.shape);

    hbarz = domain.ifft(domain.diff(domain.fft(hbar), 1));
    hbarz2 = domain.ifft(domain.diff(domain.fft(hbar), 2));
    hbarz3 = domain.ifft(domain.diff(domain.fft(hbar), 3));
    hbarz4 = domain.ifft(domain.diff(domain.fft(hbar), 4));

    A = zeros(N,N);

    for n = 1:N
        fhhat = zeros(N,1);
        fhhat(n) = 1;

        hhat = domain.ifft(fhhat);
        hhatz = domain.ifft(domain.diff(fhhat, 1));
        hhatz2 = domain.ifft(domain.diff(fhhat, 2));
        hhatz3 = domain.ifft(domain.diff(fhhat, 3));
        hhatz4 = domain.ifft(domain.diff(fhhat, 4));

        row = 8.*params.Re.*alpha.^2.*epsilon.*hbar.^6.*hhat./15 ...
            - 2.*alpha.^2.*epsilon.*hbar.^3.*hhat.*cot(params.theta)./3 ...
            - 2.*1i.*alpha.*hbar.^2.*hhat ...
            - 2.*beta.^2.*epsilon.*hbar.^3.*hhat.*cot(params.theta)./3 ...
            + 4.*1i.*beta.*epsilon.*hbar.^3.*cot(params.theta).*hhatz./3 ...
            + 4.*1i.*beta.*epsilon.*hbar.^2.*hhat.*cot(params.theta).*hbarz ...
            + 2.*epsilon.*hbar.^3.*cot(params.theta).*hhatz2./3 ...
            + 2.*epsilon.*hbar.^2.*hhat.*cot(params.theta).*hbarz2 ...
            + 4.*epsilon.*hbar.^2.*cot(params.theta).*hbarz.*hhatz ...
            + 4.*epsilon.*hbar.*hhat.*cot(params.theta).*hbarz.^2 ...
            - alpha.^4.*epsilon.^3.*hbar.^3.*hhat./(3.*params.C) ...
            - 2.*alpha.^2.*beta.^2.*epsilon.^3.*hbar.^3.*hhat./(3.*params.C) ...
            + 4.*1i.*alpha.^2.*beta.*epsilon.^3.*hbar.^3.*hhatz./(3.*params.C) ...
            + 1i.*alpha.^2.*beta.*epsilon.^3.*hbar.^2.*hhat.*hbarz./params.C ...
            + 2.*alpha.^2.*epsilon.^3.*hbar.^3.*hhatz2./(3.*params.C) ...
            + alpha.^2.*epsilon.^3.*hbar.^2.*hbarz.*hhatz./params.C ...
            - beta.^4.*epsilon.^3.*hbar.^3.*hhat./(3.*params.C) ...
            + 4.*1i.*beta.^3.*epsilon.^3.*hbar.^3.*hhatz./(3.*params.C) ...
            + 1i.*beta.^3.*epsilon.^3.*hbar.^2.*hhat.*hbarz./params.C ...
            + 2.*beta.^2.*epsilon.^3.*hbar.^3.*hhatz2./params.C ...
            + 3.*beta.^2.*epsilon.^3.*hbar.^2.*hbarz.*hhatz./params.C ...
            - 4.*1i.*beta.*epsilon.^3.*hbar.^3.*hhatz3./(3.*params.C) ...
            - 1i.*beta.*epsilon.^3.*hbar.^2.*hhat.*hbarz3./params.C ...
            - 3.*1i.*beta.*epsilon.^3.*hbar.^2.*hbarz.*hhatz2./params.C ...
            - epsilon.^3.*hbar.^3.*hhatz4./(3.*params.C) ...
            - epsilon.^3.*hbar.^2.*hhat.*hbarz4./params.C ...
            - epsilon.^3.*hbar.^2.*hbarz.*hhatz3./params.C ...
            - epsilon.^3.*hbar.^2.*hbarz3.*hhatz./params.C ...
            - 2.*epsilon.^3.*hbar.*hhat.*hbarz.*hbarz3./params.C;
    
        A(:,n) = domain.fft(row);
    end

end