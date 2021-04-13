function [A] = linearisedWIBL1STF(domain, hbar, F1bar, params, alpha)
    beta = 0;
    epsilon = 1;

    N = prod(domain.shape)*2;

    hbarz = domain.ifft(domain.diff(domain.fft(hbar), 1));
    hbarz2 = domain.ifft(domain.diff(domain.fft(hbar), 2));
    hbarz3 = domain.ifft(domain.diff(domain.fft(hbar), 3));
    hbarz4 = domain.ifft(domain.diff(domain.fft(hbar), 4));

    A = zeros(N, N);

    for n = 1:N
        y = zeros(N,1);
        y(n) = 1;
        
        fhhat = y(1:end/2);
        fF1hat = y(1+end/2:end);

        hhat = domain.ifft(fhhat);
        hhatz = domain.ifft(domain.diff(fhhat, 1));
        hhatz2 = domain.ifft(domain.diff(fhhat, 2));
        hhatz3 = domain.ifft(domain.diff(fhhat, 3));
        hhatz4 = domain.ifft(domain.diff(fhhat, 4));

        F1hat = domain.ifft(fF1hat);

        A1 = -1i.*alpha.*F1hat - 2.*beta.^2.*hbar.^3.*hhat.*cot(params.theta)./3 + 4.*1i.*beta.*hbar.^3.*cot(params.theta).*hhatz./3 + 4.*1i.*beta.*hbar.^2.*hhat.*cot(params.theta).*hbarz + 2.*hbar.^3.*cot(params.theta).*hhatz2./3 + 2.*hbar.^2.*hhat.*cot(params.theta).*hbarz2 + 4.*hbar.^2.*cot(params.theta).*hbarz.*hhatz + 4.*hbar.*hhat.*cot(params.theta).*hbarz.^2 - alpha.*beta.^3.*hbar.^3.*hhat./(3.*params.C.*epsilon) + 1i.*alpha.*beta.^2.*hbar.^3.*hhatz./(params.C.*epsilon) + 1i.*alpha.*beta.^2.*hbar.^2.*hhat.*hbarz./(params.C.*epsilon) + alpha.*beta.*hbar.^3.*hhatz2./(params.C.*epsilon) + 2.*alpha.*beta.*hbar.^2.*hbarz.*hhatz./(params.C.*epsilon) - 1i.*alpha.*hbar.^3.*hhatz3./(3.*params.C.*epsilon) - 1i.*alpha.*hbar.^2.*hbarz.*hhatz2./(params.C.*epsilon) - beta.^4.*hbar.^3.*hhat./(3.*params.C.*epsilon) + 4.*1i.*beta.^3.*hbar.^3.*hhatz./(3.*params.C.*epsilon) + 1i.*beta.^3.*hbar.^2.*hhat.*hbarz./(params.C.*epsilon) + 2.*beta.^2.*hbar.^3.*hhatz2./(params.C.*epsilon) + 3.*beta.^2.*hbar.^2.*hbarz.*hhatz./(params.C.*epsilon) - 4.*1i.*beta.*hbar.^3.*hhatz3./(3.*params.C.*epsilon) - 1i.*beta.*hbar.^2.*hhat.*hbarz3./(params.C.*epsilon) - 3.*1i.*beta.*hbar.^2.*hbarz.*hhatz2./(params.C.*epsilon) - hbar.^3.*hhatz4./(3.*params.C.*epsilon) - hbar.^2.*hhat.*hbarz4./(params.C.*epsilon) - hbar.^2.*hbarz.*hhatz3./(params.C.*epsilon) - hbar.^2.*hbarz3.*hhatz./(params.C.*epsilon) - 2.*hbar.*hhat.*hbarz.*hbarz3./(params.C.*epsilon);

        A2 = 9.*1i.*alpha.*F1bar.^2.*hhat./(7.*hbar.^2) - 17.*1i.*alpha.*F1bar.*F1hat./(7.*hbar) - 5.*1i.*alpha.*hbar.*hhat.*cot(params.theta)./(3.*params.Re) + 5.*F1bar.*hhat./(params.Re.*epsilon.*hbar.^3) - 5.*F1hat./(2.*params.Re.*epsilon.*hbar.^2) + 5.*hhat./(3.*params.Re.*epsilon) - 5.*1i.*alpha.^3.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*1i.*alpha.^2.*beta.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*alpha.^2.*hbar.*hhatz./(6.*params.C.*params.Re.*epsilon);

        A(1:N,n) = [domain.fft(A1); ...
            domain.fft(A2)];
    end

end
