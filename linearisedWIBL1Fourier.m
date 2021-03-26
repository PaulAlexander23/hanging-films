function [A] = linearisedWIBL1Fourier(domain, hbar, F1bar, F2bar, params, alpha)
    beta = 0;
    epsilon = 1;

    N = prod(domain.shape)*3;

    hbarz = domain.ifft(domain.diff(domain.fft(hbar), 1));
    hbarz3 = domain.ifft(domain.diff(domain.fft(hbar), 3));
    F1barz = domain.ifft(domain.diff(domain.fft(F1bar), 1));
    F2barz = domain.ifft(domain.diff(domain.fft(F2bar), 1));

    A = zeros(N/2, N/2);

    for n = 1:N/2
        y = zeros(N/2,1);
        y(n) = 1;
        
        fhhat = y(1:end/3);
        fF1hat = y(1+end/3:2/3*end);
        fF2hat = y(1+2*end/3:end);

        hhat = domain.ifft(fhhat);
        hhatz = domain.ifft(domain.diff(fhhat, 1));
        hhatz2 = domain.ifft(domain.diff(fhhat, 2));
        hhatz3 = domain.ifft(domain.diff(fhhat, 3));
        hhatz4 = domain.ifft(domain.diff(fhhat, 4));

        F1hat = domain.ifft(fF1hat);
        F1hatz = domain.ifft(domain.diff(fF1hat, 1));

        F2hat = domain.ifft(fF2hat);
        F2hatz = domain.ifft(domain.diff(fF2hat, 1));

        A1 = -1i.*alpha.*F1hat - 1i.*beta.*F2hat - F2hatz;

        A2 = 9.*1i.*alpha.*F1bar.^2.*hhat./(7.*hbar.^2) - 17.*1i.*alpha.*F1bar.*F1hat./(7.*hbar) + 9.*1i.*beta.*F1bar.*F2bar.*hhat./(7.*hbar.^2) - 8.*1i.*beta.*F1bar.*F2hat./(7.*hbar) - 9.*1i.*beta.*F1hat.*F2bar./(7.*hbar) + 9.*F1bar.*F2bar.*hhatz./(7.*hbar.^2) - 18.*F1bar.*F2bar.*hhat.*hbarz./(7.*hbar.^3) + 9.*F1bar.*F2hat.*hbarz./(7.*hbar.^2) - 8.*F1bar.*F2hatz./(7.*hbar) + 8.*F1bar.*hhat.*F2barz./(7.*hbar.^2) + 9.*F1hat.*F2bar.*hbarz./(7.*hbar.^2) - 8.*F1hat.*F2barz./(7.*hbar) - 9.*F2bar.*F1hatz./(7.*hbar) + 9.*F2bar.*hhat.*F1barz./(7.*hbar.^2) - 9.*F2hat.*F1barz./(7.*hbar) - 5.*1i.*alpha.*hbar.*hhat.*cot(params.theta)./(3.*params.Re) + 5.*F1bar.*hhat./(params.Re.*epsilon.*hbar.^3) - 5.*F1hat./(2.*params.Re.*epsilon.*hbar.^2) + 5.*hhat./(3.*params.Re.*epsilon) - 5.*1i.*alpha.^3.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*1i.*alpha.^2.*beta.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*alpha.^2.*hbar.*hhatz./(6.*params.C.*params.Re.*epsilon);

        A3 = 9.*1i.*alpha.*F1bar.*F2bar.*hhat./(7.*hbar.^2) - 9.*1i.*alpha.*F1bar.*F2hat./(7.*hbar) - 8.*1i.*alpha.*F1hat.*F2bar./(7.*hbar) + 9.*1i.*beta.*F2bar.^2.*hhat./(7.*hbar.^2) - 17.*1i.*beta.*F2bar.*F2hat./(7.*hbar) + 9.*F2bar.^2.*hhatz./(7.*hbar.^2) - 18.*F2bar.^2.*hhat.*hbarz./(7.*hbar.^3) + 18.*F2bar.*F2hat.*hbarz./(7.*hbar.^2) - 17.*F2bar.*F2hatz./(7.*hbar) + 17.*F2bar.*hhat.*F2barz./(7.*hbar.^2) - 17.*F2hat.*F2barz./(7.*hbar) - 5.*1i.*beta.*hbar.*hhat.*cot(params.theta)./(3.*params.Re) - 5.*hbar.*cot(params.theta).*hhatz./(3.*params.Re) - 5.*hhat.*cot(params.theta).*hbarz./(3.*params.Re) + 5.*F2bar.*hhat./(params.Re.*epsilon.*hbar.^3) - 5.*F2hat./(2.*params.Re.*epsilon.*hbar.^2) - 5.*1i.*alpha.*beta.^2.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*alpha.*beta.*hbar.*hhatz./(3.*params.C.*params.Re.*epsilon) + 5.*1i.*alpha.*hbar.*hhatz2./(6.*params.C.*params.Re.*epsilon) - 5.*1i.*beta.^3.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*beta.^2.*hbar.*hhatz./(2.*params.C.*params.Re.*epsilon) + 5.*1i.*beta.*hbar.*hhatz2./(2.*params.C.*params.Re.*epsilon) + 5.*hbar.*hhatz3./(6.*params.C.*params.Re.*epsilon) + 5.*hhat.*hbarz3./(6.*params.C.*params.Re.*epsilon);

        A(1:N/2,n) = [domain.fft(A1); ...
            domain.fft(A2); ...
            domain.fft(A3)];
    end

end
