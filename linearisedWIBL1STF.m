function [A] = linearisedWIBL1STF(domain, hbar, F1bar, params, alpha)
    beta = 0;
    epsilon = 1;

    I = eye(prod(domain.shape));
    Z = zeros(prod(domain.shape));

    Dz = domain.diffMat(1);
    Dzz = domain.diffMat(2);
    Dz3 = domain.diffMat(3);
    Dz4 = domain.diffMat(4);

    hhat = [I, Z];
    hhatz = [Dz, Z];
    hhatz2 = [Dzz, Z];
    hhatz3 = [Dz3, Z];
    hhatz4 = [Dz4, Z];

    F1hat = [Z, I];

    A1 = -1i.*alpha.*F1hat - 2.*beta.^2.*hbar.^3.*hhat.*cot(params.theta)./3 + 4.*1i.*beta.*hbar.^3.*cot(params.theta).*hhatz./3 + 4.*1i.*beta.*hbar.^2.*hhat.*cot(params.theta).*domain.diff(hbar, 1) + 2.*hbar.^3.*cot(params.theta).*hhatz2./3 + 2.*hbar.^2.*hhat.*cot(params.theta).*domain.diff(hbar, 2) + 4.*hbar.^2.*cot(params.theta).*domain.diff(hbar, 1).*hhatz + 4.*hbar.*hhat.*cot(params.theta).*domain.diff(hbar, 1).^2 - alpha.*beta.^3.*hbar.^3.*hhat./(3.*params.C.*epsilon) + 1i.*alpha.*beta.^2.*hbar.^3.*hhatz./(params.C.*epsilon) + 1i.*alpha.*beta.^2.*hbar.^2.*hhat.*domain.diff(hbar, 1)./(params.C.*epsilon) + alpha.*beta.*hbar.^3.*hhatz2./(params.C.*epsilon) + 2.*alpha.*beta.*hbar.^2.*domain.diff(hbar, 1).*hhatz./(params.C.*epsilon) - 1i.*alpha.*hbar.^3.*hhatz3./(3.*params.C.*epsilon) - 1i.*alpha.*hbar.^2.*domain.diff(hbar, 1).*hhatz2./(params.C.*epsilon) - beta.^4.*hbar.^3.*hhat./(3.*params.C.*epsilon) + 4.*1i.*beta.^3.*hbar.^3.*hhatz./(3.*params.C.*epsilon) + 1i.*beta.^3.*hbar.^2.*hhat.*domain.diff(hbar, 1)./(params.C.*epsilon) + 2.*beta.^2.*hbar.^3.*hhatz2./(params.C.*epsilon) + 3.*beta.^2.*hbar.^2.*domain.diff(hbar, 1).*hhatz./(params.C.*epsilon) - 4.*1i.*beta.*hbar.^3.*hhatz3./(3.*params.C.*epsilon) - 1i.*beta.*hbar.^2.*hhat.*domain.diff(hbar, 3)./(params.C.*epsilon) - 3.*1i.*beta.*hbar.^2.*domain.diff(hbar, 1).*hhatz2./(params.C.*epsilon) - hbar.^3.*hhatz4./(3.*params.C.*epsilon) - hbar.^2.*hhat.*domain.diff(hbar, 4)./(params.C.*epsilon) - hbar.^2.*domain.diff(hbar, 1).*hhatz3./(params.C.*epsilon) - hbar.^2.*domain.diff(hbar, 3).*hhatz./(params.C.*epsilon) - 2.*hbar.*hhat.*domain.diff(hbar, 1).*domain.diff(hbar, 3)./(params.C.*epsilon);

    A2 = 9.*1i.*alpha.*F1bar.^2.*hhat./(7.*hbar.^2) - 17.*1i.*alpha.*F1bar.*F1hat./(7.*hbar) - 5.*1i.*alpha.*hbar.*hhat.*cot(params.theta)./(3.*params.Re) + 5.*F1bar.*hhat./(params.Re.*epsilon.*hbar.^3) - 5.*F1hat./(2.*params.Re.*epsilon.*hbar.^2) + 5.*hhat./(3.*params.Re.*epsilon) - 5.*1i.*alpha.^3.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*1i.*alpha.^2.*beta.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*alpha.^2.*hbar.*hhatz./(6.*params.C.*params.Re.*epsilon);

    A = [A1; A2];

end
