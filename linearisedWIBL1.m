function [A] = linearisedWIBL1(domain, hbar, F1bar, F2bar, params, alpha)
    beta = 0;
    epsilon = 1;

    I = eye(prod(domain.shape));
    Z = zeros(prod(domain.shape));

    Dz = domain.diffMat(1);
    Dzz = domain.diffMat(2);
    Dz3 = domain.diffMat(3);

    hhat = [I, Z, Z];
    hhatz = [Dz, Z, Z];
    hhatz2 = [Dzz, Z, Z];
    hhatz3 = [Dz3, Z, Z];

    F1hat = [Z, I, Z];
    F1hatz = [Z, Dz, Z];

    F2hat = [Z, Z, I];
    F2hatz = [Z, Z, Dz];

    A1 = -1i.*alpha.*F1hat - 1i.*beta.*F2hat - F2hatz;

    A2 = 9.*1i.*alpha.*F1bar.^2.*hhat./(7.*hbar.^2) - 17.*1i.*alpha.*F1bar.*F1hat./(7.*hbar) + 9.*1i.*beta.*F1bar.*F2bar.*hhat./(7.*hbar.^2) - 8.*1i.*beta.*F1bar.*F2hat./(7.*hbar) - 9.*1i.*beta.*F1hat.*F2bar./(7.*hbar) + 9.*F1bar.*F2bar.*hhatz./(7.*hbar.^2) - 18.*F1bar.*F2bar.*hhat.*domain.diff(hbar, 1)./(7.*hbar.^3) + 9.*F1bar.*F2hat.*domain.diff(hbar, 1)./(7.*hbar.^2) - 8.*F1bar.*F2hatz./(7.*hbar) + 8.*F1bar.*hhat.*domain.diff(F2bar, 1)./(7.*hbar.^2) + 9.*F1hat.*F2bar.*domain.diff(hbar, 1)./(7.*hbar.^2) - 8.*F1hat.*domain.diff(F2bar, 1)./(7.*hbar) - 9.*F2bar.*F1hatz./(7.*hbar) + 9.*F2bar.*hhat.*domain.diff(F1bar, 1)./(7.*hbar.^2) - 9.*F2hat.*domain.diff(F1bar, 1)./(7.*hbar) - 5.*1i.*alpha.*hbar.*hhat.*cot(params.theta)./(3.*params.Re) + 5.*F1bar.*hhat./(params.Re.*epsilon.*hbar.^3) - 5.*F1hat./(2.*params.Re.*epsilon.*hbar.^2) + 5.*hhat./(3.*params.Re.*epsilon) - 5.*1i.*alpha.^3.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*1i.*alpha.^2.*beta.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*alpha.^2.*hbar.*hhatz./(6.*params.C.*params.Re.*epsilon);

    A3 = 9.*1i.*alpha.*F1bar.*F2bar.*hhat./(7.*hbar.^2) - 9.*1i.*alpha.*F1bar.*F2hat./(7.*hbar) - 8.*1i.*alpha.*F1hat.*F2bar./(7.*hbar) + 9.*1i.*beta.*F2bar.^2.*hhat./(7.*hbar.^2) - 17.*1i.*beta.*F2bar.*F2hat./(7.*hbar) + 9.*F2bar.^2.*hhatz./(7.*hbar.^2) - 18.*F2bar.^2.*hhat.*domain.diff(hbar, 1)./(7.*hbar.^3) + 18.*F2bar.*F2hat.*domain.diff(hbar, 1)./(7.*hbar.^2) - 17.*F2bar.*F2hatz./(7.*hbar) + 17.*F2bar.*hhat.*domain.diff(F2bar, 1)./(7.*hbar.^2) - 17.*F2hat.*domain.diff(F2bar, 1)./(7.*hbar) - 5.*1i.*beta.*hbar.*hhat.*cot(params.theta)./(3.*params.Re) - 5.*hbar.*cot(params.theta).*hhatz./(3.*params.Re) - 5.*hhat.*cot(params.theta).*domain.diff(hbar, 1)./(3.*params.Re) + 5.*F2bar.*hhat./(params.Re.*epsilon.*hbar.^3) - 5.*F2hat./(2.*params.Re.*epsilon.*hbar.^2) - 5.*1i.*alpha.*beta.^2.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*alpha.*beta.*hbar.*hhatz./(3.*params.C.*params.Re.*epsilon) + 5.*1i.*alpha.*hbar.*hhatz2./(6.*params.C.*params.Re.*epsilon) - 5.*1i.*beta.^3.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*beta.^2.*hbar.*hhatz./(2.*params.C.*params.Re.*epsilon) + 5.*1i.*beta.*hbar.*hhatz2./(2.*params.C.*params.Re.*epsilon) + 5.*hbar.*hhatz3./(6.*params.C.*params.Re.*epsilon) + 5.*hhat.*domain.diff(hbar, 3)./(6.*params.C.*params.Re.*epsilon);

    A = [A1; A2; A3];

end
