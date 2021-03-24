function A = linearisedWIBL2(domain, hbar, F1bar, F2bar, params, alpha)
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
    F1hatz2 = [Z, Dzz, Z];

    F2hat = [Z, Z, I];
    F2hatz = [Z, Z, Dz];
    F2hatz2 = [Z, Z, Dzz];

    A1 = -1i.*alpha.*F1hat - 1i.*beta.*F2hat - domain.diff(F2hat, 1);

    A2 = -17.*1i.*alpha.*F1bar.*F1hat./(7.*hbar) + 9.*1i.*alpha.*hhat.*F1bar.^2./(7.*hbar.^2) - 8.*1i.*beta.*F1bar.*F2hat./(7.*hbar) - 9.*1i.*beta.*F1hat.*F2bar./(7.*hbar) + 9.*1i.*beta.*hhat.*F1bar.*F2bar./(7.*hbar.^2) - 8.*F1bar.*domain.diff(F2hat, 1)./(7.*hbar) - 8.*F1hat.*domain.diff(F2bar, 1)./(7.*hbar) - 9.*F2bar.*domain.diff(F1hat, 1)./(7.*hbar) - 9.*F2hat.*domain.diff(F1bar, 1)./(7.*hbar) + 8.*hhat.*F1bar.*domain.diff(F2bar, 1)./(7.*hbar.^2) + 9.*hhat.*F2bar.*domain.diff(F1bar, 1)./(7.*hbar.^2) + 9.*F1bar.*F2bar.*hhatz./(7.*hbar.^2) + 9.*F1bar.*F2hat.*domain.diff(hbar, 1)./(7.*hbar.^2) + 9.*F1hat.*F2bar.*domain.diff(hbar, 1)./(7.*hbar.^2) - 18.*hhat.*F1bar.*F2bar.*domain.diff(hbar, 1)./(7.*hbar.^3) - 9.*alpha.^2.*epsilon.*F1hat./(2.*params.Re) + 6.*alpha.^2.*epsilon.*hhat.*F1bar./(params.Re.*hbar) - 7.*alpha.*beta.*epsilon.*F2hat./(2.*params.Re) + 73.*alpha.*beta.*epsilon.*hhat.*F2bar./(16.*params.Re.*hbar) + 7.*1i.*alpha.*epsilon.*domain.diff(F2hat, 1)./(2.*params.Re) - 13.*1i.*alpha.*epsilon.*hhat.*domain.diff(F2bar, 1)./(16.*params.Re.*hbar) - 73.*1i.*alpha.*epsilon.*F2bar.*hhatz./(16.*params.Re.*hbar) - 43.*1i.*alpha.*epsilon.*F2hat.*domain.diff(hbar, 1)./(16.*params.Re.*hbar) + 13.*1i.*alpha.*epsilon.*hhat.*F2bar.*domain.diff(hbar, 1)./(4.*params.Re.*hbar.^2) - 5.*1i.*alpha.*hbar.*hhat.*cot(params.theta)./(3.*params.Re) - beta.^2.*epsilon.*F1hat./params.Re + 23.*beta.^2.*epsilon.*hhat.*F1bar./(16.*params.Re.*hbar) + 2.*1i.*beta.*epsilon.*domain.diff(F1hat, 1)./params.Re - 1i.*beta.*epsilon.*hhat.*domain.diff(F1bar, 1)./(params.Re.*hbar) - 23.*1i.*beta.*epsilon.*F1bar.*hhatz./(8.*params.Re.*hbar) - 1i.*beta.*epsilon.*F1hat.*domain.diff(hbar, 1)./(params.Re.*hbar) + 3.*1i.*beta.*epsilon.*hhat.*F1bar.*domain.diff(hbar, 1)./(2.*params.Re.*hbar.^2) + epsilon.*domain.diff(F1hat, 2)./params.Re - 23.*epsilon.*F1bar.*hhatz2./(16.*params.Re.*hbar) - 23.*epsilon.*F1hat.*domain.diff(hbar, 2)./(16.*params.Re.*hbar) - epsilon.*domain.diff(hbar, 1).*domain.diff(F1hat, 1)./(params.Re.*hbar) - epsilon.*hhatz.*domain.diff(F1bar, 1)./(params.Re.*hbar) + 23.*epsilon.*hhat.*F1bar.*domain.diff(hbar, 2)./(16.*params.Re.*hbar.^2) + epsilon.*hhat.*domain.diff(hbar, 1).*domain.diff(F1bar, 1)./(params.Re.*hbar.^2) + 3.*epsilon.*F1bar.*domain.diff(hbar, 1).*hhatz./(2.*params.Re.*hbar.^2) + 3.*epsilon.*F1hat.*domain.diff(hbar, 1).^2./(4.*params.Re.*hbar.^2) - 3.*epsilon.*hhat.*F1bar.*domain.diff(hbar, 1).^2./(2.*params.Re.*hbar.^3) + 5.*hhat./(3.*params.Re.*epsilon) - 5.*F1hat./(2.*params.Re.*epsilon.*hbar.^2) + 5.*hhat.*F1bar./(params.Re.*epsilon.*hbar.^3) - 5.*1i.*alpha.^3.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*1i.*alpha.*beta.^2.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*alpha.*beta.*hbar.*hhatz./(3.*params.C.*params.Re.*epsilon) + 5.*1i.*alpha.*hbar.*hhatz2./(6.*params.C.*params.Re.*epsilon);

    A3 = -9.*1i.*alpha.*F1bar.*F2hat./(7.*hbar) - 8.*1i.*alpha.*F1hat.*F2bar./(7.*hbar) + 9.*1i.*alpha.*hhat.*F1bar.*F2bar./(7.*hbar.^2) - 17.*1i.*beta.*F2bar.*F2hat./(7.*hbar) + 9.*1i.*beta.*hhat.*F2bar.^2./(7.*hbar.^2) - 17.*F2bar.*domain.diff(F2hat, 1)./(7.*hbar) - 17.*F2hat.*domain.diff(F2bar, 1)./(7.*hbar) + 17.*hhat.*F2bar.*domain.diff(F2bar, 1)./(7.*hbar.^2) + 9.*F2bar.^2.*hhatz./(7.*hbar.^2) + 18.*F2bar.*F2hat.*domain.diff(hbar, 1)./(7.*hbar.^2) - 18.*hhat.*F2bar.^2.*domain.diff(hbar, 1)./(7.*hbar.^3) - alpha.^2.*epsilon.*F2hat./params.Re + 23.*alpha.^2.*epsilon.*hhat.*F2bar./(16.*params.Re.*hbar) - 7.*alpha.*beta.*epsilon.*F1hat./(2.*params.Re) + 73.*alpha.*beta.*epsilon.*hhat.*F1bar./(16.*params.Re.*hbar) + 7.*1i.*alpha.*epsilon.*domain.diff(F1hat, 1)./(2.*params.Re) - 43.*1i.*alpha.*epsilon.*hhat.*domain.diff(F1bar, 1)./(16.*params.Re.*hbar) - 73.*1i.*alpha.*epsilon.*F1bar.*hhatz./(16.*params.Re.*hbar) - 13.*1i.*alpha.*epsilon.*F1hat.*domain.diff(hbar, 1)./(16.*params.Re.*hbar) + 13.*1i.*alpha.*epsilon.*hhat.*F1bar.*domain.diff(hbar, 1)./(4.*params.Re.*hbar.^2) - 9.*beta.^2.*epsilon.*F2hat./(2.*params.Re) + 6.*beta.^2.*epsilon.*hhat.*F2bar./(params.Re.*hbar) + 9.*1i.*beta.*epsilon.*domain.diff(F2hat, 1)./params.Re - 9.*1i.*beta.*epsilon.*hhat.*domain.diff(F2bar, 1)./(2.*params.Re.*hbar) - 12.*1i.*beta.*epsilon.*F2bar.*hhatz./(params.Re.*hbar) - 9.*1i.*beta.*epsilon.*F2hat.*domain.diff(hbar, 1)./(2.*params.Re.*hbar) + 8.*1i.*beta.*epsilon.*hhat.*F2bar.*domain.diff(hbar, 1)./(params.Re.*hbar.^2) - 5.*1i.*beta.*hbar.*hhat.*cot(params.theta)./(3.*params.Re) + 9.*epsilon.*domain.diff(F2hat, 2)./(2.*params.Re) - 6.*epsilon.*F2bar.*hhatz2./(params.Re.*hbar) - 6.*epsilon.*F2hat.*domain.diff(hbar, 2)./(params.Re.*hbar) - 9.*epsilon.*domain.diff(hbar, 1).*domain.diff(F2hat, 1)./(2.*params.Re.*hbar) - 9.*epsilon.*hhatz.*domain.diff(F2bar, 1)./(2.*params.Re.*hbar) + 6.*epsilon.*hhat.*F2bar.*domain.diff(hbar, 2)./(params.Re.*hbar.^2) + 9.*epsilon.*hhat.*domain.diff(hbar, 1).*domain.diff(F2bar, 1)./(2.*params.Re.*hbar.^2) + 8.*epsilon.*F2bar.*domain.diff(hbar, 1).*hhatz./(params.Re.*hbar.^2) + 4.*epsilon.*F2hat.*domain.diff(hbar, 1).^2./(params.Re.*hbar.^2) - 8.*epsilon.*hhat.*F2bar.*domain.diff(hbar, 1).^2./(params.Re.*hbar.^3) - 5.*hbar.*cot(params.theta).*hhatz./(3.*params.Re) - 5.*hhat.*cot(params.theta).*domain.diff(hbar, 1)./(3.*params.Re) - 5.*F2hat./(2.*params.Re.*epsilon.*hbar.^2) + 5.*hhat.*F2bar./(params.Re.*epsilon.*hbar.^3) - 5.*1i.*alpha.^2.*beta.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*alpha.^2.*hbar.*hhatz./(6.*params.C.*params.Re.*epsilon) - 5.*1i.*beta.^3.*hbar.*hhat./(6.*params.C.*params.Re.*epsilon) - 5.*beta.^2.*hbar.*hhatz./(2.*params.C.*params.Re.*epsilon) + 5.*1i.*beta.*hbar.*hhatz2./(2.*params.C.*params.Re.*epsilon) + 5.*hbar.*hhatz3./(6.*params.C.*params.Re.*epsilon) + 5.*hhat.*domain.diff(hbar, 3)./(6.*params.C.*params.Re.*epsilon);

    A = [A1; A2; A3];

end
