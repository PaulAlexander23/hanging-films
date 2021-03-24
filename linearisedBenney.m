function [A] = linearisedBenney(domain, hbar, params, alpha)
    beta = 0;
    epsilon = 1;

    I = eye(prod(domain.shape));

    hhat = I;
    hhatz = domain.diffMat(1);
    hhatz2 = domain.diffMat(2);
    hhatz3 = domain.diffMat(3);
    hhatz4 = domain.diffMat(4);

    A = 8.*params.Re.*alpha.^2.*epsilon.*hbar.^6.*hhat./15 - 2.*alpha.^2.*epsilon.*hbar.^3.*hhat.*cot(params.theta)./3 - 2.*1i.*alpha.*hbar.^2.*hhat - 2.*beta.^2.*epsilon.*hbar.^3.*hhat.*cot(params.theta)./3 + 4.*1i.*beta.*epsilon.*hbar.^3.*cot(params.theta).*hhatz./3 + 4.*1i.*beta.*epsilon.*hbar.^2.*hhat.*cot(params.theta).*domain.diff(hbar, 1) + 2.*epsilon.*hbar.^3.*cot(params.theta).*hhatz2./3 + 2.*epsilon.*hbar.^2.*hhat.*cot(params.theta).*domain.diff(hbar, 2) + 4.*epsilon.*hbar.^2.*cot(params.theta).*domain.diff(hbar, 1).*hhatz + 4.*epsilon.*hbar.*hhat.*cot(params.theta).*domain.diff(hbar, 1).^2 - alpha.^4.*epsilon.^3.*hbar.^3.*hhat./(3.*params.C) - 2.*alpha.^2.*beta.^2.*epsilon.^3.*hbar.^3.*hhat./(3.*params.C) + 4.*1i.*alpha.^2.*beta.*epsilon.^3.*hbar.^3.*hhatz./(3.*params.C) + 1i.*alpha.^2.*beta.*epsilon.^3.*hbar.^2.*hhat.*domain.diff(hbar, 1)./params.C + 2.*alpha.^2.*epsilon.^3.*hbar.^3.*hhatz2./(3.*params.C) + alpha.^2.*epsilon.^3.*hbar.^2.*domain.diff(hbar, 1).*hhatz./params.C - beta.^4.*epsilon.^3.*hbar.^3.*hhat./(3.*params.C) + 4.*1i.*beta.^3.*epsilon.^3.*hbar.^3.*hhatz./(3.*params.C) + 1i.*beta.^3.*epsilon.^3.*hbar.^2.*hhat.*domain.diff(hbar, 1)./params.C + 2.*beta.^2.*epsilon.^3.*hbar.^3.*hhatz2./params.C + 3.*beta.^2.*epsilon.^3.*hbar.^2.*domain.diff(hbar, 1).*hhatz./params.C - 4.*1i.*beta.*epsilon.^3.*hbar.^3.*hhatz3./(3.*params.C) - 1i.*beta.*epsilon.^3.*hbar.^2.*hhat.*domain.diff(hbar, 3)./params.C - 3.*1i.*beta.*epsilon.^3.*hbar.^2.*domain.diff(hbar, 1).*hhatz2./params.C - epsilon.^3.*hbar.^3.*hhatz4./(3.*params.C) - epsilon.^3.*hbar.^2.*hhat.*domain.diff(hbar, 4)./params.C - epsilon.^3.*hbar.^2.*domain.diff(hbar, 1).*hhatz3./params.C - epsilon.^3.*hbar.^2.*domain.diff(hbar, 3).*hhatz./params.C - 2.*epsilon.^3.*hbar.*hhat.*domain.diff(hbar, 1).*domain.diff(hbar, 3)./params.C;

end
