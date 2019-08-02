% addpath discretisationMethods
addpath ../discretisation-methods/

x = setupX(32,32,64,64);
params = [1, 7 * pi / 8, 1, 0.01]; % delta, theta, Re, C
domain = PSDomain(x);
pdefun = @(t, domain, y) fwibl1(domain, y, params);

vars = load('wibl1DataForRuben.mat','h'); % Note: h is y and F combine
h = vars.h;

y = domain.fftn(h(1:end/2,:));
F = domain.fftn(h(end/2+1:end,:));

%% Smoothing
y(abs(y)<1e-1) = 0;
F(abs(F)<1e0) = 0;

%%
dhdt = fwibl1(domain, [y; F], params);
dydt = dhdt(1:end/2,:);
dFdt = dhdt(end/2+1:end,:);

%%
dydx = domain.diff(y, [1, 0]');

c = apprioximateC(domain, dydt, dydx)
residual(domain, c, dydt, dydx)

figure
surf(domain.ifftn(dydt)) % You can see some numerical error here.

figure
surf(- c * domain.ifftn(dydx))

%%
dFdx = domain.diff(F, [1, 0]');

c = apprioximateC(domain, dFdt, dFdx)
residual(domain, c, dFdt, dFdx)

figure
surf(domain.ifftn(dFdt)) % You can see some numerical error here.

figure
surf(- c * domain.ifftn(dFdx))

%%

dhdx = [dydx;dFdx];

c = 2.7;
residual(domain, c, dhdt, dhdx)

function c = apprioximateC(domain, dydt, dydx)
    c = mean(domain.ifftn(dydt)./-domain.ifftn(dydx),'all');
end

function r = residual(domain, c, dydt, dydx)
    r = norm(domain.ifftn(dydt) + c * domain.ifftn(dydx));
end