% addpath discretisationMethods
addpath ../discretisation-methods/

x = setupX(32,32,64,64);
params = [1, 7 * pi / 8, 1, 0.01]; % delta, theta, Re, C
domain = FDDomain(x, [1, 0; 0, 1; 2, 0; 0, 2]', 4);
% domain = PSDomain(x);
pdefun = @(t, domain, y) fwibl1(domain, y, params);

vars = load('wibl1DataForRuben.mat','h'); % Note: h is y and F combine
h = vars.h;

y = (h(1:end/2,:));
F = (h(end/2+1:end,:));
% F = 2/3 * ones(size(F)); % set F1 to zero instead of using numerics.

%% Smoothing
fy = fftn(y);
fy(abs(fy)<1e0) = 0;
y = real(ifftn(fy));

fF = fftn(F);
fF(abs(fF)<1e-3) = 0;
F = real(ifftn(fF));

%%
dhdt = fwibl1(domain, [y; F], params);
dydt = dhdt(1:end/2,:);
dFdt = dhdt(end/2+1:end,:);

%%
dydx = domain.diff(y, [1, 0]');

c = mean(dydt./-dydx,'all');

disp(c)
disp(norm(dydt + c * dydx))

figure
surf(dydt) % You can see some numerical error here.

figure
surf(- c * dydx)

%%
dFdx = domain.diff(F, [1, 0]');

c = mean(dFdt./-dFdx,'all');

disp(c)
disp(norm(dFdt + c * dFdx))

figure
surf(dFdt) % You can see some numerical error here.

figure
surf(- c * dFdx)

%%

dhdx = [dydx;dFdx];

c = mean(dhdt./-dhdx,'all');

disp(c)
disp(norm(dhdt + c * dhdx))

figure
surf(dhdt) % You can see some numerical error here.

figure
surf(- c * dhdx)
