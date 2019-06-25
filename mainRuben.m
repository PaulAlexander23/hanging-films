addpath discretisationMethods


x = setupX(32,32,64,64);
params = [1, 7 * pi / 8, 1, 0.01]; % delta, theta, Re, C
domain = FDDomain(x, [1, 0; 0, 1; 2, 0; 0, 2]', 4); % Finite Difference Domain
% domain = PSDomain(x); % Finite Difference Domain
pdefun = @(t, domain, y) fwibl1(domain, y, params);

% h = load('wibl1DataForRuben.mat','h'); % Note: h is y and F combine

y = h(1:end/2,:);
F = h(end/2+1:end,:);
% F = zeros(size(F)); % set F1 to zero instead of using numerics.

dhdt = fwibl1(domain, h, params);
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
