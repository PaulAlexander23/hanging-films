%MAIN

import discretisationMethods.*

%% SETUP x
dim = 2;
xL = [32, 32]; % [x length, y length] - Change here to 16 for other one.
xN = [64, 64];
xS = xL ./ xN;
x = cell(dim, 1);
for n = 1:dim
    x{n} = linspace(xS(n), xL(n), xN(n))';
end

%% SETUP h

% h = load('wibl1DataForRuben.mat','h'); % Note: h is y and F combine
y = h(1:end/2,:);
F = h(end/2+1:end,:);
% F = zeros(size(F)); % set F1 to zero instead of using numerics.

% surf(y)

%% SETUP problem
params = [1, 7 * pi / 8, 1, 0.01]; % delta, theta, Re, C

% SETUP differentiation method
problemDeg = [1, 0; 0, 1; 2, 0; 0, 2]';
D = init_fd(x, problemDeg, 4);
getD = @(deg) get_fd(deg, D, problemDeg);
diffMethod = @(x, y, deg) diff_fd(x, y, deg, D, problemDeg);

pdefun = @(t, x, y, method) fwibl1(x, y, params, method);

%%

dhdt = fwibl1(x, h, params, diffMethod);
dydt = dhdt(1:end/2,:);
dFdt = dhdt(end/2+1:end,:);

%%
dydx = cell2mat(diffMethod(x, y, [1, 0]'));

c = mean(dydt./-dydx,'all');

disp(c)
disp(norm(dydt + c * dydx))

figure
surf(dydt) % You can see some numerical error here.

figure
surf(- c * dydx)

%%
dFdx = cell2mat(diffMethod(x, F, [1, 0]'));

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
