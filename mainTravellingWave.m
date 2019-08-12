addpath discretisationMethods

% method = "finite-difference";
method = "pseudo-spectral";

x = setupX(32,32,64,64);
params = [1, 7 * pi / 8, 1, 0.01]; % delta, theta, Re, C
if method == "finite-difference"
    domain = FDDomain(x, [1, 0; 0, 1; 2, 0; 0, 2]', 4);
elseif method == "pseudo-spectral"
    domain = PSDomain(x);
end

pdefun = @(t, domain, y) fwibl1(domain, y, params);

%%
vars = load('wibl1DataForRuben.mat','h');
h = vars.h;

y = (h(1:end/2,:));
F = (h(end/2+1:end,:));
% F = 2/3 * ones(size(F)); % set F1 to zero instead of using numerics.

%% Smoothing
Nx = length(x{1})/2;
Ny = length(x{2})/2;
kx = [0:Nx-1, 0, 1-Nx:-1]';
ky = [0:Ny-1, 0, 1-Ny:-1]';

fy = fftn(y);
% figure; surf(log10(abs(real(fy))));
fy(:,abs(ky)>3/4*Ny) = 0;
fy(abs(kx)>1*Nx,:) = 0;
fy(abs(fy)<1e-5) = 0;
% figure; surf(log10(abs(real(fy))));
y = ifftn(fy, 'symmetric');

fF = fftn(F);
% figure; surf(log10(abs(real(fF))))
fF(:,abs(ky)>7/8*Ny) = 0;
fF(abs(kx)>1*Nx,:) = 0;
fF(abs(fF)<1e-3) = 0;
% figure; surf(log10(abs(real(fF))))
F = ifftn(fF, 'symmetric');

%%
if method == "pseudo-spectral"
    y = domain.fftn(y);
    F = domain.fftn(F);
end

h = [y; F];

%%
dhdt = fwibl1(domain, h, params);
dydt = dhdt(1:end/2,:);
dFdt = dhdt(end/2+1:end,:);

%%
dydx = domain.diff(y, [1, 0]');

c = apprioximateC(domain, dydt, dydx, method);
residual(domain, c, dydt, dydx, method)

plotSurface(domain, dydt, method);
plotSurface(domain, - c * dydx, method);
%%
dFdx = domain.diff(F, [1, 0]');

% c = apprioximateC(domain, dFdt, dFdx, method);
residual(domain, c, dFdt, dFdx, method)

plotSurface(domain, dFdt, method);
plotSurface(domain, - c * dFdx, method);

%%

travellingWaveResidual(domain, h, params)
% x = fsolve(@(x) , x0);


%%
function c = apprioximateC(domain, dydt, dydx, method)
    if method == "pseudo-spectral"
        dydt = domain.ifftn(dydt);
        dydx = domain.ifftn(dydx);
    end
    c = -mean(dydt./dydx,'all');
    display(c)
end

function r = residual(domain, c, dydt, dydx, method)
    if method == "pseudo-spectral"
        dydt = domain.ifftn(dydt);
        dydx = domain.ifftn(dydx);
    end
    r = norm(dydt + c * dydx);
end

function plotSurface(domain, y, method)
    figure
    if method == "pseudo-spectral"
        y = domain.ifftn(y);
    end
    surf(y)
end
