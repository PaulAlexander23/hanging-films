%MAIN

addpath("../time-stepping-methods/", ...
    "../optimisation-methods/", ...
    "../discretisation-methods/")

%% SETUP x
dim = 2;
xL = [32, 32];
xN = [64, 64];
xS = xL ./ xN;
x = cell(dim, 1);
for n = 1:dim
    x{n} = linspace(xS(n), xL(n), xN(n))';
end

%% SETUP t
tL = 10;
tN = 100;
tS = tL / tN;
t = 0:tS:tL;

%% SETUP y0

interface = @(x) iload(x, 'test');
y0 = irivulet(x);

plot_surface(x, y0')

%% SETUP problem
params = [1, 7 * pi / 8, 1, 0.01]; % delta, theta, Re, C

% SETUP differentiation method
problemDeg = [1, 0; 0, 1; 2, 0; 0, 2]';
D = init_fd(x, problemDeg, 2);
getD = @(deg) get_fd(deg, D, problemDeg);
method = @(x, y, deg) diff_fd(x, y, deg, D, problemDeg);
% method = @diff_ps_2d;

% pdefun = @(t,x,y,method) fbenney(x,y,params,method);
pdefun = @(t, x, y, method) fwibl1(x, y, params, method);
y0 = cat(1, y0, 0*y0);

%% SETUP timestepping method

% options = struct( ...
%     'Jacobian', @(t, y) jbenney(x, y, params, method, getD) ...
%     );
% timestepper = @(odefun,t,y0) bdf(odefun,t,y0,options);

odeopt = odeset( ...
    ... 'Jacobian', @(t, y) jbenney(x, y, params, method, getD), ...
    'Event', @event_dewetted, ...
    'Vectorized', 'on' ...
    ... 'BDF','on',... % Backward differentiation formulae
    );
timestepper = @(odefun, t, y0) ode15s(odefun, t, y0, odeopt);
%timestepper = @(odefun,t,y0) ode45(odefun,t,y0,odeopt);

%% Solve
tic
[y, t] = solver(pdefun, t, x, y0, method, timestepper);
timeTaken = toc

%% Crop solution if method is wibl1
if func2str(pdefun) == "@(t,x,y,method)fwibl1(x,y,params,method)"
    F = y(end/2+1:end, :, :);
    y = y(1:end/2, :, :);
end
