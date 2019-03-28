%MAIN

addpath("../time-stepping-methods/",...
    "../optimisation-methods/",...
    "../discretisation-methods/")

%% SETUP x
dim = 2;
xL = [32,32];
xN = [24,40];
xS = xL./xN;
x = cell(dim,1);
for n = 1:dim
    x{n} = linspace(xS(n),xL(n),xN(n))';
end

%% SETUP t

tL = 1;
tN = 10;
tS = tL/tN;
t = 0:tS:tL;

%% SETUP y0
% A = 1e-1;
% y0 = 1 + A * (-cos(2*pi/xL(1) * x{1}) + 0*x{2}');
% y0 = 1 + A * (-cos(2*pi/xL(1) * x{1}) + sin(2*pi/xL(2) * x{2}'));

y0 = h0;

% V = zeros(1,1,25,2);
% V(1,1,:,:) = combvec(1:5,1:5)';
% R1 = rand(1,1,25);
% R2 = rand(1,1,25);
% R3 = rand(1,1,25);
% R4 = rand(1,1,25);
% y0 = 1 ...
%     + sum(1e-4 * R1 .* cos(V(1,1,:,1) .* x{1}' + V(1,1,:,2) .* x{2} + 2*pi*R2),3)...
%     + sum(1e-4 * R3 .* cos(V(1,1,:,1) .* x{1}' - V(1,1,:,2) .* x{2} + 2*pi*R4),3);

% plot_surface(x,y0')

%% SETUP problem
params = [1,7*pi/8,1.0,0,0.01]; % delta, theta, Re, We, C
problemDeg = [1,0;0,1;2,0;0,2]';

% SETUP differentiation method
D = init_fd(x, problemDeg, 2);
method = @(x,y,deg) diff_fd(x,y,deg,D,problemDeg);
getD = @(deg) get_fd(deg,D,problemDeg);

%method = @diff_ps_2d;

%problem = @fbenney;
%problem = @(x,y,params,method) fbenney2(x,y,params,method,getD);
pdefun = @(t,x,y,method) benney(x,y,params,method,getD);

%% SETUP timestepping method

% options = optimoptions('fsolve',...
%    'SpecifyObjectiveGradient',true,...
%    'Display','off');
% optimmethod = @(fun,x0) fsolve(fun,x0,options);
optimmethod = @newton;
timestepper = @(odefun,t,y0) bdf3(odefun,t,y0,optimmethod);

evnt = @(~,H) is_dewetted(H);

RelTol = 1e-3;
options = odeset(...
    ...'Vectorized','on',...
    ...'BDF','on',... % Backward differentiation formulae
    'Event',@(t,y) event(t,y),...
    ...'OutputFcn','odeprint',...
    'OutputFcn',@outputfunction,...
    'RelTol',RelTol...
    ); % Default: 1e-3

%% Solve
tic
[y, t] = solver(pdefun, t, x, y0, method, timestepper, options);
toc

%% Plot overview

plot_surface(x,real(y(:,:,end))');

figure
plot(t, squeeze(sum((y-1).^2 * x{1}(1) * x{2}(1),[1,2])))

%% Plot specific

plot_surface(x,real(y(:,:,3))');

%plot_surface(x,real(y(:,:,end)-y0)');


%% Save

filename = replace(sprintf('data-%g-%g-%g-%g-%g.mat',params),'.','_');
save(filename,'x','t','y','params');

%% Functions
function [value, isterminal, direction] = event(t,y)
    Y = reshape(y,shape);
    value = 2*evnt(t,Y) - 1;
    if any(isnan(y),'all')
        value = 1;
    end
    isterminal = 1; % Terminal
    direction = 0; % Any approach direction
end

function [status] = outputfunction(t,y,flags)
    %whos;
    status = 0;
end