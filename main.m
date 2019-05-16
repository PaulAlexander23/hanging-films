%MAIN

addpath("../time-stepping-methods/",...
    "../optimisation-methods/",...
    "../discretisation-methods/")

%% SETUP x
dim = 2;
xL = [32,32];
xN = [64,64];
xS = xL./xN;
x = cell(dim,1);
for n = 1:dim
    x{n} = linspace(xS(n),xL(n),xN(n))';
end

%% SETUP t
tL = 10;
tN = 100;
tS = tL/tN;
t = 0:tS:tL;

%% SETUP y0

interface = @(x) iload(x,'test');
y0 = irivulet(x);

plot_surface(x,y0')

%% SETUP problem
params = [1,7*pi/8,1,0.01]; % delta, theta, Re, C

% SETUP differentiation method
problemDeg = [1,0;0,1;2,0;0,2]';
D = init_fd(x, problemDeg, 2);
getD = @(deg) get_fd(deg,D,problemDeg);
method = @(x,y,deg) diff_fd(x,y,deg,D,problemDeg);
% method = @diff_ps_2d;

% pdefun = @(t,x,y,method) fbenney(x,y,params,method);
pdefun = @(t,x,y,method) fwibl1(x,y,params,method);
y0 = cat(1, y0, 0*y0);

%% SETUP timestepping method

% options = struct( ...
%     'Jacobian', @(t, y) jbenney(x, y, params, method, getD) ...
%     );
% timestepper = @(odefun,t,y0) bdf(odefun,t,y0,options);

odeopt = odeset( ...
    ...'Jacobian', @(t, y) jbenney(x, y, params, method, getD), ...
    'Event', @event_dewetted ...
    ...'Vectorized','on',...
    ...'BDF','on',... % Backward differentiation formulae
    );
timestepper = @(odefun,t,y0) ode15s(odefun,t,y0,odeopt);
%timestepper = @(odefun,t,y0) ode45(odefun,t,y0,odeopt);

%% Solve
tic
[y, t] = solver(pdefun, t, x, y0, method, timestepper);
timeTaken = toc

%% Crop solution if method is wibl1
if func2str(pdefun) == '@(t,x,y,method)fwibl1(x,y,params,method)'
    F = y(end/2+1:end,:,:);
    y = y(1:end/2,:,:);
end
%% Plot overview
figure
plot_surface(x,real(y(:,:,end))');

figure
plot_log_fourier(x,real(y(:,:,end-1))');

figure
plot(t,log10(squeeze(energy(x,y))));

figure
plot(t,squeeze(min(y,[],[1,2])));

% Minimum location
[temp,tempI] = min(y,[],1);
[miny,yI] = min(temp,[],2);
xI = zeros(size(yI));
for m = 1:length(tempI)
    xI(1,1,m) = tempI(1,yI(m),m);
end
figure
plot(t,squeeze(x{1}(xI)));
figure
plot(t,squeeze(x{2}(yI)));

%% H2 norm - Check implemented correctly
ycell = mat2cell(reshape(y,numel(y0),length(t))',ones(length(t),1));
figure
plot(t,log10(cellfun(@(y) norm_H2(x,y',method),ycell)));

%% Construct other variables
% y_snapshot = y0;
y_snapshot = y(:,:,end);

p = construct_pressure(x,y_snapshot,params,method);

[z,u,v,w] = construct_velocity(x,y_snapshot,params,method);

slice = 32;

figure
surf(repmat(x{1},1,size(z(:,slice,:),3)),squeeze(u(:,slice,:)),squeeze(z(:,slice,:)));
shading interp;

figure
surf(repmat(x{1},1,size(z(:,slice,:),3)),squeeze(v(:,slice,:)),squeeze(z(:,slice,:)));
shading interp;

figure
surf(repmat(x{1},1,size(z(:,slice,:),3)),squeeze(w(:,slice,:)),squeeze(z(:,slice,:)));
shading interp;

figure
xstep = 1;
zstep = 10;
quiver(repmat(x{1}(1:xstep:end),1,floor(size(z,3)/zstep))',...
    squeeze(z(1:xstep:end,slice,1:zstep:end))',...
    squeeze(w(1:xstep:end,slice,1:zstep:end))',...
    squeeze(v(1:xstep:end,slice,1:zstep:end))',...
    1); % Scaling

%% Animate

addpath('~/Repositories/plotting/')
figure
% F = animate(@(x,y) plot_surface(x,y'),t,x,y);
F = animateSmooth(@(x,y) plot_surface(x,y'),t,x,y,t(end)/60);

% movie(F);