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
tL = 100;
tN = 1000;
tS = tL/tN;
t = 0:tS:tL;

%% SETUP y0
A = 2e-1;
r = 0.05;
% y0 = 1 + A * (-cos(2*pi/xL(1) * x{1}) + 0*x{2}');
y0 = 1 + A * (-r*cos(2*pi/xL(1) * x{1}) - cos(2*pi/xL(2) * x{2}'));

% y0 = h0;

% V = zeros(1,1,25,2);
% V(1,1,:,:) = combvec(1:5,1:5)';
% R1 = rand(1,1,25);
% R2 = rand(1,1,25);
% R3 = rand(1,1,25);
% R4 = rand(1,1,25);
% y0 = 1 ...
%     + sum(A * R1 .* cos(V(1,1,:,1) .* x{1} + V(1,1,:,2) .* x{2}' + 2*pi*R2),3)...
%     + sum(A * R3 .* cos(V(1,1,:,1) .* x{1} - V(1,1,:,2) .* x{2}' + 2*pi*R4),3);

plot_surface(x,y0')

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

options = optimoptions('fsolve',...
   'SpecifyObjectiveGradient',true,...
   'Display','off');
optimmethod = @(fun,x0) fsolve(fun,x0,options);
% optimmethod = @newton;
timestepper = @(odefun,t,y0) bdf3(odefun,t,y0,optimmethod);

% evnt = @(~,H) is_dewetted(H);
% 
% options = odeset(...
%     ...'Vectorized','on',...
%     ...'BDF','on',... % Backward differentiation formulae
%     'Event',@(t,y) event(t,y),...
%     ...'OutputFcn','odeprint',...
%     'OutputFcn',@outputfunction,...
%     'RelTol',1e-3...% Default: 1e-3
%     ); 

%% Solve
tic
[y, t] = solver(pdefun, t, x, y0, method, timestepper);
timeTaken = toc

%% Plot overview

plot_surface(x,real(y(:,:,end))');

figure
plot_log_fourier(x,real(y(:,:,end))');

figure
plot(t,log10(squeeze(energy(x,y))));

figure
plot(t,squeeze(min(y,[],[1,2])));

[temp,xI] = min(y);
[miny,yI] = min(temp);
xI = xI(yI);

figure
plot(t,squeeze(x{1}(xI)),t,squeeze(x{2}(yI)));


%% H2 norm - Check implemented correctly
ycell = mat2cell(reshape(y,128^2,409)',ones(409,1));

plot(t,log10(cellfun(@(y) norm_H2(x,y',method),ycell)));

%% Construct other variables
y_snapshot = y0;
% y_snapshot = y(:,:,end);

p = construct_pressure(x,y_snapshot,params,method);

[z,u,v,w] = construct_velocity(x,y_snapshot,params,method);

slice = 1;

surf(repmat(x{1},1,size(z(:,slice,:),3)),squeeze(u(:,slice,:)),squeeze(z(:,slice,:)));
figure
surf(repmat(x{1},1,size(z(:,slice,:),3)),squeeze(v(:,slice,:)),squeeze(z(:,slice,:)));
figure
surf(repmat(x{1},1,size(z(:,slice,:),3)),squeeze(w(:,slice,:)),squeeze(z(:,slice,:)));
figure
quiver(squeeze(w(:,slice,:))',squeeze(v(:,slice,:))');

%% Animate

addpath('../plotting/')

F = animate(@(x,y) plot_surface(x,y'),t,x,y);

% movie(F);

%% Save

filename = replace(sprintf('data-theta-%f-Re-%f-We-%f-C-%f-xL-%f-yL-%f-T-%f',[params(2:end),xL,tL]),'.','_');
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