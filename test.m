%TEST Testing code for the pseudo-spectral Benney-type equation solver.

reset;
clc;

%% Setup
xN = 2^7;
xL = 10*pi;
xS = xL/xN;
x = linspace(xS,xL,xN)';

yN = 2^5;
yL = 1*pi;
yS = yL/yN;
y = linspace(yS,yL,yN);

L = [xL, yL];

H0 = 1 + ...
    0.1*rand() * (cos(4*y + pi*rand()) + 1.1) + ...
    0.02*rand() * (cos(x + pi*rand()) + 1.1) + ...
    0.05*rand() * (cos(x + 4*y + rand()) + 1.1) + ...
    0.05*rand() * (cos(x - 4*y + rand()) + 1.1);

%%

plot_surface(x,y,H0)

gradH = grad(H0,L);

plot_surface(x,y,gradH(:,:,1));

plot_surface(x,y,gradH(:,:,2));

lapH = lap(H0,L);

plot_surface(x,y,lapH);

RH = R(H0,L);

plot_surface(x,y,RH);

F = f_benney(H0,L,[1,1,1,1,1]);

plot_surface(x,y,F);

%%

tic
[H, ~, t] = solver(@f_benney, params, H0, tFinal, L, [xN, yN], @(~,H) is_dewetted(H), 1e-3);
toc