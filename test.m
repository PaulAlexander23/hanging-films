%TEST Testing code for the pseudo-spectral Benney-type equation solver.

clear;
clc;

setup;

L = [xL, yL];

H1 = 1 + sin(2*pi*x/xL + 4*pi*y/yL);

%%

plot_surface(x,y,H0)

plot_surface(x,y,H1)

%% Grad

gradH0 = grad(H0,L);

plot_surface(x,y,gradH0(:,:,1));
plot_surface(x,y,gradH0(:,:,2));

gradH1 = grad(H1,L);

plot_surface(x,y,gradH1(:,:,1));
plot_surface(x,y,gradH1(:,:,2));

%% Laplacian

lapH0 = lap(H0,L);

plot_surface(x,y,lapH0);

lapH1 = lap(H1,L);

plot_surface(x,y,lapH1);

%% R

RH0 = R(H0,L);

plot_surface(x,y,RH0);

RH1 = R(H1,L);

plot_surface(x,y,RH1);

%% F

F0 = f_benney(H0,L,[1,1,1,1,1]);

plot_surface(x,y,F0);

F1 = f_benney(H1,L,[1,1,1,1,1]);

plot_surface(x,y,F1);

%%% Linear
%
%tic
%[H, ~, t] = solver(@f_benney, params, H0, tFinal, L, [xN, yN], @(~,H) is_dewetted(H), 1e-3);
%toc
%
%%% Hanging film is independent of Re
%
%%% Lyapanov expontent
%
%
%
%%% Volume conservation
%
%%% Resolution decrease, increases accuracy
