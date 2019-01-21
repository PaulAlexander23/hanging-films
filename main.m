

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
% plot_surface(x,y,H0)
% 
% gradH = grad(H0,L);
% 
% plot_surface(x,y,gradH(:,:,1));
% 
% plot_surface(x,y,gradH(:,:,2));
% 
% lapH = lap(H0,L);
% 
% plot_surface(x,y,lapH);
% 
% RH = R(H0,L);
% 
% plot_surface(x,y,RH);
% 
% F = f_benney(H0,L,[1,1,1,1,1]);
% 
% plot_surface(x,y,F);

%%
params = [1,1,1,1,1];

tFinal = 10;

tic
[H, ~, t] = compute_numerical_solution2(@f_benney, params, H0, tFinal, L, [xN, yN], 1e-3);
toc

%%
plot_surface(x,y,H(:,:,1));

plot_surface(x,y,H(:,:,end));


%%

plot_surface(x,y,H(:,:,90));


%%
filename = replace(sprintf('data-%g-%g-%g-%g-%g.mat',params),'.','_');
save(filename,'H','H0','params','t','x','y');

%%
function plot_surface(x,y,H)
    figure
    [X, Y] = meshgrid(x,y);
    surf(X,Y,H');
    xlabel('x')
    ylabel('y')
    shading interp
end 