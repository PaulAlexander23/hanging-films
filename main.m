%MAIN

%% Setup

xN = 2^7;
xL = 10*pi;
xS = xL/xN;
x = linspace(xS,xL,xN)';

yN = 2^5;
yL = 1*pi;
yS = yL/yN;
y = linspace(yS,yL,yN);

H0 = 1 ...
    + 0.10*rand() * cos(4*y + pi*rand())...
    + 0.02*rand() * cos(x + pi*rand())...
    + 0.05*rand() * cos(x + 4*y + rand())...
    + 0.05*rand() * cos(x - 4*y + rand());

% delta, theta, Re, We, C
params = [1,pi/2,1,0,1];

tFinal = 1;


%% Solve

tic
[H, ~, t] = solver(@f_benney, params, H0, tFinal, [xL, yL], [xN, yN], @(~,H) is_dewetted(H), 1e-3);
toc

%% Plot overview

plot_surface(x,y,H(:,:,1));

plot_surface(x,y,H(:,:,end));

figure
plot(t, squeeze(sum((H-1).^2 * xS * yS,[1,2])))

%% Plot specific

plot_surface(x,y,H(:,:,90));


%% Save

filename = replace(sprintf('data-%g-%g-%g-%g-%g.mat',params),'.','_');
save(filename,'H','H0','params','t','x','y');