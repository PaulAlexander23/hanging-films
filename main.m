%MAIN

setup;

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
