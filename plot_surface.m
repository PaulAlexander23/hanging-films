function plot_surface(x,y,H)
    figure
    [X, Y] = meshgrid(x,y);
    surf(X,Y,H');
    xlabel('x')
    ylabel('y')
    shading interp
end