function plot_surface(x,y)
    [X, Y] = meshgrid(x{1},x{2});
    surf(X,Y,y);
    xlabel('x')
    ylabel('y')
    shading interp
end