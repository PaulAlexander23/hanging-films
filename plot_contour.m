function plot_contour(x,y)
    [X, Y] = meshgrid(x{1},x{2});
    contourf(X,Y,y);
    xlabel('x')
    ylabel('y')
    shading interp
end