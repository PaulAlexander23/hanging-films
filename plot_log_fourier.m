function plot_log_fourier(x,y)

    % Transform into fourier space
    yf = fft2(y-1);
    
    Nx = size(yf,1)/2;
    Ny = size(yf,2)/2;

    % Determine k in matlab form
    kx = [0:Nx-1, 0, 1-Nx:-1]' * 2*pi/x{1}(end);
    ky = [0:Ny-1, 0, 1-Ny:-1]' * 2*pi/x{2}(end);

    [X,Y] = meshgrid(ky,kx);

    surf(X,Y,log10(abs(yf))/Nx/Ny)
    
    shading interp
end
