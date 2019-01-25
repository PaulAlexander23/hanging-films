function plot_fourier(H,L)

    % Transform into fourier space
    HF = fft2(H-1);
    
    Nx = size(HF,1)/2;
    Ny = size(HF,2)/2;

    % Determine k in matlab form
    kx = [0:Nx-1, 0, 1-Nx:-1]' * 2*pi/L(1);
    ky = [0:Ny-1, 0, 1-Ny:-1]' * 2*pi/L(2);

    [X,Y] = meshgrid(kx,ky);

    surf(X,Y,real(HF)/Nx/Ny)
    
    shading interp
end