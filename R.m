function out = R(H, L, suppression)
    if nargin < 3
       suppression = 1e-13;
    end
    
    % Transform into fourier space
    HF = fft2(H);
    
    Nx = size(HF,1)/2;
    Ny = size(HF,2)/2;

    % Determine k in matlab form
    kx = [0:Nx-1, 0, 1-Nx:-1]' * 2*pi/L(1);
    ky = [0:Ny-1, 0, 1-Ny:-1]' * 2*pi/L(2);
    
    % Prior suppression
    HF(abs(HF)<suppression) = 0;

    % Apply fractional Laplacian
    [Kx,Ky] = meshgrid(ky,kx);
    
    RHF = sqrt(Kx.^2 + Ky.^2) .* HF;

    % Posterior suppression
    % HFx(abs(HFx) < suppression*N*2) = 0 ;
    % HFx(abs(HFx) < suppression*max(abs(HFx))) = 0 ;
    % HFy(abs(HFy) < suppression*N*2) = 0 ;
    % HFy(abs(HFy) < suppression*max(abs(HFy))) = 0 ;
    
    % Transform back into real space
    out = ifft2(RHF);
end