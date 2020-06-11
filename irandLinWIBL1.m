function y0 = irandLinWIBL1(x, num, A, B, seed)
    %IRANDWIBL1 Interface on x with "num" random modes and amplitudes "A" and 
    % "B" for the surface and flux respectively. Uses seed for random number 
    % generation.
    if nargin < 2, num = 5; end
    if nargin < 3, A = 1e-4; end
    if nargin < 4, B = 1e-4; end
    if nargin == 5
        rng(seed);
    end

    Nx = length(x{1})/2;
    Ny = length(x{2})/2;

    % Determine k in matlab form
    kx = [0:floor(Nx-1), -ceil(Nx):-1]';
    ky = [0:floor(Ny-1), -ceil(Ny):-1];
    
    fh0 = zeros(2*Nx,2*Ny);
    
    fh0(logical((kx~=0).*(abs(kx)<num+1)), ...
        logical((ky~=0).*(abs(ky)<num+1))) = A*exp(2)/sqrt(2) * ...
        (2*rand(2*num,2*num)-1 + 2i*rand(2*num,2*num)-1i);

    fh0 = exp(-abs(ky)) .* fh0 .* exp(-abs(kx));

    fh0 = fh0 * Nx * Ny;
    
    h0 = ifft2(fh0, 'symmetric');
    h0 = 1 + circshift(circshift(h0, -1, 1), -1, 2);
    
    ff0 = zeros(2*Nx,2*Ny);
    
    ff0(logical((kx~=0).*(abs(kx)<num+1)), ...
        logical((ky~=0).*(abs(ky)<num+1))) = B*exp(2)/sqrt(2) * ...
        (2*rand(2*num,2*num)-1 + 2i*rand(2*num,2*num)-1i);
    
    ff0 = exp(-abs(ky)) .* ff0 .* exp(-abs(kx));

    ff0 = ff0 * Nx * Ny;
    
    f0 = ifft2(ff0, 'symmetric');
    f0 = 2/3 + circshift(circshift(f0, -1, 1), -1, 2);

    y0 = [h0; f0];
end
