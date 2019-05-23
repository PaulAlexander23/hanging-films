function y0 = irand(x, A, num)
    %IRAND Interface on x with "num" random modes and amplitude "A".
    
    if nargin < 2, A = 1e-4; end
    if nargin < 3, num = 5; end
    Nx = length(x{1})/2;
    Ny = length(x{2})/2;

    % Determine k in matlab form
    kx = [0:Nx-1, 0, 1-Nx:-1]';
    ky = [0:Ny-1, 0, 1-Ny:-1]';
    
    fy0 = zeros(2*Nx,2*Ny);
    
    fy0(logical((kx~=0).*(abs(kx)<num+1)), ...
        logical((ky~=0).*(abs(ky)<num+1))) = (2*rand(2*num,2*num)-1)*Nx*Ny;
    
    y0 = real(ifft2(fy0));
    y0 = 1 + A * y0/max(max(abs(y0)));
    
end
