function y0 = irandWIBL1(x, num, A, B)
    %IRANDWIBL1 Interface on x with "num" random modes and amplitudes "A" and 
    % "B" for the surface and flux respectively.
    if nargin < 2, num = 5; end
    if nargin < 3, A = 1e-4; end
    if nargin < 4, B = 1e-4; end

    Nx = length(x{1})/2;
    Ny = length(x{2})/2;

    % Determine k in matlab form
    kx = [0:Nx-1, 0, 1-Nx:-1]';
    ky = [0:Ny-1, 0, 1-Ny:-1]';
    
    fh0 = zeros(2*Nx,2*Ny);
    
    fh0(logical((kx~=0).*(abs(kx)<num+1)), ...
        logical((ky~=0).*(abs(ky)<num+1))) = (2*rand(2*num,2*num)-1)*Nx*Ny;
    
    h0 = real(ifft2(fh0));
    h0 = 1 + A * h0/max(max(abs(h0)));
    
    ff0 = zeros(2*Nx,2*Ny);
    
    ff0(logical((kx~=0).*(abs(kx)<num+1)), ...
        logical((ky~=0).*(abs(ky)<num+1))) = (2*rand(2*num,2*num)-1)*Nx*Ny;
    
    f0 = real(ifft2(ff0));
    f0 = 2/3 + B * f0/max(max(abs(f0)));

    y0 = [h0; f0];
end
