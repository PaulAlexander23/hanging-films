function y0 = iloadAmplifyModeWIBL1(x, filename, amplitude, modeX, modeY)
    if nargin < 5, modeY = 0; end
    if nargin < 4, modeX = 0; end
    if nargin < 3, amplitude = 0.25; end

    kX = fftshift(-length(x{1})/2:length(x{1})/2-1);
    kY = fftshift(-length(x{2})/2:length(x{2})/2-1);

    y = iload(x, filename);

    h0 = y(1:end/2,:);
    F0 = y(1+end/2:end,:);

    h0 = fft2(h0);
    h0(kX==modeX, kY==modeY) = (1 + amplitude) * h0(kX==modeX, kY==modeY);
    h0 = ifft2(h0, 'symmetric');

    F0 = fft2(F0);
    F0(kX==modeX, kY==modeY) = (1 + amplitude) * F0(kX==modeX, kY==modeY);
    F0 = ifft2(F0, 'symmetric');

    y0 = [h0; F0];
end
