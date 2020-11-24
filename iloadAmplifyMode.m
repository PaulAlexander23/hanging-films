function y0 = iloadAmplifyMode(x, filename, amplitude, modeX, modeY)
    if nargin < 5, modeY = 0; end
    if nargin < 4, modeX = 0; end
    if nargin < 3, amplitude = 0.25; end

    h = iload(x, filename);

    y0 = fft2(h);

    kX = fftshift(-length(x{1})/2:length(x{1})/2-1);
    kY = fftshift(-length(x{2})/2:length(x{2})/2-1);

    y0(kX==modeX, kY==modeY) = (1 + amplitude) * y0(kX==modeX, kY==modeY);

    y0 = ifft2(y0, 'symmetric');
end
