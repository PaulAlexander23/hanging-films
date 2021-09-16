function y0 = irivuletWIBL2(x, rivulet_amplitude, number_of_rivulets, noise_amplitude)
    
    if nargin < 2, rivulet_amplitude = 1e-1; end
    if nargin < 3, number_of_rivulets = 1; end
    if nargin < 4, noise_amplitude = 0; end

    h0 = 1 + 0 * x{1} - rivulet_amplitude * cos(2*pi/x{2}(end) * x{2} * number_of_rivulets);
    f0 = 2/3 * h0;
    g0 = 0 * x{1} + 0 * x{2};
    
    y0 = [h0; f0; g0];

    rng(0);
    noise = noise_amplitude * (rand(size(y0)) - 1);

    y0 = y0 + noise;
end
