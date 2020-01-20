function y0 = irivulet(x, A, r)
    
    if nargin < 2, A = 2e-1; end
    if nargin < 3, r = 0.05; end
    y0 = 1 + A * (-r*cos(2*pi/x{1}(end) * x{1}) - cos(2*pi/x{2}(end) * x{2}));
    
end