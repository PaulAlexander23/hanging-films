function y0 = ipert(x, A, r)
    
    if nargin < 2, A = 1e-2; end
    if nargin < 3, r = 1e-1; end
    y0 = 1 + A * exp(-r*((x{1}-x{1}(end)/2).^2 + (x{2}'-x{2}(end)/2).^2));
    
end