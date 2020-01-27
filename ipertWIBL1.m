function y0 = ipertWIBL1(x, A, r, B, s)
    if nargin < 2, A = 1e-2; end
    if nargin < 3, r = 1e-1; end
    if nargin < 4, B = 1e-2; end
    if nargin < 5, s = 1e-1; end
    
    y0 = [1 + A * exp(-r*((x{1}-x{1}(end)/2).^2 + (x{2}-x{2}(end)/2).^2));
        2/3 + B * exp(-s*((x{1}-x{1}(end)/2).^2 + (x{2}-x{2}(end)/2).^2))];
end