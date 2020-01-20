function y0 = icos(x, a, b)

    if nargin < 2, a = 0.25; end
    if nargin < 3, b = 0.25; end
    y0 = 1 - a * cos(2*pi/x{1}(end) * x{1}) - ...
        b * cos(2*pi/x{2}(end) * x{2});
    
end