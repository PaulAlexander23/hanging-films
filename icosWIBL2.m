function y0 = icosWIBL2(x, a, b, c, d, e, f)
    if nargin < 2, a = 0.25; end
    if nargin < 3, b = 0.25; end
    if nargin < 4, c = 0; end
    if nargin < 5, d = 0; end
    if nargin < 6, e = 0; end
    if nargin < 7, f = 0; end
    
    y0 = [1 - a * cos(2*pi/x{1}(end) * x{1}) - ...
        b * cos(2*pi/x{2}(end) * x{2});
        2/3 - c * cos(2*pi/x{1}(end) * x{1}) - ...
        d * cos(2*pi/x{2}(end) * x{2});
        - e * cos(2*pi/x{1}(end) * x{1}) - ...
        f * cos(2*pi/x{2}(end) * x{2})];
end
