function y0 = iloadWithCos(x, filename, a, b, c, d)
    if nargin < 3, a = 0.25; end
    if nargin < 4, b = 0.25; end
    if nargin < 5, c = 0.0; end
    if nargin < 6, d = 0.0; end

    load(filename,'h');

    y0 = h + [icos(x, a, b); icos(x, c, d) - 1/3];
end
