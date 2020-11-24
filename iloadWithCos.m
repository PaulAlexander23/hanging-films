function y0 = iloadWithCos(x, filename, a, b)
    if nargin < 2, a = 0.25; end
    if nargin < 3, b = 0.25; end

    load(filename,'h');
    y0 = h + icos(x, a, b) - 1;
end
