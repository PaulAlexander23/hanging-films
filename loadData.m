function [y, params, t, x, timeTaken] = loadData(params, x, tFinal, interface, AbsTol, prefix)
    if nargin < 6, prefix = ""; end
    filename = makeFilename(prefix, params, x, tFinal, interface, AbsTol);
    load(filename, 'y', 'params', 't', 'x', 'timeTaken');
end