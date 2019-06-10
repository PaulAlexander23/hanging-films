function saveData(y, params, t, x, timeTaken, tFinal, interface, AbsTol, prefix)
    if nargin < 9, prefix = ""; end
    filename = makeFilename(prefix, params, x, tFinal, interface, AbsTol);
    save(filename, 'y', 'params', 't', 'x', 'timeTaken');
end