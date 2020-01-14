function saveData(y, params, t, x, timeTaken, tFinal, interface, AbsTol, model, prefix)
    if nargin < 10, prefix = ""; end
    
    filename = makeFilename(prefix, params, x, tFinal, interface, AbsTol, model);
    
    save(filename, 'y', 'params', 't', 'x', 'timeTaken');
end