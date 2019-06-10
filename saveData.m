function saveData(y, params, t, x, timeTaken, tFinal, interface, AbsTol, prefix)
    if nargin < 9, prefix = ""; end
    filename = replace(sprintf("data" + prefix + "-theta-%g-Re-%g-C-%g-xL-%g-yL-%g-T-%g-interface-%s-xN-%g-yN-%g-AbsTol-%g", ...
        params(2), params(3), params(4), x{1}(end), x{2}(end), tFinal, func2str(interface), length(x{1}), length(x{2}), AbsTol), '.', '_');
    save(filename, 'y', 'params', 't', 'x', 'timeTaken');
end