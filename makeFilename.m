function filename = makeFilename(prefix, params, x, tFinal, interface, AbsTol, model)
    filename = replace(sprintf("data" + prefix + "-theta-%g-Re-%g-C-%g-xL-%g-yL-%g-T-%g-interface-%s-xN-%g-yN-%g-AbsTol-%g-model-%s", ...
        params(2), params(3), params(4), x{1}(end), x{2}(end), tFinal, func2str(interface), length(x{1}), length(x{2}), AbsTol, model), '.', '_');
end