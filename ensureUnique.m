function filename = ensureUnique(filename)
    n = 0;
    suffix = "";
    while exist(filename + suffix + ".mat", 'file')
        n = n + 1;
        suffix = "-" + num2str(n);
    end
    filename = filename + suffix;
end