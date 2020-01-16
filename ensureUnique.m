function filename = ensureUnique(filename)
    n = 0;
    suffix = "";
    while isfile(filename + suffix + ".mat")
        n = n + 1;
        suffix = "-" + num2str(n);
    end
    filename = filename + suffix;
end