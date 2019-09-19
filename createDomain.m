function domain = createDomain(xLength, yLength, xN, yN, method)
    x = setupX(xLength, yLength, xN, yN);

    if method == "finite-difference"
        problemDiffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        domain = FDDomain(x, problemDiffDegrees, 4);
    elseif method == "pseudo-spectral"
        domain = PSDomain(x);
    end
end
