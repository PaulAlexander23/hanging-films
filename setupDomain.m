function domain = setupDomain(args)
    
    addpath('discretisationMethods/');
    
    x = setupX(args.xLength, args.yLength, args.xN, args.yN);

    if args.method == "finite-difference"
        problemDiffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        accuracy = 4;
        domain = FDDomain(x, problemDiffDegrees, accuracy);
    elseif args.method == "pseudo-spectral"
        AA = true;
        complex = false;
        domain = PSDomain(x, AA, complex);
    end
end