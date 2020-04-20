function domain = setupDomain(args)
    
    addpath('discretisationMethods/');
    
    x = setupX(args.xLength, args.yLength, args.xN, args.yN);

    if args.method == "finite-difference"
        problemDiffDegrees = [1, 0; 2, 0; 3, 0; 4, 0; 1, 1; 2, 1; 3, 1; 1, 2; 2, 2; 1, 3; 0, 1; 0, 2; 0, 3; 0, 4]';
        accuracy = 4;
        domain = FDDomain(x, problemDiffDegrees, accuracy);
    elseif args.method == "pseudo-spectral"
        AA = true;
        complex = true;
        domain = PSDomain(x, AA, complex);
    end
end
