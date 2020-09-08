function domain = setupDomain(args)
    
    addpath('discretisationMethods/');
    
    x = setupX(args.xLength, args.yLength, args.xN, args.yN);

    if args.method == "finite-difference"
        problemDiffDegrees = [1, 0; 2, 0; 0, 1; 0, 2]';
        accuracy = 4;
        direction = "central";
        domain = FDDomain(x, problemDiffDegrees, accuracy, direction);
    elseif args.method == "pseudo-spectral"
        AA = true;
        complex = true;
        domain = PSDomain(x, AA, complex);
    end
end
