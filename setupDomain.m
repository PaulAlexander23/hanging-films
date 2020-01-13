function domain = setupDomain(args)
    
    addpath("discretisationMethods/");
    
    x = setupX(args.xLength, args.yLength, args.xN, args.yN);

    if args.method == "finite-difference"
        problemDiffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        domain = FDDomain(x, problemDiffDegrees, 4);
    elseif args.method == "pseudo-spectral"
        domain = PSDomain(x, true, false);
    end
end