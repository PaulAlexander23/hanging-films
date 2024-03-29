function domain = setupDomain(args)
    
    addpath('discretisationMethods/');
    
    x = setupX(args.xLength, args.yLength, args.xN, args.yN);

    if args.method == "finite-difference"
        problemDiffDegrees = combvec(0:3,0:3);
        accuracy = 4;
        direction = "central";
        domain = FDDomain(x, problemDiffDegrees, accuracy, direction);
    elseif args.method == "pseudo-spectral"
        AA = true;
        complex = false;
        domain = PSDomain(x, AA, complex);
    end
end
