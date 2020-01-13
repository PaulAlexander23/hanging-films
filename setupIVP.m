function [domain, odeFunction, y0, odeopt] = setupIVP(args)
    
    domain = setupDomain(args.domainArguments);
    
    odeFunction = setupODEFunction(args.model, domain, args.params);
    
    odeopt = odeset();
    if args.method == "finite-difference"
        if args.model == "benney"
            odeopt = odeset(odeopt, ...
                'Jacobian', @(t, y) jbenney2d(domain, y, args.params) ...
                );
        end
    end
    
    y0 = setupInitialCondition(args.model, domain, args.interface, args.method);
end