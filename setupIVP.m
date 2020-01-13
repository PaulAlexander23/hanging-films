function [domain, odeFunction, y0, odeopt] = setupIVP(args)
    
    domain = setupDomain(args.domainArguments);
    
    odeFunction = setupODEFunction(args.model, domain, args.params);
    
    y0 = setupInitialCondition(args.model, domain, args.interface, args.method);
    
    odeopt = setupIVPodeopt(args, domain);
end