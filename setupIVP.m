function [domain, odeFunction, y0, odeopt] = setupIVP(args)
    
    domain = setupDomain(args.domainArguments);
    
    odeFunction = @(t, y) args.odefun(domain, y, args.params);
    
    y0 = setupInitialCondition(domain, args.interface, args.domainArguments.method);
    
    odeopt = setupIVPodeopt(args, domain);
end
