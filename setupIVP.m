function [domain, odeFunction, y0, ivpOdeopt] = setupIVP(args)
    
    domain = setupDomain(args.domainArguments);
    
    odeFunction = @(t, y) args.odefun(domain, y, args.params);
    
    y0 = setupInitialCondition(domain, args.interface, args.domainArguments.method);
    
    ivpOdeopt = setupIVPodeopt(args.odejac, domain, args.params, args.domainArguments.method);
end
