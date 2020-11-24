function [domain, odeFunction, y0, ivpOdeopt] = setupIVP(args)
    
    domain = setupDomain(args.domainArguments);
    
    if isstruct(args.odefun)
        odeFunction = struct();
        fields = fieldnames(args.odefun);
        for n = 1:numel(fields)
            odefun = getfield(args.odefun,fields{n});
            odeFunction = setfield(odeFunction, fields{n}, ...
                @(t, y) odefun(domain, y, args.params)); 
        end
    else
        odeFunction = @(t, y) args.odefun(domain, y, args.params);
    end
    
    y0 = setupInitialCondition(domain, args.interface, args.domainArguments.method);
    
    ivpOdeopt = setupIVPodeopt(args.odejac, domain, args.params, args.domainArguments.method);
end
