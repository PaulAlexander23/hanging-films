function odeopt = setupIVPodeopt(args, domain)
    odeopt = odeset();
    
    if args.domainArguments.method == "finite-difference"
        odeopt = odeset(odeopt ...
            , 'Jacobian', @(t, y) args.odejac(domain, y, args.params) ...
            );
    end
end
