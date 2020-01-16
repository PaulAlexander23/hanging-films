function odeopt = setupIVPodeopt(args, domain)
    odeopt = odeset();
    
    if args.domainArguments.method == "finite-difference"
        if args.model == "benney"
            odeopt = odeset(odeopt, ...
                'Jacobian', @(t, y) jbenney2d(domain, y, args.params) ...
                );
        end
    end
end