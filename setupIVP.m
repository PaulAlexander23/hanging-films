function [domain, odeFunction, y0, ivpOdeopt] = setupIVP(args)
    
    domain = setupDomain(args.domainArguments);
    
    % Return domain and params with function handle.
    % 
    % Check if the function has parts which need to be treated differently in
    % the time stepper.
    if isstruct(args.odefun)
        odeFunction = struct();
        fields = fieldnames(args.odefun);
        for n = 1:numel(fields)
            odefun = getfield(args.odefun,fields{n});
            % odefunction has multiple parts so is passed as a structure.
            odeFunction = setfield(odeFunction, fields{n}, ...
                @(t, y) odefun(domain, y, args.params)); 
        end
    else
        % odefunction has a single handle so returned as a function.
        odeFunction = @(t, y) args.odefun(domain, y, args.params);
    end
    
    y0 = setupInitialCondition(domain, args.interface, args.domainArguments.method);
    
    ivpOdeopt = setupIVPodeopt(args.odejac, domain, args.params, args.domainArguments.method);
end
