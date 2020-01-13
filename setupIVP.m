function [odeFunction, y0, odeopt] = setupIVP(model, domain, params, interface, method)
    
    odeFunction = setupODEFunction(model, domain, params);
    
    odeopt = odeset();
    if method == "finite-difference"
        if model == "benney"
            odeopt = odeset(odeopt, ...
                'Jacobian', @(t, y) jbenney2d(domain, y, params) ...
                );
        end
    end
    
    y0 = setupInitialCondition(model, domain, interface, method);
end