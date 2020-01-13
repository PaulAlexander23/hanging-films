function [odeFunction, y0] = setupIVP(model, domain, params, interface, method)
    
    odeFunction = setupODEFunction(model, domain, params);
    
    y0 = setupInitialCondition(model, domain, interface, method);
end