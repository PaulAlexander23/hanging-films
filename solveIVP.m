function [y, t, timeTaken] = solveIVP(model, domain, params, timePointsArguments, interface, method, AbsTol, debug)
    
    [odeFunction, y0] = setupIVP(model, domain, params, interface, method);
    
    t = setupTimePoints(timePointsArguments, debug);
    
    timeStepper = setupTimeStepper(method, model, domain, params, AbsTol, debug);
    
    [y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper);

    y = postprocess(method, model, domain, y);
end