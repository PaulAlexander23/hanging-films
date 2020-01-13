function [y, t, timeTaken] = solveIVP(model, domain, params, timePointsArguments, interface, method, timeStepperArguments, debug)
    
    [odeFunction, y0, odeoptIVP] = setupIVP(model, domain, params, interface, method);
    
    t = setupTimePoints(timePointsArguments, debug);
    
    timeStepper = setupTimeStepper(timeStepperArguments, odeoptIVP, debug);
    
    [y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper);

    y = postprocess(method, model, domain, y);
end