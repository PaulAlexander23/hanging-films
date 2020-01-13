function [domain, y, t, timeTaken] = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments, debug)
    
    [domain, odeFunction, y0, odeoptIVP] = setupIVP(ivpArguments);
    
    t = setupTimePoints(timePointsArguments, debug);
    
    timeStepper = setupTimeStepper(timeStepperArguments, odeoptIVP, debug);
    
    [y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper);

    y = postprocess(ivpArguments.method, ivpArguments.model, domain, y);
end