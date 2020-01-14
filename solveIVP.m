function [domain, y, t, timeTaken] = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments)
    
    [domain, odeFunction, y0, odeoptIVP] = setupIVP(ivpArguments);
    
    t = setupTimePoints(timePointsArguments);
    
    timeStepper = setupTimeStepper(timeStepperArguments, odeoptIVP);
    
    [y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper);

    y = postprocess(ivpArguments.method, ivpArguments.model, domain, y);
end