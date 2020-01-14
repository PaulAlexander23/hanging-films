function solution = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments)
    
    [domain, odeFunction, y0, odeoptIVP] = setupIVP(ivpArguments);
    
    t = setupTimePoints(timePointsArguments);
    
    timeStepper = setupTimeStepper(timeStepperArguments, odeoptIVP);
    
    [y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper);

    y = postprocess(ivpArguments.domainArguments.method, ivpArguments.model, domain, y);
    
    solution = struct('domain', domain, 't', t, 'y', y, ...
        'timeTaken', timeTaken);
end