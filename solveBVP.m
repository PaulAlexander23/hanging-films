function solution = solveBVP(bvpArguments)
    
    [domain, odeFunction, y0, bvpOdeopt] = setupBVP(bvpArguments);
    
    [y, timeTaken] = iterateTimeStepper(odeFunction, t, y0, nonlinearSolver);

    y = postprocess(bvpArguments.domainArguments.method, domain, y);
    
    solution = struct('domain', domain, 't', t, 'y', y, ...
        'timeTaken', timeTaken);
end
