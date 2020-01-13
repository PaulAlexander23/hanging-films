function main(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method, debug)
    addpath discretisationMethods
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, AbsTol = 1e-6; end
    if nargin < 12, method = "finite-difference"; end
    if nargin < 13, debug = false; end
    
    params = paramsToStruct(theta, Re, C);
    
    domainArguments = struct('xLength', xLength, 'yLength', yLength, 'xN', xN, ...
        'yN', yN, 'method', method);
    
    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'method',method,'model',model,'interface',interface);
    
    timePointsArguments = struct('tStep', 0.2, 'tFinal', tFinal);
    
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'AbsTol', AbsTol ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        );
    timeStepperArguments = struct('timeStepper', @ode15s, 'odeopt', odeoptDefault);
    
    [domain, y, t, timeTaken] = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments, debug);
    
    saveData(y, params, t, domain.x, timeTaken, timePointsArguments.tFinal, interface, AbsTol, model)
end
