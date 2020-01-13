function main(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method, debug)
    addpath discretisationMethods
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, AbsTol = 1e-6; end
    if nargin < 12, method = "finite-difference"; end
    if nargin < 13, debug = false; end
    
    params = paramsToStruct(theta, Re, C);
    domain = createDomain(xLength, yLength, xN, yN, method);
    
    timePointsArguments = struct('tStep', 0.2, 'tFinal', tFinal);
    
    [y, t, timeTaken] = solveIVP(model, domain, params, timePointsArguments, interface, method, AbsTol, debug);
    
    saveData(y, params, t, domain.x, timeTaken, timePointsArguments.tFinal, interface, AbsTol, model)
end
