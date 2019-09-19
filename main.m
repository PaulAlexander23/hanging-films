function main(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method, debug)
    addpath discretisationMethods
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, AbsTol = 1e-6; end
    if nargin < 12, method = "finite-difference"; end
    if nargin < 13, debug = false; end
    
    params = struct('theta', theta, 'Re', Re, 'C', C);
    domain = createDomain(xLength, yLength, xN, yN, method);
    
    [y, t, timeTaken] = createData(model, domain, params, tFinal, interface, method, AbsTol, debug);
    
    saveData(y, params, t, domain.x, timeTaken, tFinal, interface, AbsTol, model)
end
