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
    
    function domain = createDomain(xLength, yLength, xN, yN, method)
        x = setupX(xLength, yLength, xN, yN);
        
        if method == "finite-difference"
            problemDiffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
            domain = FDDomain(x, problemDiffDegrees, 4);
        elseif method == "pseudo-spectral"
            domain = PSDomain(x);
        end
    end
end
