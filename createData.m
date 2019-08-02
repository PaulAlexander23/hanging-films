function createData(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method)
    addpath discretisationMethods
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, AbsTol = 1e-6; end
    if nargin < 12, method = "finite-difference"; end
    
    params = [1, theta, Re, C]; % delta, theta, Re, C
    t = setupT(tFinal, 0.2);
    x = setupX(xLength, yLength, xN, yN);
    
    if method == "finite-difference"
        problemDiffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        domain = FDDomain(x, problemDiffDegrees, 4);
    elseif method == "pseudo-spectral"
        domain = PSDomain(x);
    end
    
    if model == "benney"
        y0 = interface(domain.x);
        odeFunction = @(t, y) fbenney2d(domain, y, params);
        odeopt = odeset( ...
            'Jacobian', @(t, y) jbenney2d(domain, y, params), ...
            ...'Vectorized', 'on', ...
            'AbsTol', AbsTol ...
            ... 'BDF','on' ...
            );
    elseif model == "wibl1"
        y0 = interface(x);
        y0 = [y0; 0*y0]; % [y0; F_0]
        odeFunction = @(t, y) fwibl1(domain, y, params);
        odeopt = odeset( ...
            ... 'Vectorized', 'on', ...
            'AbsTol', AbsTol ...
            ... 'BDF','on' ...
            );
    end
    
    timeStepper = @(odefun, t, y0) ode15s(odefun, t, y0, odeopt);
    tic
    [y, t] = odeMatrixSolver(odeFunction, t, y0, timeStepper);
    timeTaken = toc;
    
    saveData(y, params, t, domain.x, timeTaken, tFinal, interface, AbsTol, model)
end