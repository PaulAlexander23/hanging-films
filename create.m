function create(theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol)
    addpath discretisationMethods
    if nargin < 8, xN = 64; end
    if nargin < 9, yN = 64; end
    if nargin < 10, AbsTol = 1e-6; end

    x = setupX(xLength, yLength, xN, yN);

    t = setupT(tFinal, 0.2);
    
    y0 = interface(x);

    params = [1, theta, Re, C]; % delta, theta, Re, C
    problemDiffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';

    domain = FDDomain(x, problemDiffDegrees, 4);
    pdeFunction = @(t, domain, y) fbenney2d(domain, y, params);

    odeopt = odeset( ...
        'Jacobian', @(t, y) jbenney2d(domain, y, params), ...
        ...'Vectorized', 'on', ...
        'AbsTol', AbsTol ...
        ... 'BDF','on' ...
        );
    timeStepper = @(odefun, t, y0) ode15s(odefun, t, y0, odeopt);

    tic
    [y, t] = pdeSolver(pdeFunction, t, domain, y0, timeStepper);
    timeTaken = toc;

    saveData(y, params, t, domain.x, timeTaken, tFinal, interface, AbsTol)
end