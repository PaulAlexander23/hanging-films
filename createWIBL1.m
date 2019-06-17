function createWIBL1(theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol)
    import discretisationMethods.*
    if nargin < 8, xN = 64; end
    if nargin < 9, yN = 64; end
    if nargin < 10, AbsTol = 1e-6; end
    
    x = setupX(xLength, yLength, xN, yN);
    
    tStep = 0.2;
    t = 0:tStep:tFinal;
    
    y0 = interface(x);
    y0 = [y0; 0*y0]; % [y0; F_0]
    
    params = [1, theta, Re, C]; % delta, theta, Re, C
    problemDiffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    
    diffMat = init_fd(x, problemDiffDegrees, 4);
    diffMethod = @(x, y, degree) diff_fd(x, y, degree, diffMat, problemDiffDegrees);
    pdeFunction = @(t, x, y, diffMethod) fwibl1(x, y, params, diffMethod);
    
    odeopt = odeset( ...
        ... 'Jacobian', @(t, y) jbenney(x, y, params, method, getD), ...
        ... 'Vectorized', 'on', ...
        'AbsTol', AbsTol ...
        ... 'BDF','on' ... % Backward differentiation formulae
        );
    timeStepper = @(odefun, t, y0) ode15s(odefun, t, y0, odeopt);
    
    tic
    [y, t] = pdeSolver(pdeFunction, t, x, y0, diffMethod, timeStepper);
    timeTaken = toc;
    
    % F = y(end/2+1:end, :, :);
    % y = y(1:end/2, :, :);
    
    saveData(y, params, t, x, timeTaken, tFinal, interface, AbsTol, "-wibl1")
end