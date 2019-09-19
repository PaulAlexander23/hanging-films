function createData(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method, debug)
    addpath discretisationMethods
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, AbsTol = 1e-6; end
    if nargin < 12, method = "finite-difference"; end
    if nargin < 13, debug = false; end
    
    params = struct('theta', theta, 'Re', Re, 'C', C);
    domain = createDomain(xLength, yLength, xN, yN, method);
    
    odeFunction = createODEFunction(model, domain, params);
    
    t = createT(tFinal, debug);
    
    y0 = createInitialCondition(model, domain, interface, method);
    
    timeStepper = createTimeStepper(method, model, domain, params, AbsTol, debug);
    
    tic
    [y, t] = odeMatrixSolver(odeFunction, t, y0, timeStepper);
    timeTaken = toc;
    
    y = postProcessing(method, model, domain, y);
    
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
    
    function odeFunction = createODEFunction(model, domain, params)
        if model == "benney"
            odeFunction = @(t, y) fbenney2d(domain, y, params);
        elseif model == "wibl1"
            odeFunction = @(t, y) fwibl1(domain, y, params);
        end
    end
    
    function t = createT(tFinal, debug)
        if ~debug
            t = setupT(tFinal, 0.2);
        else
            t = [0,tFinal];
        end
    end
    
    function y0 = createInitialCondition(model, domain, interface, method)
        if model == "benney"
            y0 = interface(domain.x);
        elseif model == "wibl1"
            y0 = interface(domain.x);
            F0 = 2/3 + 0*y0;
            y0 = [y0; F0];
        end
        
        
        if method == "pseudo-spectral"
            if model == "benney"
                y0 = domain.fft(y0);
            elseif model == "wibl1"
                y0 = [domain.fft(y0(1:end/2,:)); ...
                    domain.fft(y0(1+end/2:end,:))];
            end
        end
    end
    
    function timeStepper = createTimeStepper(method, model, domain, params, AbsTol, debug)
        odeopt = odeset( ...
            ...'Vectorized', 'on', ...
            ...'BDF','on', ...
            'AbsTol', AbsTol ...
            ...'MaxStep', 1e-3, ...
            ...'InitialStep', 1e-3 ...
            );
        if debug
            odeopt = odeset(odeopt, ...
                'OutputFcn', 'odeprint', ...
                'OutputSel', 1);
        end
        
        if method == "finite-difference"
            if model == "benney"
                odeopt = odeset(odeopt, ...
                    'Jacobian', @(t, y) jbenney2d(domain, y, params) ...
                    );
            end
        end
        
        if method == "finite-difference"
            timeStepper = @(odefun, t, y0) ode15s(odefun, t, y0, odeopt);
        elseif method == "pseudo-spectral"
            timeStepper = @(odefun, t, y0) ode45(odefun, t, y0, odeopt);
        end
    end
    
    function y = postProcessing(method, model, domain, y)
        if method == "pseudo-spectral"
            if model == "benney"
                y = domain.ifft(y);
            elseif model == "wibl1"
                y = [domain.ifft(y(1:end/2,:,:)); ...
                    domain.ifft(y(1+end/2:end,:,:))];
            end
        end
    end
end
