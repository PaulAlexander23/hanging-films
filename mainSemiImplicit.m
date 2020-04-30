function mainSemiImplicit(icFilename)

    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    xL = 32; yL = 32; xN = 64; yN = 64;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);
    
    fbenney2dExplicitVec = matFuncToVecFunc(@fbenney2dExplicit);
    explicitOdefun = @(t, y) fbenney2dExplicitVec(domain, y, params);

    t = linspace(0, 1,  20)';

    load(icFilename, 'h');
    y0 = domain.reshapeToVector(h);

    timeStepper = @bdf2si;

    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true)));

    tic
    [~, y] = timeStepper(explicitOdefun, @(t, y) implicitOdefun(t, y, domain, params), t, y0, options);
    timeTaken = toc;
    y = y';
    y = permute(y, [1, 3, 2]);
    y = domain.reshapeToDomain(y);

    save(replace(icFilename, 'ic', 'data'), 't', 'y', 'timeTaken', 'params', 'x')

    function [F, J] = implicitOdefun(~, y, domain, params)
        [f, J] = fbenney2dImplicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end
