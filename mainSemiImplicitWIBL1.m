function mainSemiImplicitWIBL1(icFilename)

    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    xL = 32; yL = 32; xN = 64; yN = 64;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    Re = replace(icFilename, 'ics/', '');
    Re = replace(Re, 'ic-WIBL1-2D-Re-', '');
    Re = replace(Re, 'ic-WIBL1-3D-Re-', '');
    Re = replace(Re, '.mat', '');
    Re = replace(Re, '_', '.');
    Re = str2num(Re)
    
    params = struct('theta', 7*pi/8, 'Re', Re, 'C', 0.01);

    fwibl1ExplicitVec = matFuncToVecFunc(@fwibl1Explicit);
    explicitOdefun = @(t, y) fwibl1ExplicitVec(domain, y, params);

    t = linspace(0, 300,  40000)';

    load(icFilename, 'h');
    y0 = domain.reshapeToVector(h);

    timeStepper = @bdf1si;

    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true)));

    tic
    [~, y] = timeStepper(explicitOdefun, @(t, y) implicitOdefun(t, y, domain, params), t, y0, options);
    timeTaken = toc;
    y = y';
    y = permute(y, [1, 3, 2]);
    y = domain.reshapeToDomain(y);

    filename = replace(replace(icFilename, 'ics/', ''), 'ic', 'data');

    save(filename, 't', 'x', 'y', 'params', 'timeTaken', '-v7.3')

    function [F, J] = implicitOdefun(~, y, domain, params)
        [f, J] = fwibl1Implicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end
