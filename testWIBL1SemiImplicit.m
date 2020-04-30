function tests = testWIBL1SemiImplicit()
    tests = functiontests(localfunctions);
end

function testExplicitFunctionEvaluationSize(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = [icos(domain.x);icos(domain.x)];
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    actual = fwibl1Explicit(domain, h, params);
    
    verifySize(testCase, actual, [2*xN, yN]);
end

function testImplicitFunctionEvaluationSize(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = [icos(domain.x);icos(domain.x)];
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    actual = fwibl1Implicit(domain, h, params);
    
    verifySize(testCase, actual, [2*xN, yN]);
end

function testSplitCombinesToMakeFWIBL1(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 128; yN = 128;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = [icos(domain.x);2/3*icos(domain.x)];
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 1e-2); 
    expected = fwibl1(domain, h, params);
    actual = fwibl1Explicit(domain, h, params) + ...
        fwibl1Implicit(domain, h, params);
    
    figure;
    surf(expected);
    figure;
    surf(actual);
    figure;
    surf(log10(abs(actual - expected)));
    figure;
    surf(log10(abs((actual - expected)./expected)));

    max(abs(actual - expected), [], 'all')

    verifyEqual(testCase, actual, expected, 'AbsTol', 5e-2, 'RelTol', 2e-2);
end

function testSplitFunctionConverges(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi;
    coarseXN = 128; coarseYN = 128;
    coarseX = {linspace(xL/coarseXN,xL,coarseXN),linspace(yL/coarseYN,yL,coarseYN)};
    coarseDomain = FDDomain(coarseX, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    coarseH = [icos(coarseDomain.x);icos(coarseDomain.x)];
    fineXN = 256; fineYN = 256;
    fineX = {linspace(xL/fineXN,xL,fineXN),linspace(yL/fineYN,yL,fineYN)};
    fineDomain = FDDomain(fineX, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    fineH = [icos(fineDomain.x);icos(fineDomain.x)];
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    actual = fwibl1Explicit(coarseDomain, coarseH, params) + ...
        fwibl1Implicit(coarseDomain, coarseH, params);
    expected = fwibl1Explicit(fineDomain, fineH, params) + ...
        fwibl1Implicit(fineDomain, fineH, params);
    expected = [coarseDomain.interp(fineDomain.x, expected(1:end/2,:)); ...
        coarseDomain.interp(fineDomain.x, expected(1+end/2:end,:))];

    figure;
    surf(expected);
    figure;
    surf(actual);
    figure;
    surf(log10(abs(actual - expected)));
    figure;
    surf(log10(abs((actual - expected)./expected)));

    verifyEqual(testCase, actual, expected, 'AbsTol', 5e-2, 'RelTol', 1e-2);
end

function testSemiImplicitJacobian(testCase)
    addpath('discretisationMethods');

    xL = 2*pi; yL = 2*pi; xN = 48; yN = 48;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    fwibl1ImplicitVec = matFuncToVecFunc(@fwibl1Implicit);
    implicitOdefun = @(y) fwibl1ImplicitVec(domain, y, params);

    y = [icos(domain.x);icos(domain.x)];
    y0 = domain.reshapeToVector(y);

    expected = jacobianNumerical(implicitOdefun, y0);
    [~, actual] = fwibl1Implicit(domain, y, params);
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1.5e-1, 'AbsTol', 1e-7);

    function J = jacobianNumerical(F, u)
        N = length(u);
        epsilon = 1e-6;
        e = eye(N);
        J = sparse(N, N);

        Fu = F(u);

        for j = 1:N
            J(:, j) = (F(e(:, j)*epsilon+u) - Fu) / epsilon;
        end
    end
end

function testSemiImplicitEvolution(testCase)
    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    xL = 32; yL = 32; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);
    
    fwibl1ExplicitVec = matFuncToVecFunc(@fwibl1Explicit);
    explicitOdefun = @(t, y) fwibl1ExplicitVec(domain, y, params);

    t = linspace(0, 1,  10)';

    y0 = domain.reshapeToVector([icos(domain.x);icos(domain.x)]);

    timeStepper = @bdf1si;

    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true)));

    tic
    [~, y] = timeStepper(explicitOdefun, @(t, y) implicitOdefun(t, y, domain, params), t, y0, options);
    timeTaken = toc;
    y = y';
    y = reshape(y, [size(y,1), 1, size(y,2)]);
    y = domain.reshapeToDomain(y);

    verifySize(testCase, y, [2*xN, yN, length(t)]);

    function [F, J] = implicitOdefun(t, y, domain, params)
        [f, J] = fwibl1Implicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end

function testSemiImplicit1Dwibl1Wave(testCase)
    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    xL = 32; yL = 32; xN = 64; yN = 64;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);
    
    fwibl1ExplicitVec = matFuncToVecFunc(@fwibl1Explicit);
    explicitOdefun = @(t, y) fwibl1ExplicitVec(domain, y, params);

    t = linspace(0, 1,  10)';

    load('data/ic-WIBL1-2D-64.mat', 'h');
    y0 = domain.reshapeToVector(h);

    timeStepper = @bdf1si;

    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true)));

    tic
    [~, y] = timeStepper(explicitOdefun, @(t, y) implicitOdefun(t, y, domain, params), t, y0, options);
    timeTaken = toc;
    y = y';
    y = permute(y, [1, 3, 2]);
    y = domain.reshapeToDomain(y);

    size(y)
    save('temp.mat')

    verifySize(testCase, y, [2*xN, yN, length(t)]);

    function [F, J] = implicitOdefun(t, y, domain, params)
        [f, J] = fwibl1Implicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end

function testSemiImplicit2DWave(testCase)
    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    xL = 32; yL = 32; xN = 64; yN = 64;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);
    
    fwibl1ExplicitVec = matFuncToVecFunc(@fwibl1Explicit);
    explicitOdefun = @(t, y) fwibl1ExplicitVec(domain, y, params);

    t = linspace(0, 1,  100)';

    load('data/ic-WIBL1-3D-64.mat', 'h');
    y0 = domain.reshapeToVector(h);

    timeStepper = @bdf1si;

    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true)));

    tic
    [~, y] = timeStepper(explicitOdefun, @(t, y) implicitOdefun(t, y, domain, params), t, y0, options);
    timeTaken = toc;
    y = y';
    y = permute(y, [1, 3, 2]);
    y = domain.reshapeToDomain(y);

    save('temp.mat')

    verifySize(testCase, y, [2*xN, yN, length(t)]);

    function [F, J] = implicitOdefun(t, y, domain, params)
        [f, J] = fwibl1Implicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end

function testSemiImplicit2DWaveBDF2SI(testCase)
    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    xL = 32; yL = 32; xN = 64; yN = 64;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);
    
    fwibl1ExplicitVec = matFuncToVecFunc(@fwibl1Explicit);
    explicitOdefun = @(t, y) fwibl1ExplicitVec(domain, y, params);

    t = linspace(0, 1,  100)';

    load('data/ic-WIBL1-3D-64.mat', 'h');
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

    save('temp.mat')

    verifySize(testCase, y, [2*xN, yN, length(t)]);

    function [F, J] = implicitOdefun(t, y, domain, params)
        [f, J] = fwibl1Implicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end
