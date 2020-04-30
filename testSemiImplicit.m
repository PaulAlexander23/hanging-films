function tests = testSemiImplicit()
    tests = functiontests(localfunctions);
end

function testLaplacianEqual(testCase)
    addpath('discretisationMethods');
    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = icos(domain.x);
    z = h.^4 / 4;

    expected = domain.diff(h, [2; 0]) + domain.diff(h, [0; 2]);
    actual = h.^(-3) .* ( ...
        domain.diff(z, [2; 0]) + domain.diff(z, [0; 2]) - ...
        3/4 * (domain.diff(z, [1; 0]).^2 + domain.diff(z, [0; 1]).^2)./z);

    % figure; surf(expected);
    % figure; surf(actual);
    % figure; surf(actual - expected);

    verifyEqual(testCase, actual, expected, 'RelTol', 3e-2, 'AbsTol', 5e-3)
end

function testExplicitFunctionEvaluationSize(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = icos(domain.x);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    actual = fbenney2dExplicit(domain, h, params);
    
    verifySize(testCase, actual, [xN, yN]);
end

function testImplicitFunctionEvaluationSize(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = icos(domain.x);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    actual = fbenney2dImplicit(domain, h, params);
    
    verifySize(testCase, actual, [xN, yN]);
end

function testSplitCombinesToMakeFBenney2d(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 128; yN = 128;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = icos(domain.x);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    expected = fbenney2d(domain, h, params);
    actual = fbenney2dExplicit(domain, h, params) + ...
        fbenney2dImplicit(domain, h, params);
    
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
    coarseH = icos(coarseDomain.x);
    fineXN = 256; fineYN = 256;
    fineX = {linspace(xL/fineXN,xL,fineXN),linspace(yL/fineYN,yL,fineYN)};
    fineDomain = FDDomain(fineX, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    fineH = icos(fineDomain.x);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    actual = fbenney2dExplicit(coarseDomain, coarseH, params) + ...
        fbenney2dImplicit(coarseDomain, coarseH, params);
    expected = fbenney2dExplicit(fineDomain, fineH, params) + ...
        fbenney2dImplicit(fineDomain, fineH, params);
    expected = coarseDomain.interp(fineDomain.x, expected);

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

    fbenney2dImplicitVec = matFuncToVecFunc(@fbenney2dImplicit);
    implicitOdefun = @(y) fbenney2dImplicitVec(domain, y, params);

    y = icos(domain.x);
    y0 = domain.reshapeToVector(y);

    expected = jacobianNumerical(implicitOdefun, y0);
    [~, actual] = fbenney2dImplicit(domain, y, params);
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1.5e-1);

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
    
    fbenney2dExplicitVec = matFuncToVecFunc(@fbenney2dExplicit);
    explicitOdefun = @(t, y) fbenney2dExplicitVec(domain, y, params);

    t = linspace(0, 1,  100)';

    y0 = domain.reshapeToVector(icos(domain.x));

    timeStepper = @bdf1si;

    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true)));

    tic
    [~, y] = timeStepper(explicitOdefun, @(t, y) implicitOdefun(t, y, domain, params), t, y0, options);
    timeTaken = toc;
    y = y';
    y = reshape(y, [size(y,1), 1, size(y,2)]);
    y = domain.reshapeToDomain(y);

    size(y)
    save('temp.mat')

    verifySize(testCase, y, [xN, yN, length(t)]);

    function [F, J] = implicitOdefun(t, y, domain, params)
        [f, J] = fbenney2dImplicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end

function testSemiImplicit1DbenneyWave(testCase)
    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    xL = 32; yL = 32; xN = 64; yN = 64;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);
    
    fbenney2dExplicitVec = matFuncToVecFunc(@fbenney2dExplicit);
    explicitOdefun = @(t, y) fbenney2dExplicitVec(domain, y, params);

    t = linspace(0, 40,  1000)';

    load('data/ic-Benney-2D-64.mat', 'h');
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

    verifySize(testCase, y, [xN, yN, length(t)]);

    function [F, J] = implicitOdefun(t, y, domain, params)
        [f, J] = fbenney2dImplicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end

function testSemiImplicit2DbenneyWave(testCase)
    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    xL = 32; yL = 32; xN = 64; yN = 64;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);
    
    fbenney2dExplicitVec = matFuncToVecFunc(@fbenney2dExplicit);
    explicitOdefun = @(t, y) fbenney2dExplicitVec(domain, y, params);

    t = linspace(0, 20,  200)';

    load('data/ic-Benney-3D-64.mat', 'h');
    %load('tempic.mat', 'h');
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

    verifySize(testCase, y, [xN, yN, length(t)]);

    function [F, J] = implicitOdefun(t, y, domain, params)
        [f, J] = fbenney2dImplicit(domain, domain.reshapeToDomain(y),params);
        
        F = domain.reshapeToVector(f);
    end
end
