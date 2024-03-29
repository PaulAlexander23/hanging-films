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

    verifyEqual(testCase, actual, expected, 'RelTol', 3e-2, 'AbsTol', 5e-3)
end

function testExplicitFunctionEvaluationSize(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = icos(domain.x);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    h = domain.reshapeToVector(h);
    actual = fbenney2dExplicit(domain, h, params);
    
    verifySize(testCase, actual, [xN * yN, 1]);
end

function testImplicitFunctionEvaluationSize(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = icos(domain.x);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    h = domain.reshapeToVector(h);
    actual = fbenney2dImplicit(domain, h, params);
    
    verifySize(testCase, actual, [xN * yN, 1]);
end

function testSplitCombinesToMakeFBenney2d(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 128; yN = 128;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = icos(domain.x);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    h = domain.reshapeToVector(h);
    expected = fbenney2d(domain, h, params);
    actual = fbenney2dExplicit(domain, h, params) + ...
        fbenney2dImplicit(domain, h, params);
    
    verifyEqual(testCase, actual, expected, 'AbsTol', 4e-1, 'RelTol', 4e-2);
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

    coarseH = coarseDomain.reshapeToVector(coarseH);
    fineH = fineDomain.reshapeToVector(fineH);
    actual = fbenney2dExplicit(coarseDomain, coarseH, params) + ...
        fbenney2dImplicit(coarseDomain, coarseH, params);
    actual = coarseDomain.reshapeToDomain(actual);
    expected = fbenney2dExplicit(fineDomain, fineH, params) + ...
        fbenney2dImplicit(fineDomain, fineH, params);
    expected = fineDomain.reshapeToDomain(expected);
    expected = coarseDomain.interp(fineDomain.x, expected);

    verifyEqual(testCase, actual, expected, 'AbsTol', 5e-1, 'RelTol', 2e-2);
end

function testSemiImplicitJacobian(testCase)
    addpath('discretisationMethods');

    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    implicitOdefun = @(y) fbenney2dImplicit(domain, y, params);

    y = icos(domain.x);
    y0 = domain.reshapeToVector(y);

    expected = jacobianNumerical(implicitOdefun, y0);
    actual = jbenney2dImplicit(domain, y0, params);
    
    verifyEqual(testCase, actual, expected, 'RelTol', 2.5e-1);

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
    
    explicitOdefun = @(t, y) fbenney2dExplicit(domain, y, params);
    implicitOdefun = @(t, y) fbenney2dImplicit(domain, y, params);
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    odejac = @(t, y) jbenney2dImplicit(domain, y, params);

    t = linspace(0, 0.1,  10)';

    y0 = domain.reshapeToVector(icos(domain.x));

    timeStepper = @bdf1si;

    myoptimoptions = optimoptions('fsolve', ...
        'Display', 'off', ...
        'SpecifyObjectiveGradient', true);
    options = odeset();
    options.optimmethod = @fsolve;
    options.optimoptions = myoptimoptions;
    options.Jacobian = odejac;

    tic
    [~, y] = timeStepper(odefun, t, y0, options);
    y = y';
    y = reshape(y, [size(y,1), 1, size(y,2)]);
    y = domain.reshapeToDomain(y);

    verifySize(testCase, y, [xN, yN, length(t)]);
end
