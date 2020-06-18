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

    h = domain.reshapeToVector(h);
    actual = fwibl1Explicit(domain, h, params);
    
    verifySize(testCase, actual, [2*xN * yN, 1]);
end

function testImplicitFunctionEvaluationSize(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = [icos(domain.x);icos(domain.x)];
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    h = domain.reshapeToVector(h);
    actual = fwibl1Implicit(domain, h, params);
    
    verifySize(testCase, actual, [2*xN * yN, 1]);
end

function testSplitCombinesToMakeFWIBL1(testCase)
    addpath('discretisationMethods')
    xL = 2*pi; yL = 2*pi; xN = 128; yN = 128;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    h = [icos(domain.x);2/3*icos(domain.x)];
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 1e-2); 
    h = domain.reshapeToVector(h);
    expected = fwibl1(domain, h, params);
    actual = fwibl1Explicit(domain, h, params) + ...
        fwibl1Implicit(domain, h, params);
    
    % figure;
    % surf(expected);
    % figure;
    % surf(actual);
    % figure;
    % surf(log10(abs(actual - expected)));
    % figure;
    % surf(log10(abs((actual - expected)./expected)));

    % max(abs(actual - expected), [], 'all')

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

    coarseH = coarseDomain.reshapeToVector(coarseH);
    fineH = fineDomain.reshapeToVector(fineH);
    actual = fwibl1Explicit(coarseDomain, coarseH, params) + ...
        fwibl1Implicit(coarseDomain, coarseH, params);
    actual = coarseDomain.reshapeToDomain(actual);
    expected = fwibl1Explicit(fineDomain, fineH, params) + ...
        fwibl1Implicit(fineDomain, fineH, params);
    expected = fineDomain.reshapeToDomain(expected);
    expected = [coarseDomain.interp(fineDomain.x, expected(1:end/2,:)); ...
        coarseDomain.interp(fineDomain.x, expected(1+end/2:end,:))];

    % figure;
    % surf(expected);
    % figure;
    % surf(actual);
    % figure;
    % surf(log10(abs(actual - expected)));
    % figure;
    % surf(log10(abs((actual - expected)./expected)));

    verifyEqual(testCase, actual, expected, 'AbsTol', 5e-2, 'RelTol', 1e-2);
end

function testSemiImplicitJacobian(testCase)
    addpath('discretisationMethods');

    xL = 2*pi; yL = 2*pi; xN = 32; yN = 32;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 2);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    implicitOdefun = @(y) fwibl1Implicit(domain, y, params);

    y = [icos(domain.x);icos(domain.x)];
    y0 = domain.reshapeToVector(y);

    expected = jacobianNumerical(implicitOdefun, y0);
    actual = jwibl1Implicit(domain, y0, params);
    
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
    
    explicitOdefun = @(t, y) fwibl1Explicit(domain, y, params);
    implicitOdefun = @(t, y) fwibl1Implicit(domain, y, params);
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    odejac = @(t, y) jwibl1Implicit(domain, y, params);

    t = linspace(0, 1,  10)';

    y0 = domain.reshapeToVector([icos(domain.x);icos(domain.x)]);

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
    timeTaken = toc;
    y = y';
    y = reshape(y, [size(y,1), 1, size(y,2)]);
    y = domain.reshapeToDomain(y);

    verifySize(testCase, y, [2*xN, yN, length(t)]);
end
