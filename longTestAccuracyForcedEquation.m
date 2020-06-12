function tests = testAccuracyForcedEquation()
    tests = functiontests(localfunctions);
end

function testBenney(testCase)
    addpath('discretisationMethods');
    addpath('timeSteppingMethods');

    t = (0:0.001:0.01)';
    timeStep = 1e-3;
    xL = 2*pi; yL = 2*pi; xN = 164; yN = 164;
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 4);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01, ...
        'a', 0.5, 'b', 0.25, 'c', 0);
    h0 = 1 - params.a * cos(domain.x{1}) - params.b * cos(domain.x{2});

    y0 = domain.reshapeToVector(h0);

    explicitOdefun = @(t, y) fbenney2dExplicit(domain, y, params) ...
          + fbenneyforcing(t, domain, y0, params);
    implicitOdefun = @(t, y) fbenney2dImplicit(domain, y, params);
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    odejac = @(t, y) jbenney2dImplicit(domain, y, params);

    f0 = explicitOdefun(t(1), y0) + implicitOdefun(t(1), y0);
    norm(f0)

    surf(domain.reshapeToDomain(f0))
    
    timeStepper = @bdf2si;

    myoptimoptions = optimoptions('fsolve', ...
        'Display', 'iter', ...
        'SpecifyObjectiveGradient', true);
    options = odeset('MaxStep', timeStep);
    options.optimmethod = @fsolve;
    options.optimoptions = myoptimoptions;
    options.Jacobian = odejac;

    tic
    solution = timeStepper(odefun, t, y0, options);
    timeTaken = toc;

    t = reshape(t, [1, 1, length(t)]);
    expected = 1 - params.a * cos(domain.x{1} - params.c * t) ...
        - params.b * cos(domain.x{2});
    actual = solution.y;
    actual = reshape(actual', [size(actual,2), 1, size(actual,1)]);
    actual = domain.reshapeToDomain(actual);

    save("temp.mat")

    verifyEqual(testCase, actual, expected, 'AbsTol', 3e-2);

    fprintf("Time step: %.2e, ", timeStep)
    fprintf("Max diff: %.2e, ", max(abs(actual-expected),[],'all'))
    fprintf("Time taken: %.2f, ", timeTaken)
    fprintf("\n")
end
