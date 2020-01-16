function tests = testMain()
    tests = functiontests(localfunctions);
end

function testMainWIBL1(testCase)
    model = "benney";
    theta = 1;
    Re = 1;
    C = 1;
    xLength = 2 * pi;
    yLength = 2 * pi;
    xN = 2^6;
    yN = 2^6;
    tFinal = 0.5;
    interface = @icos;
    method = "finite-difference";
    AbsTol = 1e-6;

    filename = "data-xLength-6.2832-yLength-6.2832-xN-64-yN-64-method-finite-difference-theta-1-Re-1-C-1-model-benney-interface-icos-tStep-0.2-tFinal-0.5-timeStepper-ode15s-AbsTol-1e-06.mat";

    if isfile(filename)
        delete(filename)
    end

    main(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method)

    load(filename, 'solution')
    
    actual = solution.y(:, :, end);
    load('data/testCreate2DBenneyEquationExpected', 'expected');

    verifyEqual(testCase, actual, expected, ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
    if isfile(filename)
        delete(filename)
    end
end