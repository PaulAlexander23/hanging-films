function tests = testMain()
    tests = functiontests(localfunctions);
end

function testMainWIBL1(testCase)
    model = "benney";
    theta = 1; Re = 1; C = 1;
    xLength = 2*pi; yLength = 2*pi; xN = 2^6; yN = 2^6;
    tFinal = 0.5;
    interface = @icos;
    method = "finite-difference";
    AbsTol = 1e-6;
    debug = false;
    
    filename = makeFilename("", ...
        struct('theta', theta, 'Re', Re, 'C', C), ...
        setupX(xLength, yLength, xN, yN), ...
        tFinal,interface,AbsTol,model);
    if isfile(filename)
        delete(filename)
    end
    
    main(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method, debug)
    
    load(filename, 'y')
    if isfile(filename)
        delete(filename)
    end
    
    actual = y(:,:,end);
    load('testCreate2DBenneyEquationExpected','expected');
    
    verifyEqual(testCase, actual, expected, ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end