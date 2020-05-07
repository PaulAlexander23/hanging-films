function tests = testSetupODEFunction()
    tests = functiontests(localfunctions);
end

function test2DBenneyFiniteDifference(testCase)
    model = "benney";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "finite-difference");
    domain = setupDomain(args);
    params = struct('theta', 1, 'Re', 2, 'C', 3);
    
    actual = setupODEFunction(model, domain, params);
    expected = @(t, y) fbenney2d(domain, y, params);
    
    testT = 1;
    testY = 1 + 0.1 * (cos(domain.x{1}) + cos(domain.x{2}));

    verifyEqual(testCase, ...
        actual(testT, domain.reshapeToVector(testY)), ...
        expected(testT, domain.reshapeToVector(testY)));
end

function test2DBenneyPseudoSpectral(testCase)
    model = "benney";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    params = struct('theta', 1, 'Re', 2, 'C', 3);
    
    actual = setupODEFunction(model, domain, params);
    expected = @(t, y) fbenney2d(domain, y, params);
    
    testT = 1;
    testY = domain.fft(1 + 0.1 * (cos(domain.x{1}) + cos(domain.x{2})));

    verifyEqual(testCase, ...
        actual(testT, domain.reshapeToVector(testY)), ...
        expected(testT, domain.reshapeToVector(testY)));
end

function testWIBL1FiniteDifference(testCase)
    model = "wibl1";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "finite-difference");
    domain = setupDomain(args);
    params = struct('theta', 1, 'Re', 2, 'C', 3);
    
    actual = setupODEFunction(model, domain, params);
    expected = @(t, y) fwibl1(domain, y, params);
    
    testT = 1;
    testY = repmat(1 + 0.1 * (cos(domain.x{1}) + cos(domain.x{2})),2,1);
    
    verifyEqual(testCase, ...
        actual(testT, domain.reshapeToVector(testY)), ...
        expected(testT, domain.reshapeToVector(testY)));
end

function testWIBL1PseudoSpectral(testCase)
    model = "wibl1";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    params = struct('theta', 1, 'Re', 2, 'C', 3);
    
    actual = setupODEFunction(model, domain, params);
    expected = @(t, y) fwibl1(domain, y, params);
    
    testT = 1;
    testY = repmat(domain.fft(1 + 0.1 * (cos(domain.x{1}) + cos(domain.x{2}))),2,1);
    
    verifyEqual(testCase, ...
        actual(testT, domain.reshapeToVector(testY)), ...
        expected(testT, domain.reshapeToVector(testY)));
end
