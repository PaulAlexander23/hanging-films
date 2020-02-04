function tests = testSetupInitialCondition()
    tests = functiontests(localfunctions);
end

function testCosBenney(testCase)
    model = "benney";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "finite-difference");
    domain = setupDomain(args);
    interface = @icos;
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = interface(domain.x);
    verifyEqual(testCase, actual, expected);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = domain.fft(interface(domain.x));
    verifyEqual(testCase, actual, expected);
end

function testCosWIBL1(testCase)
    model = "wibl1";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "finite-difference");
    domain = setupDomain(args);
    interface = @icosWIBL1;
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [icos(domain.x); ...
        0*domain.x{1} + 0*domain.x{2}+2/3];
    verifyEqual(testCase, actual, expected, 'AbsTol', 2*eps);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [domain.fft(icos(domain.x)); ...
        domain.fft(0*domain.x{1} + 0*domain.x{2}+2/3)];
    verifyEqual(testCase, actual, expected, 'AbsTol', 2*eps);
end

function testPertBenney(testCase)
    model = "benney";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "finite-difference");
    domain = setupDomain(args);
    interface = @ipert;
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = interface(domain.x);
    verifyEqual(testCase, actual, expected);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = domain.fft(interface(domain.x));
    verifyEqual(testCase, actual, expected);
end

function testPertWIBL1(testCase)
    model = "wibl1";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "finite-difference");
    domain = setupDomain(args);
    interface = @ipertWIBL1;
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [ipert(domain.x); ...
        -1/3 + ipert(domain.x)];
    verifyEqual(testCase, actual, expected, 'AbsTol', 2*eps);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [domain.fft(ipert(domain.x)); ...
        domain.fft(-1/3 + ipert(domain.x))];
    verifyEqual(testCase, actual, expected, 'AbsTol', 2*eps);
end

function testLoadBenney(testCase)
    model = "benney";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 256, 'yN', 256, 'method', "finite-difference");
    domain = setupDomain(args);
    interface = @(y) iload(y, "data/testInterfaceLoad.mat");
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = interface(domain.x);
    verifyEqual(testCase, actual, expected);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 256, 'yN', 256, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = domain.fft(interface(domain.x));
    verifyEqual(testCase, actual, expected);
end

function testLoadBenneyError(testCase)
    model = "benney";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "finite-difference");
    domain = setupDomain(args);
    interface = @(y) iload(y, "data/testInterfaceLoad.mat");
    method = "finite-difference";
    
    try
        setupInitialCondition(model, domain, interface, method);
    catch actualException
        verifyEqual(testCase, actualException.message, ...
            'Interface loaded has incorrect size. Expected size: [4, 8], Actual: [256, 256]')
    end
end

function testLoadWIBL1(testCase)
    model = "wibl1";
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 256, 'yN', 256, 'method', "finite-difference");
    domain = setupDomain(args);
    interface = @(y) iload(y, "data/testInterfaceLoadWIBL1.mat");
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [icos(domain.x); ...
        -1/3 + icos(domain.x)];
    verifyEqual(testCase, actual, expected);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 256, 'yN', 256, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [domain.fft(icos(domain.x)); ...
        domain.fft(-1/3 + icos(domain.x))];
    verifyEqual(testCase, actual, expected);
end