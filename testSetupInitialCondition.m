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
    interface = @icos;
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [interface(domain.x); ...
        2/3 + 0*interface(domain.x)];
    verifyEqual(testCase, actual, expected);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [domain.fft(interface(domain.x)); ...
        domain.fft(2/3 + 0 * interface(domain.x))];
    verifyEqual(testCase, actual, expected);
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
    interface = @ipert;
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [interface(domain.x); ...
        2/3 + 0*interface(domain.x)];
    verifyEqual(testCase, actual, expected);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [domain.fft(interface(domain.x)); ...
        domain.fft(2/3 + 0 * interface(domain.x))];
    verifyEqual(testCase, actual, expected);
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
    interface = @(y) iload(y, "data/testInterfaceLoad.mat");
    method = "finite-difference";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [interface(domain.x); ...
        2/3 + 0*interface(domain.x)];
    verifyEqual(testCase, actual, expected);
    
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 256, 'yN', 256, 'method', "pseudo-spectral");
    domain = setupDomain(args);
    method = "pseudo-spectral";
    
    actual = setupInitialCondition(model, domain, interface, method);
    expected = [domain.fft(interface(domain.x)); ...
        domain.fft(2/3 + 0 * interface(domain.x))];
    verifyEqual(testCase, actual, expected);
end