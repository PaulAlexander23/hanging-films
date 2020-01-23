function tests = testDataHandling()
    tests = functiontests(localfunctions);
end

function testStruct2str(testCase)
    subStruct = struct('sub',1);
    testStruct = struct('a', 1, 'b', true, 'c', 's', 'd', "t", 'e', 0.1, ...
        'f', [], 'g', @icos, 'h', subStruct);
    
    actual = struct2str(testStruct);
    expected = "-a-1-b-1-c-s-d-t-e-0.1-g-icos-sub-1";
    
    verifyEqual(testCase, actual, expected);
end

function testEnsureUnique(testCase)
    filename = "testTemp";
    save(filename)
    
    actual = ensureUnique(filename);
    expected = "testTemp-1";
    
    verifyEqual(testCase, actual, expected);
    
    delete(char(filename + ".mat"));
end

function testMakeFilename(testCase)
    
    method = "finite-difference";
    params = struct('theta', 1, 'Re', 2, 'C', 3);
    
    domainArguments = struct('xLength', 2, 'yLength', 3, 'xN', 16, ...
        'yN', 32, 'method', method);
    
    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'model', "benney", 'interface', @icos);
    timePointsArguments = struct('tStep', 0.2, 'tFinal', 10);
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'AbsTol', 1e-3 ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        );
    timeStepperArguments = struct('timeStepper', @ode15s, 'odeopt', odeoptDefault);
    
    actual = makeFilename("test", ivpArguments, timePointsArguments, timeStepperArguments);
    expected = "test-xLength-2-yLength-3-xN-16-yN-32-method-finite-difference-theta-1-Re-2-C-3-model-benney-interface-icos-tStep-0.2-tFinal-10";
    
    verifyEqual(testCase, actual, expected)
end

function testSaveAndLoadData(testCase)
    method = "finite-difference";
    params = struct('theta', 1, 'Re', 2, 'C', 3);
    
    domainArguments = struct('xLength', 2, 'yLength', 3, 'xN', 16, ...
        'yN', 32, 'method', method);
    
    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'model', "benney", 'interface', @icos);
    timePointsArguments = struct('tStep', 0.2, 'tFinal', 10);
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'AbsTol', 1e-3 ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        );
    timeStepperArguments = struct('timeStepper', @ode15s, 'odeopt', odeoptDefault);
    
    expected = 1;
    
    saveData(expected, ivpArguments, timePointsArguments, timeStepperArguments);
    
    actual = loadData(ivpArguments, timePointsArguments, timeStepperArguments);
    
    verifyEqual(testCase, actual, expected);
    
    filename = makeFilename("data", ivpArguments, timePointsArguments, timeStepperArguments);
    delete(char(filename + '.mat'));
end
