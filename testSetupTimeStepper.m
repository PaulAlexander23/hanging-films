function tests = testSetupTimeStepper()
    tests = functiontests(localfunctions);
end

function testOptions(testCase)
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'AbsTol', 1e-5 ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        );
    timeStepperArguments = struct('timeStepper', @ode15s, 'odeopt', odeoptDefault);
    
    odeoptIVP = odeset();
    
    timeStepper = setupTimeStepper(timeStepperArguments, odeoptIVP);
    
    testODE = @(t, y) 1;
    testT = [0,1];
    testy0 = 1;
    
    testactual = timeStepper(testODE, testT, testy0);
    
    verifyEqual(testCase, testactual.solver, 'ode15s');
    
    actual = testactual.extdata.options;
    expected = odeset('AbsTol', 1e-5);
    
    verifyEqual(testCase, actual, expected);
end

function testDontOverwriteOptions(testCase)
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'AbsTol', 1e-5 ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        );
    timeStepperArguments = struct('timeStepper', @ode15s, 'odeopt', odeoptDefault);
    
    odeoptIVP = odeset('AbsTol', 1e-3);
    
    timeStepper = setupTimeStepper(timeStepperArguments, odeoptIVP);
    
    testODE = @(t, y) 1;
    testT = [0,1];
    testy0 = 1;
    
    testactual = timeStepper(testODE, testT, testy0);
    
    verifyEqual(testCase, testactual.solver, 'ode15s');
    
    actual = testactual.extdata.options;
    expected = odeset('AbsTol', 1e-5);
    
    verifyEqual(testCase, actual, expected);
end

function testSetExtraOptions(testCase)
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'AbsTol', 1e-5 ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        );
    timeStepperArguments = struct('timeStepper', @ode45, 'odeopt', odeoptDefault);
    
    odeoptIVP = odeset('MaxStep', 1e-3);
    
    timeStepper = setupTimeStepper(timeStepperArguments, odeoptIVP);
    
    testODE = @(t, y) 1;
    testT = [0,1];
    testy0 = 1;
    
    testactual = timeStepper(testODE, testT, testy0);
    
    verifyEqual(testCase, testactual.solver, 'ode45');
    
    actual = testactual.extdata.options;
    expected = odeset('AbsTol', 1e-5, 'MaxStep', 1e-3);
    
    verifyEqual(testCase, actual, expected);
end