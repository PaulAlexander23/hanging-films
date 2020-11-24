function tests = testDataHandling()
    tests = functiontests(localfunctions);
end

function testStruct2str(testCase)
    subStruct = struct('sub',1);
    testStruct = struct('a', 1, 'b', true, 'c', 's', 'd', "t", 'e', 0.1, ...
        'f', [], 'g', @icos, 'h', subStruct, 'i', @subfunction);
    
    actual = struct2str(testStruct);
    expected = "-a-1-b-1-c-s-d-t-e-0.1-g-icos-sub-1-i-subfunction";
    
    verifyEqual(testCase, actual, expected);
    
    function subfunction()
        
    end
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
        'Events', @eventNan,...
        'AbsTol', 1e-3, ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        'OutputFcn', 'odeprint',...
        'OutputSel', 1 ...
        );
    timeStepperArguments = struct('timeStepper', @ode15s, ...
        'odeopt', odeoptDefault);
    
    actual = makeFilename("test", ivpArguments, timePointsArguments);
    expected = "test-xLength-2-yLength-3-xN-16-yN-32-method-finite-difference-theta-1-Re-2-C-3-model-benney-tStep-0.2-tFinal-10";
    
    verifyEqual(testCase, actual, expected)
end

function testSaveFilenameLength(testCase)
    % 'benney',7*pi/8,0.2,0.01,32,32,100,@(x)iloadWithCos(x,"ic-benney-1d-Re-0_20.mat",0,0.1),64,64,1e-6,'finite-difference',@bdf2si,0.2
    AbsTol = 1e-6;
    timeStepper = @bdf2si;
    method = "finite-difference";
    interface = @(x)iloadWithCos(x,"ic-benney-1d-Re-0_20.mat",0,0.1);
    params = struct('theta', 7*pi/8, 'Re', 0.2, 'C', 0.01);
    
    domainArguments = struct('xLength', 32, 'yLength', 32, 'xN', 64, ...
        'yN', 64, 'method', method);
    
    explicitOdefun = @fbenney2dExplicit;
    implicitOdefun = @fbenney2dImplicit;
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    odejac = @jbenney2dImplicit;
    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'odefun',odefun,'odejac',odejac,'interface',interface);

    timePointsArguments = struct('tStep', 0.2, 'tFinal', 100);
    timerID = tic;
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'Events', @eventNan,...
        'AbsTol', AbsTol, ...
        'Events', @(~, ~) eTimeout(timerID, timeout) ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        ...'OutputFcn', 'odeprint',...
        ...'OutputSel', 1 ...
        );
    myoptimoptions = optimoptions('fsolve', ...
        'Display', 'off');
    odeoptDefault.optimoptions = myoptimoptions;

    timeStepperArguments = struct('timeStepper', timeStepper, ...
        'odeopt', odeoptDefault);
    
    expected = 1;
    
    saveData(expected, ivpArguments, timePointsArguments, timeStepperArguments);

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
        'AbsTol', 1e-3, ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        'OutputFcn', 'odeprint',...
        'OutputSel', 1 ...
        );
    timeStepperArguments = struct('timeStepper', @ode15s, ...
        'odeopt', odeoptDefault);
    
    expected = 1;
    
    saveData(expected, ivpArguments, timePointsArguments, timeStepperArguments);
    
    actual = loadData(ivpArguments, timePointsArguments);
    
    verifyEqual(testCase, actual, expected);
    
    filename = makeFilename("data", ivpArguments, timePointsArguments);
    delete(char(filename + '.mat'));
end
