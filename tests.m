function tests = tests()
    tests = functiontests(localfunctions);
end

function testSetupX(testCase)
    actualX = setupX(10,pi,12,11);
    verifyEqual(testCase, actualX{1}, linspace(10/12,10,12)')
    verifyEqual(testCase, actualX{2}, linspace(pi/11,pi,11)')
end

function testMakeFilename(testCase)
    actual = makeFilename('-test',[1,2,3,4],{0.1:0.1:3.3,1:1:10},pi,@interface,1e-6);
    expected = "data-test-theta-2-Re-3-C-4-xL-3_3-yL-10-T-3_14159-interface-interface-xN-33-yN-10-AbsTol-1e-06";
    verifyEqual(testCase, actual, expected)
end

function testSaveAndLoadData(testCase)
    filename = 'data-test-theta-0_876255-Re-0_488629-C-0_407077-xL-3_3-yL-10-T-3_14159-interface-interface-xN-33-yN-10-AbsTol-1e-06.mat';
    if isfile(filename)
        delete(filename)
    end
    
    expectedY = rand(100);
    expectedParams = [0.9618, 0.876255, 0.488629, 0.407077];
    expectedT = linspace(0,0.126576313737181);
    expectedX = {0.1:0.1:3.3,1:1:10};
    expectedTimeTaken = 0.925425280986515;
    saveData(expectedY, expectedParams, expectedT, expectedX, ...
        expectedTimeTaken, pi, @interface, 1e-6, "-test");
    
    verifyTrue(testCase,isfile(filename))
    
    [actualY, actualParams, actualT, actualX, actualTimeTaken] = loadData(expectedParams,expectedX,pi,@interface,1e-6,"-test");
    
    verifyEqual(testCase,actualY,expectedY)
    verifyEqual(testCase,actualParams,expectedParams)
    verifyEqual(testCase,actualT,expectedT)
    verifyEqual(testCase,actualX,expectedX)
    verifyEqual(testCase,actualTimeTaken,expectedTimeTaken)
    
    if isfile(filename)
        delete(filename)
    end
end

function testPDESolver1DExponentialDecay(testCase)
    pdefun = @(t, x, y, diffMethod) - y;
    t = linspace(0, 1, 11);
    x = linspace(2*pi/64, 2*pi, 64)';
    y0 = cos(x);
    diffMethod = @(x, y, degree) 0 * y;
    timestepper = @ode45;
    actual = pdeSolver(pdefun, t, x, y0, diffMethod, timestepper);
    expected = y0 * exp(-t);
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolver1DTravellingWaveEquation(testCase)
    problemDiffDegrees = 1;
    function out = pdefun(t, x, y, diffMethod)
        out = diffMethod(x, y, 1);
        out = -out{1};
    end
    t = linspace(0, 1, 11);
    x = {linspace(2*pi/64, 2*pi, 64)'};
    y0 = cos(x{1});
    import discretisationMethods.*
    diffMat = init_fd(x, problemDiffDegrees, 4);
    diffMethod = @(x, y, degree) diff_fd(x, y, degree, diffMat, ...
        problemDiffDegrees);
    timestepper = @ode45;
    actual = pdeSolver(@pdefun, t, x, y0, diffMethod, timestepper);
    expected = cos(x{1}-t);
    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolver1DHeatEquation(testCase)
    function out = pdefun(t, x, y, diffMethod)
        out = diffMethod(x, y, 2);
        out = out{1};
    end
    t = linspace(0, 1, 11);
    x = {linspace(2*pi/64, 2*pi, 64)'};
    y0 = cos(x{1});
    import discretisationMethods.*
    diffMat = init_fd(x, 2, 4);
    diffMethod = @(x, y, degree) diff_fd(x, y, degree, diffMat, 2);
    timestepper = @ode45;
    actual = pdeSolver(@pdefun, t, x, y0, diffMethod, timestepper);
    expected = y0 * exp(-t);
    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolver1DHeatEquationPseudoSpectral(testCase)
    function out = pdefun(t, x, y, diffMethod)
        out = diffMethod(x, y, 2);
    end
    t = linspace(0, 1, 11);
    x = linspace(2*pi/64, 2*pi, 64)';
    y0 = cos(x);
    import discretisationMethods.*
    diffMethod = @(x, y, degree) diff_ps(x, y, degree);
    timestepper = @ode45;
    actual = pdeSolver(@pdefun, t, x, y0, diffMethod, timestepper);
    expected = y0 * exp(-t);
    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolver1DBenneyEquation(testCase)
    t = linspace(0, 1, 11);
    x = {linspace(2*pi/64, 2*pi, 64)'};
    y0 = 1+0.5*cos(x{1});
    problemDiffDegrees = [1, 2];
    import discretisationMethods.*
    diffMat = init_fd(x, problemDiffDegrees, 4);
    diffMethod = @(x, y, degree) diff_fd(x, y, degree, diffMat, problemDiffDegrees);
    params = [1,1,1,1];
    pdefun = @(t, x, y, diffMethod) fbenney1d(x, y, params, diffMethod);
    timestepper = @ode15s;
    actual = pdeSolver(pdefun, t, x, y0, diffMethod, timestepper);
    load('testPDESolver1DBenneyEquationExpected','expected')
    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testCreate2DBenneyEquation(testCase)
    filename = 'data-theta-1-Re-1-C-1-xL-6_28319-yL-6_28319-T-0_5-interface-@(x)1+0_5*cos(x{1}+x{2}'')-xN-64-yN-64-AbsTol-1e-06.mat';
    if isfile(filename)
        delete(filename)
    end
    create(1,1,1,2*pi,2*pi,0.5,@(x)1+0.5*cos(x{1}+x{2}'))
    load(filename,'y');
    actual = y;
    load('testCreate2DBenneyEquationExpected','expected')
    verifyEqual(testCase, actual(:, :, end), expected(:, :, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
    if isfile(filename)
        delete(filename)
    end
end

function testCreateWIBL1GetsToFinalTime(testCase)
    filename = 'data-wibl1-theta-2_74889-Re-1-C-0_01-xL-32-yL-32-T-1-interface-icos-xN-32-yN-32-AbsTol-1e-06.mat';
    if isfile(filename)
        delete(filename)
    end
    createWIBL1(7/8 * pi, 1, 0.01, 32, 32, 1, @icos, 32, 32, 1e-6)
    load(filename,'t')
    
    verifyEqual(testCase, t(end), 1)
    if isfile(filename)
        delete(filename)
    end
end