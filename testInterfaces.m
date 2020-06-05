function tests = testInterfaces()
    tests = functiontests(localfunctions);
end

function testInterfaceCos(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = icos(x);
    expected = 1 - 0.25 * cos(2*pi/x{1}(end)*x{1}) - ...
        0.25 * cos(2*pi/x{2}(end)*x{2});
    verifyEqual(testCase, actual, expected);

    actual = icos(x, 0.5);
    expected = 1 - 0.5 * cos(2*pi/x{1}(end)*x{1}) - ...
        0.25 * cos(2*pi/x{2}(end)*x{2});
    verifyEqual(testCase, actual, expected);

    actual = icos(x, 0.5, 0.5);
    expected = 1 - 0.5 * cos(2*pi/x{1}(end)*x{1}) - ...
        0.5 * cos(2*pi/x{2}(end)*x{2});
    verifyEqual(testCase, actual, expected);
end

function testInterfaceSingle(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    A = 1e-2;

    actual = isingle(x, A);
    expected = 1 - A * (cos(2*pi/x{1}(end)*x{1}) + cos(2*pi/x{2}(end)*x{2}));
    verifyEqual(testCase, actual, expected)
end

function testInterfaceRivulet(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    A = 2e-1;
    r = 0.05;

    actual = irivulet(x, A, r);
    expected = 1 + A * (-r * cos(2*pi/x{1}(end)*x{1}) - cos(2*pi/x{2}(end)*x{2}));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceDoubleX(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = idoublex(x);
    expected = 1 + 0.2 * (-0.05 * cos(4*pi/x{1}(end)*x{1}) - cos(2*pi/x{2}(end)*x{2}));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceDoubleY(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = idoubley(x);
    expected = 1 + 0.2 * (-0.05 * cos(2*pi/x{1}(end)*x{1}) - cos(4*pi/x{2}(end)*x{2}));
    verifyEqual(testCase, actual, expected);
end

function testInterfacePert(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = ipert(x);
    expected = 1 + 1e-2 * exp(-1e-1*((x{1} - x{1}(end) / 2).^2 + (x{2} - x{2}(end) / 2).^2));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceRand(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    A = 1e-4;
    numberOfModes = 5;

    actual = irand(x, A, numberOfModes);
    verifyTrue(testCase, max(max(abs(actual-1))) <= A);
    verifyEqual(testCase, sum(sum(abs(fft2(actual-1)) > 1e-10)), ...
        (numberOfModes * 2)^2);
end

function testInterfaceRandIsNotGridDependant(testCase)
    x = setupX(1, 1, 64, 64);
    x2 = setupX(1, 1, 75, 75);

    A = 1e-4;
    numberOfModes = 5;
    seed = 1;

    expected = irand(x, A, numberOfModes, seed);
    [X, Y] = meshgrid(x{1}, x{2});
    [X2, Y2] = meshgrid(x2{1}, x2{2});
    actual = periodicInterp2(X2, Y2, ...
        irand(x2, A, numberOfModes, seed), ...
        X, Y, 'spline');

    verifyEqual(testCase, actual, expected, 'AbsTol', 2e-8)
end

function testInterfaceRandLinIsNotGridDependant(testCase)
    x = setupX(1, 1, 64, 64);
    x2 = setupX(1, 1, 75, 75);

    A = 1e-4;
    numberOfModes = 5;
    seed = 1;

    expected = irandLin(x, A, numberOfModes, seed);
    [X, Y] = meshgrid(x{1}, x{2});
    [X2, Y2] = meshgrid(x2{1}, x2{2});
    actual = periodicInterp2(X2, Y2, ...
        irandLin(x2, A, numberOfModes, seed), ...
        X, Y, 'spline');

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-8)
end

function testInterfaceLoad(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    
    actual = iload(x, 'data/testInterfaceLoad.mat');
    expected = icos(x);
    verifyEqual(testCase, actual, expected);
end

function testInterfaceLoadWithPert(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    alpha = 0.2;

    actual = iloadWithPert(x, 'data/testInterfaceLoad.mat', alpha);
    expected = icos(x) + ipert(x, alpha) - 1;
    verifyEqual(testCase, actual, expected);
end

function testInterfaceLoadWithRand(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    alpha = 0.2;

    actual = iloadWithRand(x, 'data/testInterfaceLoad.mat', alpha);
    disturbance = actual - icos(x);

    verifyTrue(testCase, max(max(abs(disturbance))) <= alpha+1e-2);
    verifyEqual(testCase, sum(sum(abs(fft2(disturbance)) > 1e-10)), ...
        (5 * 2)^2);
end

function testInterfaceLoadWithCos(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    a = 0.2;
    b = 0.2;

    actual = iloadWithCos(x, 'data/testInterfaceLoad.mat', a, b);
    expected = icos(x) + icos(x, a, b) - 1;
    verifyEqual(testCase, actual, expected);
end

function testInterfaceLoadWithCosWIBL1(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    a = 0.1;
    b = 0.2;
    c = 0.3;
    d = 0.4;

    actual = iloadWithCosWIBL1(x, 'data/testInterfaceLoadWIBL1.mat', ...
        a, b, c, d);
    expected = [icos(x, 0.25, 0.25) + icos(x, a, b) - 1; ...
        icos(x, 0.25, 0.25) - 1/3 + icos(x, c, d) - 1];

    verifyEqual(testCase, actual, expected, 'AbsTol', eps);
end

function testInterfaceLoadAmplifyMode(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    amplitude = 1;
    xMode = 1;
    yMode = 0;
    
    expected = icos(x, 0.5, 0.25);
    actual = iloadAmplifyMode(x, 'data/testInterfaceLoad.mat', amplitude, xMode, yMode);

    verifyEqual(testCase, actual, expected, 'AbsTol', 4*eps)
end

function testInterfaceLoadAmplifyModeWIBL1(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    amplitude = 1;
    xMode = 1;
    yMode = 0;
    
    expected = [icos(x, 0.5, 0.25);icos(x,0.5,0.25) - 1/3];
    actual = iloadAmplifyModeWIBL1(x, 'data/testInterfaceLoadWIBL1.mat', amplitude, xMode, yMode);

    verifyEqual(testCase, actual, expected, 'AbsTol', 4*eps)
end
