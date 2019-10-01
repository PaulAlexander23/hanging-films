function tests = testInterfaces()
    tests = functiontests(localfunctions);
end

function testInterfaceCos(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = icos(x);
    expected = 1 - 0.25 * cos(2*pi/x{1}(end)*x{1}) - ...
        0.25 * cos(2*pi/x{2}(end)*x{2}');
    verifyEqual(testCase, actual, expected);

    actual = icos(x, 0.5);
    expected = 1 - 0.5 * cos(2*pi/x{1}(end)*x{1}) - ...
        0.25 * cos(2*pi/x{2}(end)*x{2}');
    verifyEqual(testCase, actual, expected);

    actual = icos(x, 0.5, 0.5);
    expected = 1 - 0.5 * cos(2*pi/x{1}(end)*x{1}) - ...
        0.5 * cos(2*pi/x{2}(end)*x{2}');
    verifyEqual(testCase, actual, expected);
end

function testInterfaceSingle(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    A = 1e-2;

    actual = isingle(x, A);
    expected = 1 - A * (cos(2*pi/x{1}(end)*x{1}) + cos(2*pi/x{2}(end)*x{2}'));
    verifyEqual(testCase, actual, expected)
end

function testInterfaceRivulet(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    A = 2e-1;
    r = 0.05;

    actual = irivulet(x, A, r);
    expected = 1 + A * (-r * cos(2*pi/x{1}(end)*x{1}) - cos(2*pi/x{2}(end)*x{2}'));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceDoubleX(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = idoublex(x);
    expected = 1 + 0.2 * (-0.05 * cos(4*pi/x{1}(end)*x{1}) - cos(2*pi/x{2}(end)*x{2}'));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceDoubleY(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = idoubley(x);
    expected = 1 + 0.2 * (-0.05 * cos(2*pi/x{1}(end)*x{1}) - cos(4*pi/x{2}(end)*x{2}'));
    verifyEqual(testCase, actual, expected);
end

function testInterfacePert(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = ipert(x);
    expected = 1 + 1e-2 * exp(-1e-1*((x{1} - x{1}(end) / 2).^2 + (x{2}' - x{2}(end) / 2).^2));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceRand(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    A = 1e-4;
    numberOfModes = 5;

    actual = irand(x, A, numberOfModes);
    verifyTrue(testCase, max(max(abs(actual-1))) <= A);
    verifyEqual(testCase, sum(abs(fft2(actual-1)) > 1e-10, [1, 2]), ...
        (numberOfModes * 2)^2);
end

function testInterfaceLoad(testCase)
    x = setupX(1, 1, 2^8, 2^8);

    actual = iload(x, 'data/testInterfaceLoad');
    expected = icos(x);
    verifyEqual(testCase, actual, expected);
end

function testInterfaceLoadWithPert(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    alpha = 0.2;

    actual = iloadWithPert(x, 'data/testInterfaceLoad', alpha);
    expected = icos(x) + ipert(x, alpha) - 1;
    verifyEqual(testCase, actual, expected);
end

function testInterfaceLoadWithRand(testCase)
    x = setupX(1, 1, 2^8, 2^8);
    alpha = 0.2;

    actual = iloadWithRand(x, 'data/testInterfaceLoad', alpha);
    disturbance = actual - icos(x);
    verifyTrue(testCase, max(max(abs(disturbance))) <= alpha+1e-13);
    verifyEqual(testCase, sum(abs(fft2(disturbance)) > 1e-10, [1, 2]), ...
        (5 * 2)^2);
end