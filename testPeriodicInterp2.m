function tests = testPeriodicInterp2()
    tests = functiontests(localfunctions);
end

function testVectorLinear(testCase)
    X = linspace(1/20, 1, 20)';
    Y = linspace(1/20, 1, 20);
    V = sin(2*pi*X) + sin(2*pi*Y);
    Xq = linspace(1/40, 1, 40)';
    Yq = linspace(1/40, 1, 40);
    method = 'linear';

    expected = sin(2*pi*Xq) + sin(2*pi*Yq);
    actual = periodicInterp2(X, Y, V, Xq, Yq, method);

    verifyEqual(testCase, actual, expected, "AbsTol", 3e-2);
end

function testVectorSpline(testCase)
    X = linspace(1/20, 1, 20)';
    Y = linspace(1/20, 1, 20);
    V = sin(2*pi*X) + sin(2*pi*Y);
    Xq = linspace(1/40, 1, 40)';
    Yq = linspace(1/40, 1, 40);
    method = 'spline';

    expected = sin(2*pi*Xq) + sin(2*pi*Yq);
    actual = periodicInterp2(X, Y, V, Xq, Yq, method);

    verifyEqual(testCase, actual, expected, "AbsTol", 2e-4);
end

function testMatrix(testCase)
    X = linspace(1/20, 1, 20)';
    Y = linspace(1/20, 1, 20);
    [X, Y] = meshgrid(X, Y);
    V = sin(2*pi*X) + sin(2*pi*Y);
    Xq = linspace(1/40, 1, 40)';
    Yq = linspace(1/40, 1, 40);
    [Xq, Yq] = meshgrid(Xq, Yq);
    method = 'spline';

    expected = sin(2*pi*Xq) + sin(2*pi*Yq);
    actual = periodicInterp2(X, Y, V, Xq, Yq, method);

    verifyEqual(testCase, actual, expected, "AbsTol", 2e-4);
end
