function tests = testFunctionEvalution2()
    tests = functiontests(localfunctions);
end

function testFunctionOutputSize(testCase)
    addpath discretisationMethods
    N = 8;
    domain = FDDomain(setupX(1, 1, N, N), combvec(0:4,0:4), 4);
    y = 1 + 0.1 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}));
    F1 = 2 / 3 + 0.1 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}));
    F2 = 0 / 3 + 0.1 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}));
    Y = domain.reshapeToVector([y; F1; F2]);
    params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);

    actual = fwibl1(domain, Y, params);
    expectedSize = [3*N*N, 1];

    verifySize(testCase, actual, expectedSize)

    actual = fwibl2STF(domain, Y, params);
    expectedSize = [3*N*N, 1];

    verifySize(testCase, actual, expectedSize)

    actual = fwibl2(domain, Y, params);
    expectedSize = [3*N*N, 1];

    verifySize(testCase, actual, expectedSize)
end

