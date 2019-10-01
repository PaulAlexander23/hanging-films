function tests = testSetup()
    tests = functiontests(localfunctions);
end

function testSetupX(testCase)
    actual = setupX(10, pi, 12, 11);
    verifyEqual(testCase, actual{1}, linspace(10/12, 10, 12)')
    verifyEqual(testCase, actual{2}, linspace(pi/11, pi, 11)')
end

function testSetupT(testCase)
    actual = setupT(10, 0.2);
    verifyEqual(testCase, actual(end), 10);
    verifyTrue(testCase, all(abs(diff(actual(1:end-1))-0.2) < 1e-14));
end

function testParamsToStruct(testCase)
    actual = paramsToStruct(1, 1, 1);
    expected = struct('theta', 1, 'Re', 1, 'C', 1);
    verifyEqual(testCase, actual, expected);
end

function testCreateDomainFiniteDifference(testCase)
    xLength = 32;
    yLength = 32;
    xN = 32;
    yN = 32;
    method = "finite-difference";
    actual = createDomain(xLength, yLength, xN, yN, method);
    expected = FDDomain(setupX(xLength, yLength, xN, yN), ...
        [1, 0; 0, 1; 2, 0; 0, 2]', 4);

    verifyEqual(testCase, actual, expected)
end

function testCreateDomainPseudoSpectral(testCase)
    xLength = 32;
    yLength = 32;
    xN = 32;
    yN = 32;
    method = "pseudo-spectral";
    actual = createDomain(xLength, yLength, xN, yN, method);
    expected = PSDomain(setupX(xLength, yLength, xN, yN), true, false);

    verifyEqual(testCase, actual, expected)
end