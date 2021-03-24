function tests = testProfileInterfaces()
    tests = functiontests(localfunctions);
end

function testBenney(testCase)
    addpath("discretisationMethods");

    M = 72;
    N = 60;
    x = setupX(40, 32, M, N);
    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    actual = iloadProfileEigBenney(x, "data/ic-rivulet-benney.mat", 1e-3, params);

    verifySize(testCase, actual, [M,N]);
end

function testWIBL1STF(testCase)
    addpath("discretisationMethods");

    M = 72;
    N = 60;
    x = setupX(40, 32, M, N);
    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    actual = iloadProfileEigWIBL1STF(x, "data/ic-rivulet-wibl1-stf.mat", 1e-3, params);

    verifySize(testCase, actual, [2*M,N]);
end

function testWIBL1(testCase)
    addpath("discretisationMethods");

    M = 72;
    N = 60;
    x = setupX(40, 32, M, N);
    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    actual = iloadProfileEigWIBL1(x, "data/ic-rivulet-wibl1.mat", 1e-3, params);

    verifySize(testCase, actual, [3*M,N]);
end

function testWIBL2STF(testCase)
    addpath("discretisationMethods");

    M = 72;
    N = 60;
    x = setupX(50.9, 32, M, N);
    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    actual = iloadProfileEigWIBL2STF(x, "data/ic-rivulet-wibl2-stf.mat", 1e-3, params);

    verifySize(testCase, actual, [3*M,N]);
end

function testWIBL2(testCase)
    addpath("discretisationMethods");

    M = 72;
    N = 60;
    x = setupX(40, 32, M, N);
    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    actual = iloadProfileEigWIBL2(x, "data/ic-rivulet-wibl2.mat", 1e-3, params);

    verifySize(testCase, actual, [3*M,N]);
end

