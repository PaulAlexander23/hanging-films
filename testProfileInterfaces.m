function tests = testProfileInterfaces()
    tests = functiontests(localfunctions);
end

function testBenney(testCase)
    addpath("discretisationMethods");

    M = 64;
    N = 64;
    x = setupX(1, 1, M, N);

    actual = iloadProfileBenney(x, "data/ic-rivulet-benney.mat");

    verifySize(testCase, actual, [M,N]);
end

function testWIBL1STF(testCase)
    addpath("discretisationMethods");

    M = 64;
    N = 64;
    x = setupX(1, 1, M, N);

    actual = iloadProfileWIBL1STF(x, "data/ic-rivulet-wibl1-stf.mat");

    verifySize(testCase, actual, [2*M,N]);
end

function testWIBL1(testCase)
    addpath("discretisationMethods");

    M = 64;
    N = 64;
    x = setupX(1, 1, M, N);

    actual = iloadProfileWIBL1(x, "data/ic-rivulet-wibl1.mat");

    verifySize(testCase, actual, [3*M,N]);
end

function testWIBL2STF(testCase)
    addpath("discretisationMethods");

    M = 64;
    N = 64;
    x = setupX(1, 1, M, N);

    actual = iloadProfileWIBL2STF(x, "data/ic-rivulet-wibl2-stf.mat");

    verifySize(testCase, actual, [3*M,N]);
end

function testWIBL2(testCase)
    addpath("discretisationMethods");

    M = 64;
    N = 64;
    x = setupX(1, 1, M, N);

    actual = iloadProfileWIBL2(x, "data/ic-rivulet-wibl2.mat");

    verifySize(testCase, actual, [3*M,N]);
end
