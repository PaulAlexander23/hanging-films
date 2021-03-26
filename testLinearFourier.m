function tests = testLinearFourier()
    tests = functiontests(localfunctions);
end

function testBenneyOutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-benney.mat");
    hbar = data.y;

    alpha = 1;

    actual = linearisedBenneyFourier(domain, hbar, params, alpha);

    verifySize(testCase, actual, [zN/2,zN/2]);

end

function testBenneyDROutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-benney.mat");
    hbar = data.y;

    alpha = linspace(1e-3,1)';
    modes = zN/2;

    actual = dispersionRelationBenneyFourier(domain, hbar, params, alpha, modes)

    verifySize(testCase, actual, [100,zN/2]);

end

function testWIBL1STFOutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-wibl1-stf.mat");
    hbar = data.y(1:end/2);
    F1bar = data.y(1+end/2:end);

    alpha = 1;

    actual = linearisedWIBL1STFFourier(domain, hbar, F1bar, params, alpha);

    verifySize(testCase, actual, 2*[zN/2,zN/2]);

end

function testWIBL1STFDROutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-wibl1-stf.mat");
    hbar = data.y(1:end/2);
    F1bar = data.y(1+end/2:end);

    alpha = linspace(1e-3,1)';
    modes = zN/2;

    actual = dispersionRelationWIBL1STFFourier(domain, hbar, F1bar, params, alpha, modes)

    verifySize(testCase, actual, [100,zN/2]);

end

function testWIBL1OutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-wibl1.mat");
    hbar = data.y(1:end/3);
    F1bar = data.y(1+end/3:2/3*end);
    F2bar = data.y(1+2*end/3:end);


    alpha = 1;

    actual = linearisedWIBL1Fourier(domain, hbar, F1bar, F2bar, params, alpha);

    verifySize(testCase, actual, 3*[zN/2,zN/2]);

end

function testWIBL1DROutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-wibl1.mat");
    hbar = data.y(1:end/3);
    F1bar = data.y(1+end/3:2/3*end);
    F2bar = data.y(1+2*end/3:end);

    alpha = linspace(1e-3,1)';
    modes = zN/2;

    actual = dispersionRelationWIBL1Fourier(domain, hbar, F1bar, F2bar, params, alpha, modes)

    verifySize(testCase, actual, [100,zN/2]);

end

function testWIBL2STFOutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-wibl2-stf.mat");
    hbar = data.y(1:end/3);
    F1bar = data.y(1+end/3:2/3*end);
    F2bar = data.y(1+2*end/3:end);


    alpha = 1;

    actual = linearisedWIBL2STFFourier(domain, hbar, F1bar, F2bar, params, alpha);

    verifySize(testCase, actual, 3*[zN/2,zN/2]);

end

function testWIBL2STFDROutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-wibl2-stf.mat");
    hbar = data.y(1:end/3);
    F1bar = data.y(1+end/3:2/3*end);
    F2bar = data.y(1+2*end/3:end);

    alpha = linspace(1e-3,1)';
    modes = zN/2;

    actual = dispersionRelationWIBL2STFFourier(domain, hbar, F1bar, F2bar, params, alpha, modes)

    verifySize(testCase, actual, [100,zN/2]);

end

function testWIBL2OutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-wibl2.mat");
    hbar = data.y(1:end/3);
    F1bar = data.y(1+end/3:2/3*end);
    F2bar = data.y(1+2*end/3:end);


    alpha = 1;

    actual = linearisedWIBL2Fourier(domain, hbar, F1bar, F2bar, params, alpha);

    verifySize(testCase, actual, 3*[zN/2,zN/2]);

end

function testWIBL2DROutputSize(testCase)
    addpath("discretisationMethods");

    zN = 64;
    zL = 32;
    z = {linspace(zL/zN,zL,zN)'};

    domain = PSDomain(z, true, false);

    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);

    data = load("data/ic-rivulet-wibl2.mat");
    hbar = data.y(1:end/3);
    F1bar = data.y(1+end/3:2/3*end);
    F2bar = data.y(1+2*end/3:end);

    alpha = linspace(1e-3,1)';
    modes = zN/2;

    actual = dispersionRelationWIBL2Fourier(domain, hbar, F1bar, F2bar, params, alpha, modes)

    verifySize(testCase, actual, [100,zN/2]);

end
