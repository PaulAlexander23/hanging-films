function tests = testODESolver()
    tests = functiontests(localfunctions);
end

%% odeMatrixSolver
function testPDESolver1DExponentialDecay(testCase)
    addpath discretisationMethods
    odefun = @(t, y) - y;
    t = linspace(0, 1, 11);
    x = linspace(2*pi/64, 2*pi, 64)';
    domain = Domain({x});
    y0 = cos(domain.x{1});
    timestepper = @ode45;
    actual = odeMatrixSolver(odefun, t, y0, timestepper);
    expected = y0 * exp(-t);
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolver1DTravellingWaveEquation(testCase)
    addpath discretisationMethods
    t = linspace(0, 1, 11);
    domain = FDDomain({linspace(2*pi/64, 2*pi, 64)'}, 1, 4);
    odefun = @(t, y) - domain.diff(y, 1);
    y0 = cos(domain.x{1});
    timestepper = @ode45;

    actual = odeMatrixSolver(odefun, t, y0, timestepper);
    expected = cos(domain.x{1}-t);

    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolver1DHeatEquation(testCase)
    addpath discretisationMethods
    t = linspace(0, 1, 11);
    x = {linspace(2*pi/64, 2*pi, 64)'};
    domain = FDDomain(x, 2, 4);
    odefun = @(t, y) domain.diff(y, 2);
    y0 = cos(domain.x{1});
    timestepper = @ode45;
    actual = odeMatrixSolver(odefun, t, y0, timestepper);
    expected = y0 * exp(-t);
    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolver1DHeatEquationPseudoSpectral(testCase)
    addpath discretisationMethods
    t = linspace(0, 1, 11);
    x = {linspace(2*pi/64, 2*pi, 64)'};
    domain = PSDomain(x);
    odefun = @(t, y) domain.diff(y, 2);
    y0 = cos(domain.x{1});
    timestepper = @ode45;
    actual = real(domain.ifft(odeMatrixSolver(odefun, t, domain.fft(y0), timestepper)));
    expected = y0 * exp(-t);
    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolver1DBenneyEquation(testCase)
    addpath discretisationMethods
    t = linspace(0, 1, 11);
    x = {linspace(2*pi/64, 2*pi, 64)'};
    domain = FDDomain(x, [1, 2], 4);
    y0 = 1 + 0.5 * cos(domain.x{1});
    params = struct('theta', 1, 'Re', 1, 'C', 1);
    odefun = @(t, y) fbenney1d(domain, y, params);
    timestepper = @ode15s;
    actual = odeMatrixSolver(odefun, t, y0, timestepper);
    load('data/testPDESolver1DBenneyEquationExpected', 'expected')
    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolverBurgers(testCase)
    addpath discretisationMethods /
    addpath ../time-stepping-methods /
    N = 2^4;
    x = {linspace(2*pi/N, 2*pi, N)'};
    domain = PSDomain(x, true, false);
    y0 = domain.fft(cos(domain.x{1}));
    params = struct('nu', 10);
    tFinal = 1e-4;
    odefun = @(t, y) fburgers(domain, y, params);

    y = ab1(odefun, 0:1e-5:tFinal, y0);
    t = 0:1e-5:tFinal;
    %     plot(t)
    %     figure; plot(log10(abs(y(1:end/2,1:10:end))));
    %     figure; plot(log10(abs(y(1:end/2,end))));
    %     figure; plot(domain.ifft(y(:,end)));

    verifyEqual(testCase, t(end), tFinal)
end

function testPDESolverKDV(testCase)
    addpath discretisationMethods /
    addpath ../time-stepping-methods/
    N = 2^8;
    x = {linspace(2*pi/N, 2*pi, N)'};
    domain = PSDomain(x, true, false);
    y0 = domain.fft(cos(domain.x{1}));
    params = struct('nu', 1);

    odefun = @(t, y) fkdv(domain, y, params);
    tFinal = 70e-6;
    y = ab1(odefun, 0:1e-6:tFinal, y0);
    t = 0:1e-5:tFinal;
    %     plot(t)
    %     figure; plot(log10(abs(y(1:end/2,1:10:end))));
    %     figure; plot(log10(abs(y(1:end/2,end))));
    %     figure; plot(domain.ifft(y(:,end)));
    verifyEqual(testCase, t(end), tFinal)
end
