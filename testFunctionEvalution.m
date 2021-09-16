function tests = testFunctionEvalution()
    tests = functiontests(localfunctions);
end

%% Burgers function evaluation
function testFiniteDifferenceFburgersSize(testCase)
    addpath discretisationMethods
    N = 4;
    domain = tFDDomain(N);
    y = cos(domain.x{1});
    params = struct('nu', 0.1);

    actual = fburgers(domain, y, params);
    expectedSize = [N, 1];

    verifySize(testCase, actual, expectedSize)
end

function testFiniteDifferenceFburgers1dResolution(testCase)
    addpath discretisationMethods

    expected = eval(2^10);
    error = zeros(5, 1);
    for n = 1:6
        actual = eval(2^(n + 3));
        error(n) = max(max(abs(actual-expected)));
    end

    %     figure; plot(4:9, log10(error));

    verifyEqual(testCase, error(6), zeros(1, 1), 'AbsTol', 1e-4);

    function f = eval(N)
        domain = tFDDomain(N);
        y = cos(domain.x{1});
        params = struct('nu', 0.1);

        f = fburgers(domain, y, params);

        f = interp1(domain.x{1}, f, linspace(2*pi/2^10, 2*pi, 2^10)');
    end
end

function testPseudoSpectralFburgersSize(testCase)
    addpath discretisationMethods
    N = 2^6;
    x = {linspace(2*pi/N, 2*pi, N)'};
    domain = PSDomain(x);
    y = domain.fft(cos(domain.x{1}));
    params = struct('nu', 0.1);

    actual = fburgers(domain, y, params);
    expectedSize = [N, 1];

    verifySize(testCase, actual, expectedSize)
end

function testPseudoSpectralFburgers1dResolution(testCase)
    addpath discretisationMethods

    expected = eval(2^10);
    error = zeros(6, 1);
    for n = 1:6
        actual = eval(2^(n + 3));
        error(n) = max(abs(actual-expected));
    end

    %     figure; plot(4:9, log10(error));
    %         figure; plot(log10(abs(actual(1:end/2)))); hold on; plot(log10(abs(expected(1:end/2))));
    %         figure; plot(log10(abs(actual(1:end/2)-expected(1:end/2))));

    verifyEqual(testCase, error(6), zeros(1, 1), 'AbsTol', 1e-2);

    function f = eval(N)
        x = {linspace(2*pi/N, 2*pi, N)'};
        domain = PSDomain(x);
        y = domain.fft(cos(domain.x{1}));
        params = struct('nu', 0.1);

        f = fburgers(domain, y, params);

        f = domain.ifft(f);
        f = interp1(domain.x{1}, f, linspace(2*pi/2^10, 2*pi, 2^10)');
    end
end

%% Benney function evaluation
function testFiniteDifferenceFbenney1dSize(testCase)
    addpath discretisationMethods
    N = 4;
    domain = tFDDomain(4);
    y = 1 + 0.25 * cos(domain.x{1});
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    actual = fbenney1d(domain, y, params);
    expectedSize = [N, 1];

    verifySize(testCase, actual, expectedSize)
end

function testFiniteDifferenceFbenney1dResolution(testCase)
    addpath discretisationMethods

    expected = eval(2^10);
    error = zeros(5, 1);
    for n = 1:6
        actual = eval(2^(n + 3));
        error(n) = max(max(abs(actual-expected)));
    end

    %     figure; plot(4:9, log10(error));

    verifyEqual(testCase, error(6), zeros(1, 1), 'AbsTol', 1e-3);

    function f = eval(N)
        domain = tFDDomain(N);
        y = 1 + 0.25 * cos(domain.x{1});
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

        f = fbenney1d(domain, y, params);

        f = interp1(domain.x{1}, f, linspace(2*pi/2^10, 2*pi, 2^10)');
    end
end

function testFiniteDifferenceFbenney2dSize(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    domain = FDDomain(setupX(1, 1, 2^8, 2^8), diffDegrees, 4);
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2});
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    y = domain.reshapeToVector(y);
    actual = fbenney2d(domain, y, params);

    expectedSize = [2^8 * 2^8, 1];

    verifySize(testCase, actual, expectedSize)
end

function testFiniteDifferenceFbenney2dConvergenceBetweenResolutions(testCase)
    addpath discretisationMethods

    for expectedOrder = [2, 4]
        resolutions = round(logspace(log10(50),log10(250),6));

        N = length(resolutions);
        errNorm = ones(N-1, 1);

        expected = myEval(resolutions(1), expectedOrder, resolutions(1), resolutions(end));
        for n = 2:N
            actual = myEval(resolutions(n), expectedOrder, resolutions(1),resolutions(end));

            errNorm(n-1) = max(max(abs(actual - expected)));

            expected = actual;
        end

        % hold on
        % plot(resolutions(2:end), log10(errNorm), 'o');
        % set(gca, 'Xscale', 'log')

        actualOrder = -(diff(log10(errNorm.')) ./ diff(log10(resolutions(1:N-1))));

        %  mean(actualOrder)
        %  errNorm(end)
        %  max(abs(actual),[],[1,2])
        verifyTrue(testCase, mean(actualOrder) > expectedOrder - 3e-1)
    end

    function f = myEval(N, order, minN, maxN)
        diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        domain = FDDomain(setupX(32, 32, N, N), diffDegrees, order);
        fineX = setupX(32, 32, maxN, maxN);
        coarseX = setupX(32, 32, minN, minN);
        dt = 1e-3;

        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

        % y = icos(fineX);
        % y = irand(fineX, 3e-1, 4, 1);

        %figure; surf(y)
        %y = interp2(fineX{1}, fineX{2}, ...
        %    y, ...
        %    domain.x{1}, domain.x{2});
        %figure; surf(y)

        % y = icos(domain.x);
        y = irandLin(domain.x, 3e-1, 20, 1);
        % y = iloadInterp(domain.x, "ics/ic-tw-benney-2d.mat");

        y = domain.reshapeToVector(y);
        f = dt * fbenney2d(domain, y, params);
        f = domain.reshapeToDomain(f);

        f = periodicInterp2(domain.x{1}, domain.x{2}, ...
            f, ...
            coarseX{1}, coarseX{2}, 'spline'); 
    end
end

function testPseudoSpectralFbenney2dSize(testCase)
    addpath discretisationMethods
    domain = PSDomain(setupX(1, 1, 2^8, 2^8));
    y = icos(domain.x);
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    y = domain.reshapeToVector(y);
    actual = fbenney2d(domain, y, params);

    expectedSize = [2^8 * 2^8, 1];

    verifySize(testCase, actual, expectedSize)
end

function testPseudoSpectralFbenney2dResolutionDefault(testCase)
    addpath discretisationMethods /

    p = 7;
    expected = eval(2^p);
    error = zeros(p-3, 1);
    for n = 1:p-3
        actual = eval(2^(n + 3));
        error(n) = max(max(abs(actual-expected)));
    end

    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error./max(abs(actual),[],[1,2])));

    verifyEqual(testCase, error(end), zeros(1, 1), 'AbsTol', 1e3);

    function f = eval(N)
        domain = PSDomain(setupX(1, 1, N, N));
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.reshapeToVector(domain.fft(y)), params);
        f = domain.reshapeToDomain(f);
        f = domain.ifft(f);
        f = interp2(domain.x{1}, domain.x{2}, ...
            f, ...
            linspace(2*pi/2^10, 2*pi, 2^10)', linspace(2*pi/2^10, 2*pi, 2^10));
    end
end

function testPseudoSpectralFbenney2dResolutionAntiAliasing(testCase)
    addpath discretisationMethods /

    N = 2.^(5:7);
    expected = eval(N(end), N(end));
    abserror = zeros(length(N), 1);
    relerror = zeros(length(N), 1);
    for n = 1:length(N)
        actual = eval(N(n), N(end));
        abserror(n) = max(max(abs(actual-expected)));
        relerror(n) = max(max(abs(actual-expected)./abs(expected)));
    end

    %     figure; plot(log2(N), log10(relerror));
    %     figure; plot(log2(N), log10(error./max(abs(actual),[],[1,2])));

    verifyEqual(testCase, abserror(end), zeros(1, 1), 'AbsTol', 1e3);
    function f = eval(N, Nmax)
        domain = PSDomain(setupX(1, 1, N, N), true);
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.reshapeToVector(domain.fft(y)), params);
        f = domain.reshapeToDomain(f);
        f = domain.ifft(f);
        f = interp2(domain.x{1}, domain.x{2}, ...
            f, ...
            linspace(2*pi/Nmax, 2*pi, Nmax)', linspace(2*pi/Nmax, 2*pi, Nmax));
    end
end

function testPseudoSpectralFbenney2dResolutionAntiAliasingReal(testCase)
    addpath discretisationMethods /

    N = 2.^(5:7);
    expected = eval(N(end), N(end));
    error = zeros(length(N), 1);
    for n = 1:length(N)
        actual = eval(N(n), N(end));
        error(n) = max(max(abs(actual-expected)));
    end

    %     figure; plot(log2(N), log2(error));
    %     figure; plot(4:9, log10(error./max(abs(actual),[],[1,2]))); 
    verifyEqual(testCase, error(end), zeros(1, 1), 'AbsTol', 1e3);
    function f = eval(N, Nmax)
        domain = PSDomain(setupX(1, 1, N, N), true, false);
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.reshapeToVector(domain.fft(y)), params);
        f = domain.reshapeToDomain(f);
        f = domain.ifft(f);
        f = interp2(domain.x{1}, domain.x{2}, ...
            f, ...
            linspace(2*pi/Nmax, 2*pi, Nmax)', linspace(2*pi/Nmax, 2*pi, Nmax));
    end
end

function testFiniteDifferenceEqualsPseudoSpectralFbenney2d(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    M = 2^7;
    N = 2^7;
    domain = FDDomain(setupX(1, 1, M, N), diffDegrees, 4);
    domainPS = PSDomain(setupX(1, 1, M, N));
    y = icos(domain.x);

    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    actual = fbenney2d(domain, domain.reshapeToVector(y), params);
    actual = domain.reshapeToDomain(actual);
    expected = domainPS.ifft( ...
        domain.reshapeToDomain(fbenney2d(domainPS, domain.reshapeToVector(domainPS.fft(y)), params)));

    %     figure; surf(log10(abs(fft2(actualFD - actualPS))))
    %     figure; surf(log10(abs(fft2(actualFD - actualPS))/max(abs(actualFD),[],[1,2])))

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-1, 'AbsTol', 1e-2)
end

function testFiniteDifferenceEqualsPseudoSpectralFbenney2dDiagonal(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    domain = FDDomain(setupX(1, 1, 2^8, 2^8), diffDegrees, 4);
    domainPS = PSDomain(setupX(1, 1, 2^8, 2^8));
    y = 1 + 0.5 * cos(2*pi*domain.x{1}+2*pi*domain.x{2});

    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    actual = fbenney2d(domain, domain.reshapeToVector(y), params);
    actual = domain.reshapeToDomain(actual);
    expected = domainPS.ifft( ...
        domain.reshapeToDomain(fbenney2d(domainPS, domain.reshapeToVector(domainPS.fft(y)), params)));

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-2, 'AbsTol', 2e-1)
end

%% WIBL1STF
function testFiniteDifferenceWIBL1Size(testCase)
    addpath discretisationMethods

    domain = FDDomain(setupX(1, 1, 2^6, 2^6), [1, 0; 0, 1; 2, 0; 0, 2; 3, 0; 1, 2]', 4);
    y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2});
    f = 2 / 3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2});
    Y = [y; f];
    params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);

    Y = domain.reshapeToVector(Y);
    actual = fwibl1STF(domain, Y, params);
    expectedSize = [2^7 * 2^6, 1];

    verifySize(testCase, actual, expectedSize)
end

function testFiniteDifferenceWIBL1Resolution(testCase)
    addpath discretisationMethods

    p = 7;
    expected = eval(2^p);
    error = zeros(p-3, 1);
    for n = 1:p-3
        actual = eval(2^(n + 3));
        error(n) = max(max(abs(actual-expected)));
    end

    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error/max(abs(actual),[],[1,2])));

    verifyEqual(testCase, error(end), zeros(1, 1), 'AbsTol', 1e2);

    function f = eval(N)
        diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2; 3, 0; 1, 2]';
        domain = FDDomain(setupX(1, 1, N, N), diffDegrees, 4);
        y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2});
        F1 = 2 / 3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2});
        Y = [y; F1];
        params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);

        Y = domain.reshapeToVector(Y);
        f = fwibl1STF(domain, Y, params);
        f = domain.reshapeToDomain(f);

        newy = myinterp(domain, f(1:end/2, :));
        newF1 = myinterp(domain, f(end/2+1:end, :));

        f = [newy; newF1];
    end

    function f = myinterp(domain, f)
        f = interp2(domain.x{1}, domain.x{2}, ...
            f, ...
            linspace(2*pi/2^10, 2*pi, 2^10)', linspace(2*pi/2^10, 2*pi, 2^10));
    end
end

function testFiniteDifferenceWIBL1ConvergenceBetweenResolutions(testCase)
    addpath discretisationMethods

    for expectedOrder = [2, 4]
        resolutions = 2.^(5:8);

        N = length(resolutions);
        errNorm = ones(N-1, 1);

        expected = myEval(resolutions(1), expectedOrder, resolutions(1), resolutions(end));
        for n = 2:N
            actual = myEval(resolutions(n), expectedOrder, resolutions(1),resolutions(end));

            errNorm(n-1) = max(max(abs(actual - expected)));

            expected = actual;
        end

        actualOrder = -(gradient(log10(errNorm), log10(resolutions(1:N-1))));

        verifyTrue(testCase,all(actualOrder > expectedOrder - 9e-1))
    end

    function f = myEval(N, order, minN, maxN)
        diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2; 3, 0; 1, 2]';
        domain = FDDomain(setupX(1, 1, N, N), diffDegrees, order);
        fineX = setupX(1, 1, maxN, maxN);

        y = icosWIBL1(fineX);
        % y = irandWIBL1(fineX, 3e-1, 4, 1);

        y = [interp2(fineX{1}, fineX{2}, ...
            y(1:end/2,:), ...
            domain.x{1}, domain.x{2}); ...
            interp2(fineX{1}, fineX{2}, ...
            y(1+end/2:end,:), ...
            domain.x{1}, domain.x{2})]; 

        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

        y = domain.reshapeToVector(y);
        f = fwibl1STF(domain, y, params);
        f = domain.reshapeToDomain(f);

        coarseX = setupX(1, 1, minN, minN);

        f = [interp2(domain.x{1}, domain.x{2}, ...
            f(1:end/2,:), ...
            coarseX{1}, coarseX{2}); ...
            interp2(domain.x{1}, domain.x{2}, ...
            f(1+end/2:end,:), ...
            coarseX{1}, coarseX{2})]; 
    end
end

function testPseudoSpectralWIBL1Size(testCase)
    addpath discretisationMethods
    domain = PSDomain(setupX(1, 1, 2^6, 2^6));
    y = domain.fft(1+0.1*cos(4*pi*domain.x{1})+0.1*cos(4*pi*domain.x{2}));
    f = domain.fft(2/3+0.1*cos(4*pi*domain.x{1})+0.1*cos(4*pi*domain.x{2}));
    Y = [y; f];
    params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);

    Y = domain.reshapeToVector(Y);
    actual = fwibl1STF(domain, Y, params);
    expectedSize = [2^7 * 2^6, 1];

    verifySize(testCase, actual, expectedSize)
end

function testPseudoSpectralWIBL1Resolution(testCase)
    addpath discretisationMethods


    expected = eval(2^8);
    error = zeros(3, 1);
    for n = 1:4
        actual = eval(2^(n + 3));
        error(n) = max(max(abs(actual-expected)));
    end

    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error/max(abs(actual),[],[1,2])));

    verifyEqual(testCase, error(end), zeros(1, 1), 'AbsTol', 5e3);

    function f = eval(N)
        domain = PSDomain(setupX(1, 1, N, N), true, false);
        y = domain.fft(1+0.1*cos(4*pi*domain.x{1})+0.1*cos(4*pi*domain.x{2}));
        F1 = domain.fft(2/3+0.1*cos(4*pi*domain.x{1})+0.1*cos(4*pi*domain.x{2}));
        Y = [y; F1];
        params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);

        Y = domain.reshapeToVector(Y);
        f = fwibl1STF(domain, Y, params);
        f = domain.reshapeToDomain(f);

        f = domain.ifft(f);
        newy  = interp2(domain.x{1}, domain.x{2}, ...
            f(1:end/2, :), ...
            linspace(2*pi/2^8, 2*pi, 2^8)', linspace(2*pi/2^8, 2*pi, 2^8));
        newF1= interp2(domain.x{1}, domain.x{2}, ...
            f(end/2+1:end, :), ...
            linspace(2*pi/2^8, 2*pi, 2^8)', linspace(2*pi/2^8, 2*pi, 2^8));

        f = [newy; newF1];
    end
end

function testFiniteDifferenceEqualsPseudoSpectralFwibl12d(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2; 3, 0; 1, 2]';
    domain = FDDomain(setupX(1, 1, 2^8, 2^8), diffDegrees, 4);
    domainPS = PSDomain(setupX(1, 1, 2^8, 2^8));
    y = 1 + 0.1 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}));
    F = 2 / 3 + 0.1 * 2 / 3 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}));
    Y = [y; F];
    fY = [domainPS.fft(y); domainPS.fft(F)];
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    Y = domain.reshapeToVector(Y);
    actualFD = fwibl1STF(domain, Y, params);
    actualFD = domain.reshapeToDomain(actualFD);
    fY = domain.reshapeToVector(fY);
    Z = fwibl1STF(domainPS, fY, params);
    Z = domain.reshapeToDomain(Z);
    actualPS = [domainPS.ifft(Z(1:end/2, :)); ...
        domainPS.ifft(Z(1+end/2:end, :))];

    verifyEqual(testCase, actualFD, actualPS, 'RelTol', 6e-4, 'AbsTol', 1e-3)
end

%% Hybrid

function testHybridEqualsBenney(testCase)
    addpath discretisationMethods

    domain = FDDomain(setupX(1, 1, 2^6, 2^6), [1, 0; 0, 1; 2, 0; 0, 2]', 4);
    y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2});

    params = struct('theta', pi/4, 'Re', 1, 'C', 0.01, ...
        'epsilon', 1, 'delta', 0);
    P = 2 * cot(params.theta) * y - ...
        (domain.diff(y, [2; 0]) + domain.diff(y, [0; 2])) / params.C;
    f = 2 / 3 * domain.multiply(y, y, [2, 1]) - ...
        1 / 3 * domain.multiply(y, domain.diff(P, [1; 0]), [3, 1]) + ...
        8 * params.Re / 15 * domain.multiply(y, domain.diff(y, [1; 0]), [6, 1]);

    Y = [y; f];

    Y = domain.reshapeToVector(Y);
    expected = fbenney2d(domain, domain.reshapeToVector(y), params);
    actual = fhybrid(domain, Y, params);

    verifyEqual(testCase, actual(1:end/2,:), expected, 'RelTol', 1e-11)
end

function testHybridEqualsWIBL1(testCase)
    addpath discretisationMethods

    domain = FDDomain(setupX(1, 1, 2^6, 2^6), [1, 0; 0, 1; 2, 0; 0, 2; 3, 0; 1, 2]', 4);
    y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2});
    f = 2 / 3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2});

    params = struct('theta', pi/4, 'Re', 1, 'C', 0.01, ...
        'epsilon', 0, 'delta', 1);

    Y = [y; f];

    Y = domain.reshapeToVector(Y);
    expected = fwibl1STF(domain, Y, params);
    actual = fhybrid(domain, Y, params);

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-5, 'RelTol', 1e-4)
end

function domain = tFDDomain(N)
    if nargin < 1, N = 2^6; end
    x = {linspace(2*pi/N, 2*pi, N)'};
    domain = FDDomain(x, [1, 2], 4);
end
