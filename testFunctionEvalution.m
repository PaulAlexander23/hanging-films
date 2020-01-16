function tests = testFunctionEvalution()
    tests = functiontests(localfunctions);
end

%% Burgers function evaluation
function testFiniteDifferenceFburgersSize(testCase)
    addpath discretisationMethods
    N = 2^6;
    x = {linspace(2*pi/N, 2*pi, N)'};
    domain = FDDomain(x, [1, 2], 4);
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
        error(n) = max(abs(actual-expected), [], [1, 2]);
    end

    %     figure; plot(4:9, log10(error));

    verifyEqual(testCase, error(6), zeros(1, 1), 'AbsTol', 1e-4);

    function f = eval(N)
        x = {linspace(2*pi/N, 2*pi, N)'};
        domain = FDDomain(x, [1, 2], 4);
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

        f = domain.zeropad(f, 2^10/N);
    end
end

%% Benney function evaluation
function testFiniteDifferenceFbenney1dSize(testCase)
    addpath discretisationMethods
    N = 2^6;
    x = {linspace(2*pi/N, 2*pi, N)'};
    domain = FDDomain(x, [1, 2], 4);
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
        error(n) = max(abs(actual-expected), [], [1, 2]);
    end

    %     figure; plot(4:9, log10(error));

    verifyEqual(testCase, error(6), zeros(1, 1), 'AbsTol', 1e-3);

    function f = eval(N)
        x = {linspace(2*pi/N, 2*pi, N)'};
        domain = FDDomain(x, [1, 2], 4);
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
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2}');
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    actual = fbenney2d(domain, y, params);

    expectedSize = [2^8, 2^8];

    verifySize(testCase, actual, expectedSize)
end

function testFiniteDifferenceFbenney2dResolution(testCase)
    addpath discretisationMethods

    expected = eval(2^10);
    error = zeros(5, 1);
    for n = 1:6
        actual = eval(2^(n + 3));
        error(n) = max(abs(actual-expected), [], [1, 2]);
    end

    %     figure; plot(4:9, log10(error));
    %         figure; plot(4:9, log10(error/max(abs(actual),[],[1,2])));

    verifyEqual(testCase, error(6), zeros(1, 1), 'AbsTol', 1e1);

    function f = eval(N)
        diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        domain = FDDomain(setupX(1, 1, N, N), diffDegrees, 4);
        y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2}');
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

        f = fbenney2d(domain, y, params);

        f = interp2(domain.x{1}, domain.x{2}, ...
            f, ...
            linspace(2*pi/2^10, 2*pi, 2^10)', linspace(2*pi/2^10, 2*pi, 2^10)');
    end
end

function testPseudoSpectralFbenney2dSize(testCase)
    addpath discretisationMethods
    domain = PSDomain(setupX(1, 1, 2^8, 2^8));
    y = icos(domain.x);
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    actual = fbenney2d(domain, y, params);

    expectedSize = [2^8, 2^8];

    verifySize(testCase, actual, expectedSize)
end

function testPseudoSpectralFbenney2dResolutionDefault(testCase)
    addpath discretisationMethods /

    expected = eval(2^10);
    error = zeros(5, 1);
    for n = 1:6
        actual = eval(2^(n + 3));
        error(n) = max(abs(actual-expected), [], [1, 2]);
    end

    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error./max(abs(actual),[],[1,2])));

    verifyEqual(testCase, error(6), zeros(1, 1), 'AbsTol', 1e3);
    function f = eval(N)
        domain = PSDomain(setupX(1, 1, N, N));
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.fft(y), params);
        f = domain.zeropad(f, 2^10/N);
    end
end

function testPseudoSpectralFbenney2dResolutionAntiAliasing(testCase)
    addpath discretisationMethods /

    N = 2.^(5:8);
    expected = eval(N(end), N(end));
    abserror = zeros(length(N), 1);
    relerror = zeros(length(N), 1);
    for n = 1:length(N)
        actual = eval(N(n), N(end));
        abserror(n) = max(abs(actual-expected), [], [1, 2]);
        relerror(n) = max(abs(actual-expected)./abs(expected), [], [1, 2]);
    end

    %     figure; plot(log2(N), log10(relerror));
    %     figure; plot(log2(N), log10(error./max(abs(actual),[],[1,2])));

    verifyEqual(testCase, abserror(end), zeros(1, 1), 'AbsTol', 1e3);
    function f = eval(N, Nmax)
        domain = PSDomain(setupX(1, 1, N, N), true);
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.fft(y), params);
        f = domain.zeropad(f, Nmax/N);
    end
end

function testPseudoSpectralFbenney2dResolutionAntiAliasingReal(testCase)
    addpath discretisationMethods /

    N = 2.^(5:8);
    expected = eval(N(end), N(end));
    error = zeros(length(N), 1);
    for n = 1:length(N)
        actual = eval(N(n), N(end));
        error(n) = max(abs(actual-expected), [], [1, 2]);
    end

    %     figure; plot(log2(N), log2(error));
    %     figure; plot(4:9, log10(error./max(abs(actual),[],[1,2])));

    verifyEqual(testCase, error(end), zeros(1, 1), 'AbsTol', 1e3);
    function f = eval(N, Nmax)
        domain = PSDomain(setupX(1, 1, N, N), true, false);
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.fft(y), params);
        f = domain.zeropad(f, Nmax/N);
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

    actual = fbenney2d(domain, y, params);
    expected = domainPS.ifft( ...
        fbenney2d(domainPS, domainPS.fft(y), params));

    %     figure; surf(log10(abs(fft2(actualFD - actualPS))))
    %     figure; surf(log10(abs(fft2(actualFD - actualPS))/max(abs(actualFD),[],[1,2])))

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-1, 'AbsTol', 1e-2)
end

function testFiniteDifferenceEqualsPseudoSpectralFbenney2dDiagonal(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    domain = FDDomain(setupX(1, 1, 2^8, 2^8), diffDegrees, 4);
    domainPS = PSDomain(setupX(1, 1, 2^8, 2^8));
    y = 1 + 0.5 * cos(2*pi*domain.x{1}+2*pi*domain.x{2}');

    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    actual = fbenney2d(domain, y, params);
    expected = domainPS.ifft( ...
        fbenney2d(domainPS, domainPS.fft(y), params));

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-2, 'AbsTol', 2e-1)
end

%% WIBL1
function testFiniteDifferenceWIBL1Size(testCase)
    addpath discretisationMethods

    domain = FDDomain(setupX(1, 1, 2^6, 2^6), [1, 0; 0, 1; 2, 0; 0, 2]', 4);
    y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
    f = 2 / 3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
    Y = [y; f];
    params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);

    actual = fwibl1(domain, Y, params);
    expectedSize = [2^7, 2^6];

    verifySize(testCase, actual, expectedSize)
end

function testFiniteDifferenceWIBL1Resolution(testCase)
    addpath discretisationMethods

    expected = eval(2^10);
    error = zeros(5, 1);
    for n = 1:6
        actual = eval(2^(n + 3));
        error(n) = max(abs(actual-expected), [], [1, 2]);
    end

    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error/max(abs(actual),[],[1,2])));

    verifyEqual(testCase, error(6), zeros(1, 1), 'AbsTol', 1e2);

    function f = eval(N)
        diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        domain = FDDomain(setupX(1, 1, N, N), diffDegrees, 4);
        y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
        F1 = 2 / 3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
        Y = [y; F1];
        params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);

        f = fwibl1(domain, Y, params);

        newy = myinterp(domain, f(1:end/2, :));
        newF1 = myinterp(domain, f(end/2+1:end, :));

        f = [newy; newF1];
    end

    function f = myinterp(domain, f)
        f = interp2(domain.x{1}, domain.x{2}, ...
            f, ...
            linspace(2*pi/2^10, 2*pi, 2^10)', linspace(2*pi/2^10, 2*pi, 2^10)');
    end
end

function testPseudoSpectralWIBL1Size(testCase)
    addpath discretisationMethods
    domain = PSDomain(setupX(1, 1, 2^6, 2^6));
    y = domain.fft(1+0.1*cos(4*pi*domain.x{1})+0.1*cos(4*pi*domain.x{2}'));
    f = domain.fft(2/3+0.1*cos(4*pi*domain.x{1})+0.1*cos(4*pi*domain.x{2}'));
    Y = [y; f];
    params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);

    actual = fwibl1(domain, Y, params);
    expectedSize = [2^7, 2^6];

    verifySize(testCase, actual, expectedSize)
end

function testPseudoSpectralWIBL1Resolution(testCase)
    addpath discretisationMethods


    expected = eval(2^8);
    error = zeros(3, 1);
    for n = 1:4
        actual = eval(2^(n + 3));
        error(n) = max(abs(actual-expected), [], [1, 2]);
    end

    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error/max(abs(actual),[],[1,2])));

    verifyEqual(testCase, error(end), zeros(1, 1), 'AbsTol', 5e3);

    function f = eval(N)
        domain = PSDomain(setupX(1, 1, N, N), true, false);
        y = domain.fft(1+0.1*cos(4*pi*domain.x{1})+0.1*cos(4*pi*domain.x{2}'));
        F1 = domain.fft(2/3+0.1*cos(4*pi*domain.x{1})+0.1*cos(4*pi*domain.x{2}'));
        Y = [y; F1];
        params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);

        f = fwibl1(domain, Y, params);

        newy = domain.zeropad(f(1:end/2, :), 2^10/N);
        newF1 = domain.zeropad(f(end/2+1:end, :), 2^10/N);

        f = [newy; newF1];
    end
end

function testFiniteDifferenceEqualsPseudoSpectralFwibl12d(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    domain = FDDomain(setupX(1, 1, 2^8, 2^8), diffDegrees, 4);
    domainPS = PSDomain(setupX(1, 1, 2^8, 2^8));
    y = 1 + 0.1 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}'));
    F = 2 / 3 + 0.1 * 2 / 3 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}'));
    Y = [y; F];
    fY = [domainPS.fft(y); domainPS.fft(F)];
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    actualFD = fwibl1(domain, Y, params);
    Z = fwibl1(domainPS, fY, params);
    actualPS = [domainPS.ifft(Z(1:end/2, :)); ...
        domainPS.ifft(Z(1+end/2:end, :))];

    verifyEqual(testCase, actualFD, actualPS, 'RelTol', 6e-4, 'AbsTol', 1e-3)
end