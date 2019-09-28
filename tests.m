function tests = tests()
    tests = functiontests(localfunctions);
end

%% Independent variable setup
function testSetupX(testCase)
    actual = setupX(10,pi,12,11);
    verifyEqual(testCase, actual{1}, linspace(10/12,10,12)')
    verifyEqual(testCase, actual{2}, linspace(pi/11,pi,11)')
end

function testSetupT(testCase)
    actual = setupT(10,0.2);
    verifyEqual(testCase, actual(end), 10);
    verifyTrue(testCase, all(abs(diff(actual(1:end-1)) - 0.2)<1e-14));
end

%% Interface evalution

function testInterfaceCos(testCase)
    x = setupX(1,1,2^8,2^8);
    
    actual = icos(x);
    expected = 1 - 0.25 * cos(2*pi/x{1}(end) * x{1}) - ...
        0.25 * cos(2*pi/x{2}(end) * x{2}');
    verifyEqual(testCase, actual, expected);
    
    actual = icos(x, 0.5);
    expected = 1 - 0.5 * cos(2*pi/x{1}(end) * x{1}) - ...
        0.25 * cos(2*pi/x{2}(end) * x{2}');
    verifyEqual(testCase, actual, expected);
    
    actual = icos(x, 0.5, 0.5);
    expected = 1 - 0.5 * cos(2*pi/x{1}(end) * x{1}) - ...
        0.5 * cos(2*pi/x{2}(end) * x{2}');
    verifyEqual(testCase, actual, expected);
end

function testInterfaceSingle(testCase)
    x = setupX(1,1,2^8,2^8);
    A = 1e-2;
    
    actual = isingle(x, A);
    expected = 1 - A * (cos(2*pi/x{1}(end) * x{1}) + cos(2*pi/x{2}(end) * x{2}'));
    verifyEqual(testCase, actual, expected)
end

function testInterfaceRivulet(testCase)
    x = setupX(1,1,2^8,2^8);
    A = 2e-1;
    r = 0.05;
    
    actual = irivulet(x, A, r);
    expected = 1 + A * (-r*cos(2*pi/x{1}(end) * x{1}) - cos(2*pi/x{2}(end) * x{2}'));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceDoubleX(testCase)
    x = setupX(1,1,2^8,2^8);
    
    actual = idoublex(x);
    expected = 1 + 0.2 * (-0.05*cos(4*pi/x{1}(end) * x{1}) - cos(2*pi/x{2}(end) * x{2}'));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceDoubleY(testCase)
    x = setupX(1,1,2^8,2^8);
    
    actual = idoubley(x);
    expected = 1 + 0.2 * (-0.05*cos(2*pi/x{1}(end) * x{1}) - cos(4*pi/x{2}(end) * x{2}'));
    verifyEqual(testCase, actual, expected);
end

function testInterfacePert(testCase)
    x = setupX(1,1,2^8,2^8);
    
    actual = ipert(x);
    expected = 1 + 1e-2 * exp(-1e-1*((x{1}-x{1}(end)/2).^2 + (x{2}'-x{2}(end)/2).^2));
    verifyEqual(testCase, actual, expected);
end

function testInterfaceRand(testCase)
    x = setupX(1,1,2^8,2^8);
    A = 1e-4;
    numberOfModes = 5;
    
    actual = irand(x, A, numberOfModes);
    verifyTrue(testCase, max(max(abs(actual-1))) <= A);
    verifyEqual(testCase,sum(abs(fft2(actual-1))>1e-10,[1,2]), ...
        (numberOfModes * 2)^2);
end

function testInterfaceLoad(testCase)
    x = setupX(1,1,2^8,2^8);
    
    actual = iload(x, 'testInterfaceLoad');
    expected = icos(x);
    verifyEqual(testCase, actual, expected);
end

function testInterfaceLoadWithPert(testCase)
    x = setupX(1,1,2^8,2^8);
    alpha = 0.2;
    
    actual = iloadWithPert(x, 'testInterfaceLoad', alpha);
    expected = icos(x) + ipert(x, alpha) - 1;
    verifyEqual(testCase, actual, expected);
end

function testInterfaceLoadWithRand(testCase)
    x = setupX(1,1,2^8,2^8);
    alpha = 0.2;
    
    actual = iloadWithRand(x, 'testInterfaceLoad', alpha);
    disturbance = actual - icos(x);
    verifyTrue(testCase, max(max(abs(disturbance))) <= alpha + 1e-13);
    verifyEqual(testCase,sum(abs(fft2(disturbance))>1e-10,[1,2]), ...
        (5 * 2)^2);
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
    error = zeros(5,1);
    for n = 1:6
        actual = eval(2^(n+3));
        error(n) = max(abs(actual-expected), [], [1,2]);
    end
    
    figure; plot(4:9, log10(error));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 1e-4);
    
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
    error = zeros(6,1);
    for n = 1:6
        actual = eval(2^(n+3));
        size(actual)
        error(n) = max(abs(actual-expected));
    end
    
    figure; plot(4:9, log10(error));
    %         figure; plot(log10(abs(actual(1:end/2)))); hold on; plot(log10(abs(expected(1:end/2))));
    %         figure; plot(log10(abs(actual(1:end/2)-expected(1:end/2))));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 1e-2);
    
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
    error = zeros(5,1);
    for n = 1:6
        actual = eval(2^(n+3));
        error(n) = max(abs(actual-expected), [], [1,2]);
    end
    
    figure; plot(4:9, log10(error));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 1e-3);
    
    function f = eval(N)
        x = {linspace(2*pi/N, 2*pi, N)'};
        domain = FDDomain(x, [1, 2], 4);
        y = 1 + 0.25*cos(domain.x{1});
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        
        f = fbenney1d(domain, y, params);
        
        f = interp1(domain.x{1}, f, linspace(2*pi/2^10, 2*pi, 2^10)');
    end
end

function testFiniteDifferenceFbenney2dSize(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    domain = FDDomain(setupX(1,1,2^8,2^8), diffDegrees, 4);
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2}');
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
    
    actual = fbenney2d(domain, y, params);
    
    expectedSize = [2^8, 2^8];
    
    verifySize(testCase, actual, expectedSize)
end

function testFiniteDifferenceFbenney2dResolution(testCase)
    addpath discretisationMethods
    
    expected = eval(2^10);
    error = zeros(5,1);
    for n = 1:6
        actual = eval(2^(n+3));
        error(n) = max(abs(actual-expected), [], [1,2]);
    end
    
    figure; plot(4:9, log10(error));
    %         figure; plot(4:9, log10(error/max(abs(actual),[],[1,2])));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 1e1);
    
    function f = eval(N)
        diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        domain = FDDomain(setupX(1,1,N,N), diffDegrees, 4);
        y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2}');
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        
        f = fbenney2d(domain, y, params);
        
        f = interp2(domain.x{1}, domain.x{2}, ...
            f, ...
            linspace(2*pi/2^10, 2*pi, 2^10)',linspace(2*pi/2^10, 2*pi, 2^10)');
    end
end

function testPseudoSpectralFbenney2dSize(testCase)
    addpath discretisationMethods
    domain = PSDomain(setupX(1,1,2^8,2^8));
    y = icos(domain.x);
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
    
    actual = fbenney2d(domain, y, params);
    
    expectedSize = [2^8, 2^8];
    
    verifySize(testCase, actual, expectedSize)
end

function testPseudoSpectralFbenney2dResolutionDefault(testCase)
    addpath discretisationMethods/
    
    expected = eval(2^10);
    error = zeros(5,1);
    for n = 1:6
        actual = eval(2^(n+3));
        error(n) = max(abs(actual-expected), [], [1,2]);
    end
    
    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error./max(abs(actual),[],[1,2])));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 1e3);
    function f = eval(N)
        domain = PSDomain(setupX(1,1,N,N));
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.fft(y), params);
        f = domain.zeropad(f, 2^10/N);
    end
end

function testPseudoSpectralFbenney2dResolutionNoAntiAliasing(testCase)
    addpath discretisationMethods/
    
    N = 2.^(5:11);
    expected = eval(N(end), N(end));
    error = zeros(length(N),1);
    for n = 1:length(N)
        actual = eval(N(n), N(end));
        error(n) = max(abs(actual-expected), [], [1,2]);
    end
    
    figure; plot(log2(N), log2(error));
    %     figure; plot(4:9, log10(error./max(abs(actual),[],[1,2])));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 1e3);
    function f = eval(N, Nmax)
        domain = PSDomain(setupX(1,1,N,N), nan, false, 2);
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.fft(y), params);
        f = domain.zeropad(f, Nmax/N);
    end
end

function testPseudoSpectralFbenney2dResolutionNoShortwaveFilter(testCase)
    addpath discretisationMethods/
    
    expected = eval(2^10);
    error = zeros(5,1);
    for n = 1:6
        actual = eval(2^(n+3));
        error(n) = max(abs(actual-expected), [], [1,2]);
    end
    
    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error./max(abs(actual),[],[1,2])));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 1e3);
    function f = eval(N)
        domain = PSDomain(setupX(1,1,N,N), 1e-13, true, 1);
        params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
        y = icos(domain.x);
        f = fbenney2d(domain, domain.fft(y), params);
        f = domain.zeropad(f, 2^10/N);
    end
end

function testFiniteDifferenceEqualsPseudoSpectralFbenney2d(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    M = 2^7; N = 2^7;
    domain = FDDomain(setupX(1,1,M,N), diffDegrees, 4);
    domainPS = PSDomain(setupX(1,1,M,N));
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
    domain = FDDomain(setupX(1,1,2^8,2^8), diffDegrees, 4);
    domainPS = PSDomain(setupX(1,1,2^8,2^8));
    y = 1 + 0.5*cos(2*pi*domain.x{1} + 2*pi*domain.x{2}');
    
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
    
    actual = fbenney2d(domain, y, params);
    expected = domainPS.ifft( ...
        fbenney2d(domainPS, domainPS.fft(y), params));
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-2, 'AbsTol', 2e-1)
end

%% WIBL1 function evaluation
function testFiniteDifferenceWIBL1Size(testCase)
    addpath discretisationMethods
    
    domain = FDDomain(setupX(1,1,2^6,2^6),[1,0;0,1;2,0;0,2]',4);
    y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
    f = 2/3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
    Y = [y;f];
    params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);
    
    actual = fwibl1(domain, Y, params);
    expectedSize = [2^7, 2^6];
    
    verifySize(testCase, actual, expectedSize)
end

function testFiniteDifferenceWIBL1Resolution(testCase)
    addpath discretisationMethods
    
    expected = eval(2^10);
    error = zeros(5,1);
    for n = 1:6
        actual = eval(2^(n+3));
        error(n) = max(abs(actual-expected), [], [1,2]);
    end
    
    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error/max(abs(actual),[],[1,2])));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 1e2);
    
    function f = eval(N)
        diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
        domain = FDDomain(setupX(1,1,N,N), diffDegrees, 4);
        y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
        F1 = 2/3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
        Y = [y;F1];
        params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);
        
        f = fwibl1(domain, Y, params);
        
        newy = myinterp(domain, f(1:end/2,:));
        newF1 = myinterp(domain, f(end/2+1:end,:));
        
        f = [newy; newF1];
    end
    
    function f = myinterp(domain, f)
        f = interp2(domain.x{1}, domain.x{2}, ...
            f, ...
            linspace(2*pi/2^10, 2*pi, 2^10)',linspace(2*pi/2^10, 2*pi, 2^10)');
    end
end

function testPseudoSpectralWIBL1Size(testCase)
    addpath discretisationMethods
    domain = PSDomain(setupX(1,1,2^6,2^6));
    y = domain.fft(1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}'));
    f = domain.fft(2/3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}'));
    Y = [y;f];
    params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);
    
    actual = fwibl1(domain, Y, params);
    expectedSize = [2^7, 2^6];
    
    verifySize(testCase, actual, expectedSize)
end

function testPseudoSpectralWIBL1Resolution(testCase)
    addpath discretisationMethods
    
    
    expected = eval(2^10);
    error = zeros(5,1);
    for n = 1:6
        actual = eval(2^(n+3));
        error(n) = max(abs(actual-expected), [], [1,2]);
    end
    
    %     figure; plot(4:9, log10(error));
    %     figure; plot(4:9, log10(error/max(abs(actual),[],[1,2])));
    
    verifyEqual(testCase, error(6), zeros(1,1), 'AbsTol', 2e3);
    
    function f = eval(N)
        domain = PSDomain(setupX(1,1,N,N));
        y = domain.fft(1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}'));
        F1 = domain.fft(2/3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}'));
        Y = [y;F1];
        params = struct('theta', pi/4, 'Re', 1, 'C', 0.01);
        
        f = fwibl1(domain, Y, params);
        
        newy = domain.zeropad(f(1:end/2,:), 2^10/N);
        newF1 = domain.zeropad(f(end/2+1:end,:), 2^10/N);
        
        f = [newy; newF1];
    end
end

function testFiniteDifferenceEqualsPseudoSpectralFwibl12d(testCase)
    addpath discretisationMethods
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    domain = FDDomain(setupX(1,1,2^8,2^8), diffDegrees, 4);
    domainPS = PSDomain(setupX(1,1,2^8,2^8));
    y = 1 + 0.1 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}'));
    F = 2/3 + 0.1*2/3 * (cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}'));
    Y = [y;F];
    fY = [domainPS.fft(y);domainPS.fft(F)];
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
    
    actualFD = fwibl1(domain, Y, params);
    Z = fwibl1(domainPS, fY, params);
    actualPS = [domainPS.ifft(Z(1:end/2,:)); ...
        domainPS.ifft(Z(1+end/2:end,:))];
    
    verifyEqual(testCase, actualFD, actualPS, 'RelTol', 6e-4, 'AbsTol', 1e-3)
end

%% Data handling
function testMakeFilename(testCase)
    params = struct('theta', 2, 'Re', 3, 'C', 4);
    actual = makeFilename('-test',params,{0.1:0.1:3.3,1:1:10},pi,@interface,1e-6,"benney");
    expected = "data-test-theta-2-Re-3-C-4-xL-3_3-yL-10-T-3_14159-interface-interface-xN-33-yN-10-AbsTol-1e-06-model-benney";
    verifyEqual(testCase, actual, expected)
end

function testSaveAndLoadData(testCase)
    filename = 'data-test-theta-0_876255-Re-0_488629-C-0_407077-xL-3_3-yL-10-T-3_14159-interface-interface-xN-33-yN-10-AbsTol-1e-06-model-benney.mat';
    if isfile(filename)
        delete(filename)
    end
    
    expectedY = rand(100);
    expectedParams = struct('theta', 0.876255, 'Re', 0.488629, 'C', 0.407077);
    expectedT = linspace(0,0.126576313737181);
    expectedX = {0.1:0.1:3.3,1:1:10};
    expectedTimeTaken = 0.925425280986515;
    saveData(expectedY, expectedParams, expectedT, expectedX, ...
        expectedTimeTaken, pi, @interface, 1e-6, "benney", "-test");
    
    verifyTrue(testCase,isfile(filename))
    
    [actualY, actualParams, actualT, actualX, actualTimeTaken] = loadData(expectedParams,expectedX,pi,@interface,1e-6,"benney","-test");
    
    verifyEqual(testCase,actualY,expectedY)
    verifyEqual(testCase,actualParams,expectedParams)
    verifyEqual(testCase,actualT,expectedT)
    verifyEqual(testCase,actualX,expectedX)
    verifyEqual(testCase,actualTimeTaken,expectedTimeTaken)
    
    if isfile(filename)
        delete(filename)
    end
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
    odefun = @(t, y) -domain.diff(y, 1);
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
    y0 = 1+0.5*cos(domain.x{1});
    params = struct('theta', 1, 'Re', 1, 'C', 1);
    odefun = @(t, y) fbenney1d(domain, y, params);
    timestepper = @ode15s;
    actual = odeMatrixSolver(odefun, t, y0, timestepper);
    load('testPDESolver1DBenneyEquationExpected','expected')
    verifyEqual(testCase, actual(:, end), expected(:, end), ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testPDESolverBurgers(testCase)
    addpath discretisationMethods/
    addpath ../time-stepping-methods/
    N = 2^6;
    x = {linspace(2*pi/N, 2*pi, N)'};
    domain = PSDomain(x, nan, false, 1);
    y0 = domain.fft(cos(domain.x{1}));
    params = struct('nu', 0.1);
    
    odefun = @(t, y) fburgers(domain, y, params);
    
    [t, y] = ab1(odefun, 0:1e-3:20, y0);
    
    plot(t)
    figure; plot(log10(abs(y(1:end/2,1:10:end))));
    figure; plot(log10(abs(y(1:end/2,end))));
    figure; plot(domain.ifft(y(:,end)));
    verifyEqual(testCase, t(end), 1)
end

function testPDESolverKDV(testCase)
    addpath discretisationMethods/
    addpath ../time-stepping-methods/
    N = 2^8;
    x = {linspace(2*pi/N, 2*pi, N)'};
    domain = PSDomain(x, nan, false, 1);
    y0 = domain.fft(cos(domain.x{1}));
    params = struct('nu', 1);
    
    odefun = @(t, y) fkdv(domain, y, params);
    
    [t, y] = ab1(odefun, 0:1e-6:70e-6, y0);
    
    plot(t)
    figure; plot(log10(abs(y(1:end/2,1:10:end))));
    figure; plot(log10(abs(y(1:end/2,end))));
    figure; plot(domain.ifft(y(:,end)));
    %     verifyEqual(testCase, t(end), 1)
end

function testSemiImplicitBenney(testCase)
    addpath discretisationMethods/
    addpath ../time-stepping-methods/
    N = 2^6;
    x = {linspace(2*pi/N, 2*pi, N)'; linspace(2*pi/N, 2*pi, N)'};
    domain = FDDomain(x, [1,0;2,0;3,0;4,0;0,1;0,2;0,3;0,4;2,1;1,2;2,2]', ...
        2);
    %     domain = PSDomain(x, nan, false, 1);
    y0 = icos(domain.x);
    t = linspace(0, 0.1, 20);
    params = struct("theta", 7*pi/8, "Re", 1, "C", 0.01);
    
    y0 = domain.reshapeToVector(y0);
    odefunLinear = @(t) -fbenney2dLinear(domain, params);
    odefunNonlinear = @(t, y) domain.reshapeToVector( ...
        fbenney2d(domain, domain.reshapeToDomain(y), params)) - ...
        fbenney2dLinear(domain, params) * domain.reshapeToVector(y);
    
    options = struct("Linear", odefunLinear);
    
    y = ab3cn(odefunNonlinear, t, y0, options);
    y = domain.reshapeToDomain(y);
    
    
    save("temp5")
end

%% Create Data
function testCreateDataBenneyFiniteDifference(testCase)
    model = "benney";
    domain = createDomain(2*pi, 2*pi, 2^6, 2^6, "finite-difference");
    params = struct('theta', 1, 'Re', 1, 'C', 1);
    tFinal = 0.5;
    interface = @icos;
    method = "finite-difference";
    AbsTol = 1e-6;
    debug = false;
    
    [y, ~, ~] = createData(model, domain, params, tFinal, interface, method, AbsTol, debug);
    
    actual = y(:,:,end);
    load('testCreate2DBenneyEquationExpected','expected')
    
    verifyEqual(testCase, actual, expected, ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testCreateDataBenneyPseudoSpectral(testCase)
    model = "benney";
    domain = createDomain(2*pi, 2*pi, 96, 48, "pseudo-spectral");
    params = struct('theta', 1, 'Re', 1, 'C', 1);
    tFinal = 0.02;
    interface = @icos;
    method = "pseudo-spectral";
    AbsTol = 1e-5;
    debug = true;
    
    [y, t, ~] = createData(model, domain, params, tFinal, interface, method, AbsTol, debug);
    
    figure; plot(t)
    figure; surf(log10(abs(domain.fft(y(:,:,end)))));
    
    save("temp5")
    %     actual = y(:,:,end);
    %     load('temp','expected')
    %
    %     verifyEqual(testCase, actual, expected, ...
    %         'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testCreateDataWIBL1FiniteDifference(testCase)
    model = "wibl1";
    domain = createDomain(2*pi, 2*pi, 2^5, 2^5, "finite-difference");
    params = struct('theta', 1, 'Re', 1, 'C', 1);
    tFinal = 0.5;
    interface = @icos;
    method = "finite-difference";
    AbsTol = 1e-6;
    debug = false;
    
    [y, ~, ~] = createData(model, domain, params, tFinal, interface, method, AbsTol, debug);
    
    actual = y(:,:,end);
    load('testCreateWIBL1EquationExpected','expected')
    
    verifyEqual(testCase, actual, expected, ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

function testCreateDataWIBL1PseudoSpectral(testCase)
    model = "wibl1";
    domain = createDomain(2*pi, 2*pi, 2^5, 2^5, "pseudo-spectral");
    params = struct('theta', 1, 'Re', 1, 'C', 1);
    tFinal = 0.5;
    interface = @icos;
    method = "pseudo-spectral";
    AbsTol = 1e-6;
    debug = false;
    
    [y, ~, ~] = createData(model, domain, params, tFinal, interface, method, AbsTol, debug);
    
    actual = y(:,:,end);
    load('testCreateWIBL1EquationExpected','expected')
    
    verifyEqual(testCase, actual, expected, ...
        'RelTol', 2e-3, 'AbsTol', 1e-4)
end

%% Main

function testMainWIBL1(testCase)
    model = "benney";
    theta = 1; Re = 1; C = 1;
    xLength = 2*pi; yLength = 2*pi; xN = 2^6; yN = 2^6;
    tFinal = 0.5;
    interface = @icos;
    method = "finite-difference";
    AbsTol = 1e-6;
    debug = false;
    
    filename = makeFilename("", ...
        struct('theta', theta, 'Re', Re, 'C', C), ...
        setupX(xLength, yLength, xN, yN), ...
        tFinal,interface,AbsTol,model);
    if isfile(filename)
        delete(filename)
    end
    
    main(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method, debug)
    
    load(filename, 'y')
    if isfile(filename)
        delete(filename)
    end
    
    actual = y(:,:,end);
    load('testCreate2DBenneyEquationExpected','expected');
    
    verifyEqual(testCase, actual, expected, ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end
