function tests = testCreateData()
    tests = functiontests(localfunctions);
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
    load('data/testCreate2DBenneyEquationExpected','expected')
    
    verifyEqual(testCase, actual, expected, ...
        'RelTol', 1e-3, 'AbsTol', 1e-6)
end

% function testCreateDataBenneyPseudoSpectral(testCase)
%     model = "benney";
%     domain = createDomain(2*pi, 2*pi, 96, 48, "pseudo-spectral");
%     params = struct('theta', 1, 'Re', 1, 'C', 1);
%     tFinal = 0.02;
%     interface = @icos;
%     method = "pseudo-spectral";
%     AbsTol = 1e-5;
%     debug = true;
%     
%     [y, t, ~] = createData(model, domain, params, tFinal, interface, method, AbsTol, debug);
%     
%     figure; plot(t)
%     figure; surf(log10(abs(domain.fft(y(:,:,end)))));
%     
%     save("temp5")
%     %     actual = y(:,:,end);
%     %     load('temp','expected')
%     %
%     %     verifyEqual(testCase, actual, expected, ...
%     %         'RelTol', 1e-3, 'AbsTol', 1e-6)
% end

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
    load('data/testCreateWIBL1EquationExpected','expected')
    
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
    load('data/testCreateWIBL1EquationExpected','expected')
    
    verifyEqual(testCase, actual, expected, ...
        'RelTol', 1e-1, 'AbsTol', 1e-1)
end