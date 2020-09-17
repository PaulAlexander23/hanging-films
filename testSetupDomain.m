function tests = testSetupDomain()
    tests = functiontests(localfunctions);
end

function testSetupX(testCase)
    actual = setupX(10, pi, 12, 11);
    verifyEqual(testCase, actual{1}, linspace(10/12, 10, 12)')
    verifyEqual(testCase, actual{2}, linspace(pi/11, pi, 11))
end

function testFiniteDifference(testCase)
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "finite-difference");
    
    actual = setupDomain(args);
    expected = FDDomain(setupX(args.xLength, args.yLength, args.xN, args.yN),...
        [1, 0; 2, 0; 0, 1; 0, 2]', 4);
    
    verifyEqual(testCase, actual, expected);
end

function testPseudoSpectral(testCase)
    args = struct('xLength', 1, 'yLength', 2, ...
        'xN', 4, 'yN', 8, 'method', "pseudo-spectral");
    
    actual = setupDomain(args);
    expected = PSDomain(setupX(args.xLength, args.yLength, args.xN, args.yN),...
        true, true);
    
    verifyEqual(testCase, actual, expected);
end
