function tests = testSetupIVPodeopt()
    tests = functiontests(localfunctions);
end

function testPseudoSpectral(testCase)
    args = struct("method", "pseudo-spectral");
    domain = 0;
    
    actual = setupIVPodeopt(args, domain);
    expected = odeset();
    
    verifyEqual(testCase, actual, expected);
end

function testFiniteDifferenceWIBL1(testCase)
    args = struct("method", "finite-difference", "model", "wibl1");
    domain = 0;
    
    actual = setupIVPodeopt(args, domain);
    expected = odeset();
    
    verifyEqual(testCase, actual, expected);
end

function testFiniteDifferenceBenney(testCase)
    args = struct("method", "finite-difference", "model", "benney");
    domain = 0;
    
    actual = setupIVPodeopt(args, domain);
    expected = odeset();
    
    verifyEqual(testCase, func2str(actual.Jacobian), '@(t,y)jbenney2d(domain,y,args.params)')
    actual.Jacobian = [];
    verifyEqual(testCase, actual, expected);
end