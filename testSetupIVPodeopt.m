function tests = testSetupIVPodeopt()
    tests = functiontests(localfunctions);
end

function testPseudoSpectral(testCase)
    domainArguments = struct('method', 'pseudo-spectral');
    args = struct('domainArguments', domainArguments);
    domain = 0;
    
    actual = setupIVPodeopt(args, domain);
    expected = odeset();
    
    verifyEqual(testCase, actual, expected);
end

function testFiniteDifferenceWIBL1(testCase)
    domainArguments = struct('method', 'finite-difference');
    args = struct('domainArguments', domainArguments, 'model', 'wibl1');
    domain = 0;
    
    actual = setupIVPodeopt(args, domain);
    expected = odeset();
    
    verifyEqual(testCase, actual, expected);
end

function testFiniteDifferenceBenney(testCase)
    domainArguments = struct('method', 'finite-difference');
    args = struct('domainArguments', domainArguments, 'model', 'benney');
    domain = 0;
    
    actual = setupIVPodeopt(args, domain);
    expected = odeset();
    
    verifyEqual(testCase, func2str(actual.Jacobian), '@(t,y)jbenney2d(domain,y,args.params)')
    actual.Jacobian = [];
    verifyEqual(testCase, actual, expected);
end