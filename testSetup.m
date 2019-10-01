function tests = testSetup()
    tests = functiontests(localfunctions);
end

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
