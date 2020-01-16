function tests = testSetupTimePoints()
    tests = functiontests(localfunctions);
end

function testNoRemainder(testCase)
    timePointsArguments = struct('tStep', 0.2, 'tFinal', 10);
    
    actual = setupTimePoints(timePointsArguments);
    
    verifyEqual(testCase, actual(end), 10);
    verifyTrue(testCase, all(abs(diff(actual(1:end-1))-0.2) < 1e-14));
end

function testRemainder(testCase)
    timePointsArguments = struct('tStep', 0.3, 'tFinal', 10);
    
    actual = setupTimePoints(timePointsArguments);
    
    verifyEqual(testCase, actual(end), 10);
    verifyTrue(testCase, all(abs(diff(actual(1:end-1))-0.3) < 1e-14));
end