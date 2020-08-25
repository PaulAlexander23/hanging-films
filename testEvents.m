function tests = testEvents()
    tests = functiontests(localfunctions);
end

function testDewet(testCase)
    y = [1,2;3,4]; 

    [value, isTerminal, direction] = eDewet(y);

    verifyEqual(testCase, value, 1);
    verifyEqual(testCase, isTerminal, 1);
    verifyEqual(testCase, direction, 0);
end

function testDewetWIBL1(testCase)
    y = [1,2;0,4]; 

    [value, isTerminal, direction] = eDewetWIBL1(y);

    verifyEqual(testCase, value, 1);
    verifyEqual(testCase, isTerminal, 1);
    verifyEqual(testCase, direction, 0);
end

