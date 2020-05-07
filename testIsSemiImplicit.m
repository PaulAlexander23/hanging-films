function tests = testIsSemiImplicit()
    tests = functiontests(localfunctions);
end

function testNonSemiImplicit(testCase)
    method = "ab1";
    expected = false;
    actual = isSemiImplicit(method);
    verifyEqual(testCase, actual, expected);
end

function testbdf1si(testCase)
    method = "bdf1si";
    expected = true;
    actual = isSemiImplicit(method);
    verifyEqual(testCase, actual, expected);
end

function testbdf2si(testCase)
    method = "bdf2si";
    expected = true;
    actual = isSemiImplicit(method);
    verifyEqual(testCase, actual, expected);
end
