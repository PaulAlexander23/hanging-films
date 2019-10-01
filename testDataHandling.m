function tests = testDataHandling()
    tests = functiontests(localfunctions);
end

function testMakeFilename(testCase)
    params = struct('theta', 2, 'Re', 3, 'C', 4);
    actual = makeFilename('-test', params, {0.1:0.1:3.3, 1:1:10}, pi, @interface, 1e-6, "benney");
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
    expectedT = linspace(0, 0.126576313737181);
    expectedX = {0.1:0.1:3.3, 1:1:10};
    expectedTimeTaken = 0.925425280986515;
    saveData(expectedY, expectedParams, expectedT, expectedX, ...
        expectedTimeTaken, pi, @interface, 1e-6, "benney", "-test");

    verifyTrue(testCase, isfile(filename))

    [actualY, actualParams, actualT, actualX, actualTimeTaken] = loadData(expectedParams, expectedX, pi, @interface, 1e-6, "benney", "-test");

    verifyEqual(testCase, actualY, expectedY)
    verifyEqual(testCase, actualParams, expectedParams)
    verifyEqual(testCase, actualT, expectedT)
    verifyEqual(testCase, actualX, expectedX)
    verifyEqual(testCase, actualTimeTaken, expectedTimeTaken)

    if isfile(filename)
        delete(filename)
    end
end
