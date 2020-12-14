function tests =  testEvents()
    tests = functiontests(localfunctions);
end

function testTimeoutTriggered(testCase)

    timerID = tic;
    options = odeset('Events', @(t,y) eTimeout(timerID, 1e-1), 'MaxStep', 1e-6);

    [~,~,te,~,~] = ode45(@(t,y)-y, [0,1], 1, options);

    verifyEqual(testCase, length(te), 1)
end

function testTimeoutNotTriggered(testCase)

    timerID = tic;
    options = odeset('Events', @(t,y) eTimeout(timerID, 1e0));

    [~,~,te,~,~] = ode45(@(t,y)1, [0,1], 1, options);

    verifyEqual(testCase, length(te), 0)
end

function testNanTriggered(testCase)

    options = odeset('Events', @(t,y) eNan(y));

    [~,~,te,~,~] = ode45(@odefun, [0,1], 0.5, options);

    verifyEqual(testCase, length(te), 1)

    function f = odefun(t,y)
        if t > 0.1
            f = nan;
        else
            f = 1;
        end
    end
end

function testNanNotTriggered(testCase)

    options = odeset('Events', @(t,y) eNan(y));

    [~,~,te,~,~] = ode45(@(t,y)1, [0,1], 1, options);

    verifyEqual(testCase, length(te), 0)
end

function testDewetTriggered(testCase)

    options = odeset('Events', @(t,y) eDewet(y));

    [~,~,te,~,~] = ode45(@(t,y)-1, [0,1], 0.5, options);

    verifyEqual(testCase, length(te), 1)
end

function testDewetNotTriggered(testCase)

    options = odeset('Events', @(t,y) eDewet(y));

    [~,~,te,~,~] = ode45(@(t,y)-1, [0,1], 1.5, options);

    verifyEqual(testCase, length(te), 0)
end

function testDewetWIBL1Triggered(testCase)

    options = odeset('Events', @(t,y) eDewetWIBL1(y));

    [~,~,te,~,~] = ode45(@(t,y)[-1;-1], [0,1], [0.5;0.4], options);

    verifyEqual(testCase, length(te), 1)
end

function testDewetWIBL1NotTriggered(testCase)

    options = odeset('Events', @(t,y) eDewetWIBL1(y));

    [~,~,te,~,~] = ode45(@(t,y)[-1;-1], [0,1], [1.5;0.4], options);

    verifyEqual(testCase, length(te), 0)
end

function testDewetWIBL2Triggered(testCase)

    options = odeset('Events', @(t,y) eDewetWIBL2(y));

    [~,~,te,~,~] = ode45(@(t,y)[-1;-1;-1], [0,1], [0.5;0.4;1], options);

    verifyEqual(testCase, length(te), 1)
end

function testDewetWIBL2NotTriggered(testCase)

    options = odeset('Events', @(t,y) eDewetWIBL2(y));

    [~,~,te,~,~] = ode45(@(t,y)[-1;-1;1], [0,1], [1.5;0.4;0.1], options);

    verifyEqual(testCase, length(te), 0)
end
