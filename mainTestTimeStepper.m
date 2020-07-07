function mainTestTimeStepper(theta, Re, C, xL, yL, tFinal, xN, yN, timeStepper, timeStep, tStepOut, timeout)
    if nargin < 12
        timeout = -1;
    else
        graceTime = seconds(durationR2018('00:05:00'));
        timeout = seconds(durationR2018(timeout)) - graceTime;
    end
    addpath('timeSteppingMethods/')
    addpath('discretisationMethods');

    t = (0:tStepOut:tFinal)';
    x = {linspace(xL/xN,xL,xN),linspace(yL/yN,yL,yN)};
    domain = FDDomain(x, [1, 0; 2, 0; 0, 1; 0, 2]', 4);
    params = struct('theta', theta, 'Re', Re, 'C', C, ...
        'a', 0.5, 'b', 0.25, 'c', 2);

    h0 = 1 - params.a * cos(2*pi*domain.x{1} / xL) ...
        - params.b * cos(2 * pi * domain.x{2} / yL);

    y0 = domain.reshapeToVector(h0);

    explicitOdefun = @(t, y) fbenney2dExplicit(domain, y, params);
    implicitOdefun = @(t, y) fbenney2dImplicit(domain, y, params) ...
          + fbenneyforcing(t, domain, params);
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    odejac = @(t, y) jbenney2dImplicit(domain, y, params);

    timerID = tic;
    myoptimoptions = optimoptions('fsolve', ...
        'Display', 'off', ...
        'SpecifyObjectiveGradient', true);
    options = odeset( ...
        'Events', @(t,y) myEvents(t, y, timerID, timeout), ...
        'MaxStep', timeStep);
    options.optimmethod = @fsolve;
    options.optimoptions = myoptimoptions;
    options.Jacobian = odejac;

    tic
    solution = timeStepper(odefun, t, y0, options);
    timeTaken = toc;

    solution.y = reshape(solution.y', [size(solution.y,2), 1, size(solution.y,1)]);
    solution.y = domain.reshapeToDomain(solution.y);

    t = reshape(t, [1, 1, length(t)]);
    expected = 1 - params.a * cos(2*pi*domain.x{1}/xL - params.c * t) ...
        - params.b * cos(2*pi*domain.x{2}/yL);
    actual = solution.y;
    err = abs(actual - expected);

    save("data.mat")

    function [value, isterminal, direction] = myEvents(t, y, timerID, timeout)
        [value1, isterminal1, direction1] = eventNan(t, y);
        [value2, isterminal2, direction2] = eTimeout(timerID, timeout);

        value = [value1; value2];
        isterminal = [isterminal1; isterminal2];
        direction = [direction1; direction2];
    end

    function [value, isterminal, direction] = eventNan(~, y)
        value = double(~any(isnan(y)));
        isterminal = 1;
        direction = 0;
    end

    function secondsOut = durationR2018(stringIn)
        secondsOut = duration(cell2mat(cellfun(@str2num,split(char(stringIn),':'),'UniformOutput',false)'));
    end
end
