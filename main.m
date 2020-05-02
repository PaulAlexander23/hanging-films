function main(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method, timeStepper, timeout)
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, AbsTol = 1e-6; end
    if nargin < 12, method = "finite-difference"; end
    if nargin < 13, timeStepper = @ode15s; end
    if nargin < 14
        timeout = -1;
    else
        graceTime = seconds(durationR2018('00:05:00'));
        timeout = seconds(durationR2018(timeout)) - graceTime;
    end

    params = struct('theta', theta, 'Re', Re, 'C', C);

    domainArguments = struct('xLength', xLength, 'yLength', yLength, 'xN', xN, ...
        'yN', yN, 'method', method);

    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'model',model,'interface',interface);

    timePointsArguments = struct('tStep', 0.2, 'tFinal', tFinal);

    timerID = tic;
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'Events', @eventNan,...
        'AbsTol', AbsTol, ...
        'Events', @(~, ~) eTimeout(timerID, timeout) ...
        ...'MaxStep', 5e-6 ...
        ...'InitialStep', 1e-3 ...
        );
    odeoptOutput = odeset( ...
        'OutputFcn', 'odeprint',...
        'OutputSel', 1 ...
        );

    timeStepperArguments = struct('timeStepper', timeStepper, ...
        'odeopt', odeoptDefault, 'outputOpt', odeoptOutput);

    solution = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments);

    saveData(solution, ivpArguments, timePointsArguments, timeStepperArguments);

    function [value, isterminal, direction] = eventNan(~, y)
        value = double(~any(isnan(y)));
        isterminal = 1;
        direction = 0;
    end

    function secondsOut = durationR2018(stringIn)
        secondsOut = duration(cell2mat(cellfun(@str2num,split(char(stringIn),':'),'UniformOutput',false)'));
    end
end
