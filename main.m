function main(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, timeout)
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, RelTol = 1e-6; end
    if nargin < 12, method = 'finite-difference'; end
    if nargin < 13, timeStepper = @ode15s; end
    if nargin < 14, timeStepOut = 0.2; end
    if nargin < 15, timeStep = timeStepOut; end
    if nargin < 16
        timeout = -1;
    else
        graceTime = seconds(durationR2018('00:05:00'));
        timeout = seconds(durationR2018(timeout)) - graceTime;
    end
    addpath('timeSteppingMethods/')
    boolSemiImplicit = isSemiImplicit(timeStepper);

    params = struct('theta', theta, 'Re', Re, 'C', C);

    domainArguments = struct('xLength', xLength, 'yLength', yLength, 'xN', xN, ...
        'yN', yN, 'method', method);

    if ~boolSemiImplicit
        if model == "benney"
            odefun = @fbenney2d;
            odejac = @jbenney2d;
        elseif model == "wibl1"
            odefun = @fwibl1;
            odejac = @jwibl1;
        end
    else
        if model == "benney"
            explicitOdefun = @fbenney2dExplicit;
            implicitOdefun = @fbenney2dImplicit;
            odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
            odejac = @jbenney2dImplicit;
        elseif model == "wibl1"
            explicitOdefun = @fwibl1Explicit;
            implicitOdefun = @fwibl1Implicit;
            odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
            odejac = @jwibl1Implicit;
        end
    end
    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'odefun',odefun,'odejac',odejac,'interface',interface);

    timePointsArguments = struct('tStep', timeStepOut, 'tFinal', tFinal);


    myoptimoptions = optimoptions('fsolve', ...
        'Display', 'off');
    if model == "benney"
        myoptimoptions.StepTolerance = RelTol * sqrt(xN * yN);
    elseif model == "wibl1"
        myoptimoptions.StepTolerance = RelTol * sqrt(xN * yN) * 2;
    end

    if method == "finite-difference"
        myoptimoptions.SpecifyObjectiveGradient = true;
    end

    timerID = tic;
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'Events', @eventNan,...
        'RelTol', RelTol, ...
        'Events', @(~, ~) eTimeout(timerID, timeout), ...
        'MaxStep', timeStep ...
        ...'InitialStep', 1e-3 ...
        ...'OutputFcn', 'odeprint',...
        ...'OutputSel', 1 ...
        );
    odeoptDefault.optimoptions = myoptimoptions;

    timeStepperArguments = struct('timeStepper', timeStepper, ...
        'odeopt', odeoptDefault);

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
