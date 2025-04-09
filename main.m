function main(odefun, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, timeout)

    % Check last argument is provided otherwise default to infinite timeout
    if nargin < 16
        timeout = -1;
    else
        graceTime = seconds(durationR2018('00:05:00'));
        timeout = seconds(durationR2018(timeout)) - graceTime;
    end

    % Include packages
    addpath('timeSteppingMethods/')

    params = struct('theta', theta, 'Re', Re, 'C', C);

    domainArguments = struct('xLength', xLength, 'yLength', yLength, 'xN', xN, ...
        'yN', yN, 'method', method);

    % Add jacobian if finite-difference method is being used, default to 1.
    odejac = @(domain, y, params) 1;
    jacobianSet = false;
    if method == "finite-difference"
        if func2str(odefun) == "fbenney2d"
            odejac = @jbenney2d;
            jacobianSet = true;
        elseif func2str(odefun) == "fwibl1STF"
            odejac = @jwibl1STF;
            jacobianSet = true;
        end
    end

    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'odefun',odefun,'odejac',odejac,'interface',interface);

    timePointsArguments = struct('tStep', timeStepOut, 'tFinal', tFinal);

    myoptimoptions = optimoptions('fsolve', ...
        'Display', 'off');

    % Use jacobian if set.
    if jacobianSet
        myoptimoptions.SpecifyObjectiveGradient = true;
    end

    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'RelTol', RelTol, ...
        'MaxStep', timeStep, ...
        ...'InitialStep', 1e-3 ...
        'OutputFcn', 'odeprint',...
        'OutputSel', 1 ...
        );

    timerID = tic;
    if func2str(odefun) == "fbenney2d"
        odeoptDefault.Events = @(t,y) benneyEvents(t, y, timerID, timeout);
        myoptimoptions.StepTolerance = RelTol * sqrt(xN * yN);
    elseif func2str(odefun) == "fwibl1STF"
        odeoptDefault.Events = @(t,y) wibl1Events(t, y, timerID, timeout);
        myoptimoptions.StepTolerance = RelTol * sqrt(xN * yN) * 2;
    elseif func2str(odefun) == "fwibl1" || func2str(odefun) == "fwibl2STF" || func2str(odefun) == "fwibl2"
        odeoptDefault.Events = @(t,y) wibl2Events(t, y, timerID, timeout);
        myoptimoptions.StepTolerance = RelTol * sqrt(xN * yN) * 3;
    end
    odeoptDefault.optimoptions = myoptimoptions;

    timeStepperArguments = struct('timeStepper', timeStepper, ...
        'odeopt', odeoptDefault);

    % Solve initial value problem
    solution = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments);

    saveData(solution, ivpArguments, timePointsArguments, timeStepperArguments);

end
