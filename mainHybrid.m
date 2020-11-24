function main(theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, epsilon, delta, timeout)
    if nargin < 17
        timeout = -1;
    else
        graceTime = seconds(durationR2018('00:05:00'));
        timeout = seconds(durationR2018(timeout)) - graceTime;
    end
    addpath('timeSteppingMethods/')
    boolSemiImplicit = isSemiImplicit(timeStepper);

    params = struct('theta', theta, 'Re', Re, 'C', C, 'epsilon', epsilon, ...
        'delta', delta);

    domainArguments = struct('xLength', xLength, 'yLength', yLength, 'xN', xN, ...
        'yN', yN, 'method', method);

    if ~boolSemiImplicit
        odefun = @fhybrid;
        odejac = @jhybrid;
    else
        error("Not implemented");
    end
    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'odefun',odefun,'odejac',odejac,'interface',interface);

    timePointsArguments = struct('tStep', timeStepOut, 'tFinal', tFinal);


    myoptimoptions = optimoptions('fsolve', ...
        'Display', 'off');
    myoptimoptions.StepTolerance = RelTol * sqrt(xN * yN) * 2;

    if method == "finite-difference"
        myoptimoptions.SpecifyObjectiveGradient = true;
    end

    timerID = tic;
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        'RelTol', RelTol, ...
        'MaxStep', timeStep ...
        ...'InitialStep', 1e-3 ...
        ...'OutputFcn', 'odeprint',...
        ...'OutputSel', 1 ...
        );
    odeoptDefault.Events = @(t,y) wibl1Events(t, y, timerID, timeout);
    odeoptDefault.optimoptions = myoptimoptions;

    timeStepperArguments = struct('timeStepper', timeStepper, ...
        'odeopt', odeoptDefault);

    solution = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments);

    saveData(solution, ivpArguments, timePointsArguments, timeStepperArguments);

end
