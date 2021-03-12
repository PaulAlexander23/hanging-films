function main(tFinal, delta, timeout)
    addpath('~/Repositories/hanging-films/');
    addpath('~/Repositories/hanging-films/timeSteppingMethods/');
    addpath('~/Repositories/hanging-films/discretisationMethods/');

    if nargin < 2
        timeout = -1;
    else
        graceTime = seconds(durationR2018('00:05:00'));
        timeout = seconds(durationR2018(timeout)) - graceTime;
    end

    interface = @(x) iload(x,'ic-hybrid-Re-1-delta-0_75.mat');

    domainArguments = struct('xLength', 32, 'yLength', 32, 'xN', 64, ...
        'yN', 64, 'method', 'finite-difference');

    odefun = @fhybrid;
    odejac = @jhybrid;

    params = struct('Re', 1, 'C', 0.01, 'theta', 7/8*pi, 'epsilon', 0.25, 'delta', delta);

    ivpArguments = struct('domainArguments',domainArguments,'params',params,...
        'odefun',odefun,'odejac',odejac,'interface',interface);

    timeStepOut = 0.2;

    timePointsArguments = struct('tStep', timeStepOut, 'tFinal', tFinal);

    timerID = tic;
    odeoptDefault = odeset( ...
        ...'Vectorized', 'on', ...
        ...'BDF','on', ...
        ...'RelTol', RelTol, ...
        ...'MaxStep', timeStep ...
        ...'InitialStep', 1e-3 ...
        ...'OutputFcn', 'odeprint',...
        ...'OutputSel', 1 ...
        );

    odeoptDefault.Events = @(t,y) wibl1Events(t, y, timerID, timeout);

    timeStepperArguments = struct('timeStepper', @ode15s, ...
        'odeopt', odeoptDefault);

    solution = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments);

    saveData(solution, ivpArguments, timePointsArguments, timeStepperArguments);
end
