function main(model, theta, rho, nu, sigma, Re, xL, yL, tF, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, timeout)
    if nargin < 16
        timeout = -1;
    else
        graceTime = seconds(durationR2018('00:05:00'));
        timeout = seconds(durationR2018(timeout)) - graceTime;
    end
    addpath('timeSteppingMethods/')
    boolSemiImplicit = isSemiImplicit(timeStepper);


    g = 9.81;
    h_N = (3 * nu^2 * Re / (g * sin(theta)))^(1/3);
    mu = nu * rho;
    u_N = h_N^2 * g * sin(theta) / (2 * nu);
    t_S = h_N/u_N;
    tFinal = tF/t_S;
    ReHang = h_N * u_N / nu;
    C = u_N * mu / sigma;
    xLength = xL/h_N;
    yLength = yL/h_N;
    

    params = struct('theta', theta, 'Re', ReHang, 'C', C);

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
        'RelTol', RelTol, ...
        'MaxStep', timeStep ...
        ...'InitialStep', 1e-3 ...
        ...'OutputFcn', 'odeprint',...
        ...'OutputSel', 1 ...
        );
    if model == "benney"
        odeoptDefault.Events = @(t,y) benneyEvents(t, y, timerID, timeout);
    elseif model == "wibl1"
        odeoptDefault.Events = @(t,y) wibl1Events(t, y, timerID, timeout);
    end
    odeoptDefault.optimoptions = myoptimoptions;

    timeStepperArguments = struct('timeStepper', timeStepper, ...
        'odeopt', odeoptDefault);

    solution = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments);

    saveData(solution, ivpArguments, timePointsArguments, timeStepperArguments);

end
