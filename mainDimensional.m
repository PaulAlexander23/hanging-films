function mainDimensional(model, theta, rho, nu, sigma, Re, xL, yL, tF, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, timeout)
    if nargin < 16
        timeout = -1;
    else
        graceTime = seconds(durationR2018('00:05:00'));
        timeout = seconds(durationR2018(timeout)) - graceTime;
    end

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

    main(model, theta, ReHang, C, xLength, yLength, tFinal, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, timeout)

end
