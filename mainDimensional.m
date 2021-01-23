function mainDimensional(model, theta, rho, nu, sigma, Re, xL, yL, tF, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, timeout)

    % theta should be 0 < theta < pi for these calculations
    theta = theta - pi;

    g = 9.81;
    h_N = (3 * nu^2 * Re / (g * sin(theta)))^(1/3);
    mu = nu * rho;
    u_N = h_N^2 * g * sin(theta) / (2 * nu);

    % Scale time after OR also scale timeStepOut and timeStep
    %t_S = h_N/u_N;
    %tFinal = tF/t_S;
    ReHang = h_N * u_N / nu;
    C = u_N * mu / sigma;
    xLength = xL/h_N;
    yLength = yL/h_N;

    % theta should be pi < theta < 2*pi for the main call
    theta = theta + pi;

    main(model, theta, ReHang, C, xLength, yLength, tF, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, timeout)

end
