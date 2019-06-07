function create(theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol)

    if nargin < 8, xN = 64; end
    if nargin < 9, yN = 64; end
    if nargin < 10, AbsTol = 1e-6; end

    x = setupX(xLength, yLength, xN, yN);

    t = [0, tFinal];

    y0 = interface(x);

    params = [1, theta, Re, C]; % delta, theta, Re, C
    problemDiffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';

    diffMat = init_fd(x, problemDiffDegrees, 4);
    diffMethod = @(x, y, degree) diff_fd(x, y, degree, diffMat, problemDiffDegrees);
    getDiffMat = @(deg) get_fd(deg, diffMat, problemDiffDegrees);
    pdeFunction = @(t, x, y, diffMethod) fbenney(x, y, params, diffMethod);

    odeopt = odeset( ...
        'Jacobian', @(t, y) jbenney(x, y, params, diffMethod, getDiffMat), ...
        'Event', @event_dewetted, ...
        'Vectorized', 'on', ...
        'AbsTol', AbsTol ...
        ... 'BDF','on' ...
        );
    timeStepper = @(odefun, t, y0) ode15s(odefun, t, y0, odeopt);

    tic
    [y, t] = solver(pdeFunction, t, x, y0, diffMethod, timeStepper);
    timeTaken = toc;

    saveData(y, params, t, x, timeTaken, tFinal, interface, AbsTol)
end

function x = setupX(xLength, yLength, xN, yN)
    xL = [xLength, yLength];
    xN = [xN, yN];
    xS = xL ./ xN;
    dim = 2;
    x = cell(dim, 1);
    for n = 1:dim
        x{n} = linspace(xS(n), xL(n), xN(n))';
    end
end

function saveData(y, params, t, x, timeTaken, tFinal, interface, AbsTol)
    filename = replace(sprintf('data-theta-%g-Re-%g-C-%g-xL-%g-yL-%g-T-%g-interface-%s-xN-%g-yN-%g-AbsTol-%g', ...
        params(2), params(3), params(4), x{1}(end), x{2}(end), tFinal, func2str(interface), length(x{1}), length(x{2}), AbsTol), '.', '_');
    save(filename, 'y', 'params', 't', 'x', 'timeTaken');
end