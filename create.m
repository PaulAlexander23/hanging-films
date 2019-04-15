function create(theta,Re,C,x_length,y_length,t_final,interface,xN,yN,AbsTol)
    
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, AbsTol = 1e-6; end
    
    dim = 2;
    xL = [x_length,y_length];
    xN = [xN,yN];
    xS = xL./xN;
    x = cell(dim,1);
    for n = 1:dim
        x{n} = linspace(xS(n),xL(n),xN(n))';
    end
    
    t = [0,t_final];
    
    y0 = interface(x);
    
    params = [1,theta,Re,C]; % delta, theta, Re, C
    problemDeg = [1,0;0,1;2,0;0,2]';
    
    D = init_fd(x, problemDeg, 4);
    method = @(x,y,deg) diff_fd(x,y,deg,D,problemDeg);
    getD = @(deg) get_fd(deg,D,problemDeg);
    pdefun = @(t,x,y,method) fbenney(x,y,params,method);
    
    odeopt = odeset( ...
        'Jacobian', @(t, y) jbenney(x, y, params, method, getD), ...
        'Event', @event_dewetted, ...
        'Vectorized', 'on', ...
        'AbsTol', AbsTol ...
        ...'BDF','on' ... % Backward differentiation formulae
        );
    timestepper = @(odefun,t,y0) ode15s(odefun,t,y0,odeopt);
    
    tic
    [y, t] = solver(pdefun, t, x, y0, method, timestepper);
    timeTaken = toc;
    
    filename = replace(sprintf('data-theta-%g-Re-%g-C-%g-xL-%g-yL-%g-T-%g-interface-%s-xN-%g-yN-%g-AbsTol-%g', ...
        theta,Re,C,x_length,y_length,t_final,func2str(interface),xN(1),xN(2),AbsTol),'.','_');
    save(filename,'y','params','t','x','timeTaken');
    
end
