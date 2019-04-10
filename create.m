function create(theta,Re,We,C,x_length,y_length,t_final)
    
    dim = 2;
    xL = [x_length,y_length];
    xN = [32,64];
    xS = xL./xN;
    x = cell(dim,1);
    for n = 1:dim
        x{n} = linspace(xS(n),xL(n),xN(n))';
    end
    
    t = [0,t_final];
    
    A = 2e-1;
    r = 0.05;
    y0 = 1 + A * (-r*cos(2*pi/xL(1) * x{1}) - cos(2*pi/xL(2) * x{2}'));
    
    params = [1,theta,Re,We,C]; % delta, theta, Re, We, C
    problemDeg = [1,0;0,1;2,0;0,2]';
    
    D = init_fd(x, problemDeg, 4);
    method = @(x,y,deg) diff_fd(x,y,deg,D,problemDeg);
    getD = @(deg) get_fd(deg,D,problemDeg);
    pdefun = @(t,x,y,method) benney(x,y,params,method,getD);
    
    odeopt = odeset( ...
        'Jacobian', @(t, y) jbenney(x, y, params, method, getD), ...
        'Event', @event_dewetted, ...
        'Vectorized','on' ...
        ...'BDF','on' ... % Backward differentiation formulae
        );
    timestepper = @(odefun,t,y0) ode15s(odefun,t,y0,odeopt);
    
    tic
    [y, t] = solver(pdefun, t, x, y0, method, timestepper);
    timeTaken = toc
    
    filename = replace(sprintf('data-theta-%g-Re-%g-We-%g-C-%g-xL-%g-yL-%g-T-%g',[theta,Re,We,C,x_length,y_length,t_final]),'.','_');
    save(filename,'y','params','t','x','timeTaken');
    
    
end
