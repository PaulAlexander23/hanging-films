function [y, x, t] = solver(F, params, ic, tFinal, xL, xN, evnt, RelTol)
    %SOLVER Computes the numerical solution up to tFinal
    %   Detailed explanation goes here
    
    if nargin < 12
        RelTol = 1e-3;
    end
    
    dim = length(xL);
    xS = xL./xN;
    x = cell(dim,1);
    for k = 1:dim
        x{k} = linspace(xS(k), xL(k), xN(k))';
    end
    
    shape = size(ic);
    ic = reshape(ic, [prod(shape),1]);
    
    func = @(t,y) reshape(- F(reshape(y,shape), xL, params), [prod(shape),1]);
    
    function [value, isterminal, direction] = event(t,y)
        Y = reshape(y,shape);
        value = 2*evnt(t,Y) - 1;
        isterminal = 1; % Terminal
        direction = 0; % Any approach direction
    end
    
    options = odeset(...
        ...'Vectorized','on',...
        'BDF','on',... % Backward differentiation formulae
        'Event', @(t,y) event(t,y),...
        ...'OutputFcn','odeprint',...
        'RelTol',RelTol); % Default: 1e-3
    
    [t, y] = ode15s(func, [0,tFinal], ic, options);
    y = y';
    
    y = reshape(y,[shape,length(t)]);
    
end
