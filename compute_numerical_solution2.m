function [y, x, t] = compute_numerical_solution2(F, params, ic, tFinal, xL, xN, RelTol)
    %COMPUTE_NUMERICAL_SOLUTION Computes the numerical solution up to tFinal
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
    
    func = @(t,y) reshape(- F(reshape(y,shape), xL, params) , [prod(shape),1]);
    
    options = odeset(...
        ...'Vectorized','on',...
        'BDF','on',... % Backward differentiation formulae
        ...'Event',@(t,y) event_collision(t,y,H1,H2,xS),...
        ...'OutputFcn','odeprint',...
        'RelTol',RelTol,... % Default: 1e-3
        'AbsTol',1e-6);  % Default: 1e-6

    [t, y] = ode15s(func, [0,tFinal], ic, options);
    y = y';
    
    y = reshape(y,[shape,length(t)]);

end