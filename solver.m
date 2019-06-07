function [y, t] = solver(pdefun, t, x, y0, method, timestepper)
    %SOLVER Computes the numerical solution up to tFinal
    %   Detailed explanation goes here

    shape = size(y0);
    y0 = reshape(y0, [prod(shape), 1]);

    y = timestepper(@(t, y) callPDEFunctionWithJacobianSort(pdefun, ...
        t, x, y, method, shape), t, y0);

    if isstruct(y)
        t = y.x;
        y = y.y;
    end
    
    y = squeeze(reshape(y, [shape, length(t)]));
end

function [F, J] = callPDEFunctionWithJacobianSort(pdefun, t, x, y, method, shape)
    if nargout == 1
        f = pdefun(t, x, reshapeInput(y, shape), method);
    elseif nargout == 2
        [f, j] = pdefun(t, x, reshapeInput(y, shape), method);
        J = j;
    end
    F = reshapeOutput(y, f, shape);
end

function out = reshapeInput(y, shape)
     out = reshape(y, [shape, numel(y) / prod(shape)]);
end

function F = reshapeOutput(y, f, shape)
    F = reshape(f, [prod(shape), numel(y) / prod(shape)]);
end