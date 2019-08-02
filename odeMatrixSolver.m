function [y, t] = odeMatrixSolver(odefun, t, y0, timestepper)
    shape = size(y0);
    y0 = reshape(y0, [prod(shape), 1]);
    
    [t, y] = timestepper( ...
        @(t, y) pdefunVector(odefun, t, y, shape), ...
        t, y0);
    
    y = y.';
    y = squeeze(reshape(y, [shape, length(t)]));
end

function F = pdefunVector(odefun, t, y, shape)
    y = reshapeToShape(y, shape);
    
    f = odefun(t, y);
    
    F = reshapeToVector(f, shape);
end

function Y = reshapeToShape(y, shape)
    Y = reshape(y, [shape, numel(y) / prod(shape)]);
end

function y = reshapeToVector(Y, shape)
    y = reshape(Y, prod(shape), []);
end