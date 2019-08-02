function [y, t] = pdeSolver(pdefun, t, domain, y0, timestepper)
    shape = size(y0);
    y0 = reshape(y0, [prod(shape), 1]);
    
    [t, y] = timestepper( ...
        @(t, y) pdefunVector(pdefun, t, domain, y, shape), ...
        t, y0);
    
    y = y.';
    y = squeeze(reshape(y, [shape, length(t)]));
end

function F = pdefunVector(pdefun, t, domain, y, shape)
    y = reshapeToShape(y, shape);
    
    f = pdefun(t, domain, y);
    
    F = reshapeToVector(f, shape);
end

function Y = reshapeToShape(y, shape)
    Y = reshape(y, [shape, numel(y) / prod(shape)]);
end

function y = reshapeToVector(Y, shape)
    y = reshape(Y, prod(shape), []);
end